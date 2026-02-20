#!/bin/python

# compute zero‑inflated negative‑binomial (ZINB) eQTLs using an R package (via rpy2)
# run gene batches in parallel with pandarallel

import multiprocessing as mp
mp.set_start_method("spawn", force=True)

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os, sys, gc
from pathlib import Path
from pandarallel import pandarallel
import argparse

# rpy2 / R
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
pandas2ri.activate()
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

# Import R packages
pscl  = importr('pscl')
stats = importr('stats')
base  = importr('base')

### Input variables ###

parser = argparse.ArgumentParser(description='Compute eQTLs for a given cell type using variant and normalized count matrices')
parser.add_argument('--celltype', required=True, help='Cell type name, e.g. beta-cells_ND')
parser.add_argument('--covariates', help='CSV file with covariates; defaults to <celltype>.covariates.csv')
parser.add_argument('--variants', help='CSV variant matrix with variants as rows, samples as cols; defaults to <celltype>_variant_matrix.csv')
parser.add_argument('--normcounts', help='DESeq normalized gene count matrix; defaults to <celltype>_filtered_count_matrix_normalized.csv')
parser.add_argument('--chromosome', default=None, help='Chromosome to analyze (e.g. chr10). Pass None to scan all chromosomes')
parser.add_argument('--outdir', help='Output directory; defaults to ../<celltype>_eQTL')
parser.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4)')

args = parser.parse_args()

if args.covariates is None:
    args.covariates = f'{args.celltype}.covariates.csv'
if args.variants is None:
    args.variants = f'{args.celltype}_variant_matrix.csv'
if args.normcounts is None:
    args.normcounts = f'{args.celltype}_filtered_count_matrix_normalized.csv'
if args.outdir is None:
    args.outdir = f'../{args.celltype}_eQTL'

celltype = args.celltype
covariates_file = args.covariates
variant_file = args.variants
normcounts_file = args.normcounts
chromosome = args.chromosome
threads = args.threads
outdir = args.outdir

os.system(f'mkdir -p {outdir}')

### Parameters ###
n_PCs_gene     = 5
remove_outlier = True
EM             = True
dist           = 'negbin'
flush_every    = 15  # write results every N rows

### Functions ###

# Remove high non‑zero outliers (same function as the SCeQTL R package)
ro.r('''
removeoutlier <- function(sample.data){
  nz <- sample.data$expression[sample.data$expression != 0]
  med <- median(nz);  mad_ <- mad(nz)
  sample.data[sample.data$expression < med + 4 * mad_, ]
}
''')

def extract_genes(chrom):
    """return unique genes that have ≥1 variant on chrom"""
    if chrom is not None:
        vars_sel = [v for v in df_var_full.columns if v.split(':')[0] == chrom]
    else:
        vars_sel = list(df_var_full.columns)
    return list({v.split('_')[-1] for v in vars_sel})

def make_gene_chunks(genes, nchunks):
    return [list(chunk) for chunk in np.array_split(genes, nchunks)]

def convert_str_int(str_list):
    unique = sorted(list(set(str_list)))
    return  [int(unique.index(i)) for i in str_list]

### Read data ###
print('start reading data')
df_covar_full = pd.read_csv(covariates_file, index_col=0)
df_var_full   = pd.read_csv(variant_file,   index_col=0)
df_gene_full  = pd.read_csv(normcounts_file, index_col=0)

fixed     = ["Age", "Sex", "Dataset", "log_total_counts", "pct_counts_mt"]
genePCs   = [f"genePC{i}" for i in range(1, n_PCs_gene+1)]
cov_cols  = fixed + genePCs

df_covar_donor = df_covar_full.loc[:,["DonorID"]]
df_covar = df_covar_full.loc[:, cov_cols]

### eQTL computation using pscl.zeroinfl() ###

# R formulas
rhs = " + ".join(cov_cols)
f_full = Formula(f"expression ~ snp + {rhs} | {rhs}")
f_null = Formula(f"expression ~ {rhs} | {rhs}")

ctrl = pscl.zeroinfl_control(maxit = 300, trace = False)

def Run(run_idx, chunks):
    genes_chunk = chunks[run_idx]  # plain Python list

    # output file per chunk
    out_path = (f"{outdir}/{run_idx}.chr{chromosome}.{celltype}.eQTL.csv"
                if chromosome else f"{outdir}/{run_idx}.{celltype}.eQTL.csv")
    Path(out_path).unlink(missing_ok=True)
    header = ["Gene","Variant","pvalue","Beta","Rho_donor","Rho_cell",
          "med0_cell","med1_cell","med2_cell",
          "med0_donor","med1_donor","med2_donor"]
    pd.DataFrame(columns=header).to_csv(out_path, index=False)

    buffer = []
    variants = [v for v in df_var_full.columns if v.split('_')[-1] in genes_chunk]

    for var in variants:
        gene = var.split('_')[-1]
        if gene not in df_gene_full.index:
            continue

        # sample selection
        cells = df_var_full[df_var_full[var].isin((0, 1, 2))].index
        if len(cells) < 10:
            continue


        covars = df_covar.loc[cells, cov_cols]

        sex_levels = ['Female','Male'] if set(covars['Sex']) >= {'Female','Male'} else sorted(covars['Sex'].dropna().unique())
        dataset_levels = sorted(covars['Dataset'].dropna().unique())

        covars['Sex']     = pd.Categorical(covars['Sex'], categories=sex_levels)
        covars['Dataset'] = pd.Categorical(covars['Dataset'], categories=dataset_levels)

        covars.loc[:,'Sex'] = convert_str_int(list(covars['Sex']))
        covars.loc[:,'Dataset'] = convert_str_int(list(covars['Dataset']))
        expression = df_gene_full.loc[gene, cells].astype(float)
        snp        = df_var_full.loc[cells, var].astype(int)

        sample_df = pd.concat(
            [pd.DataFrame({'expression': expression, 'snp': snp}), covars], axis=1
        )

        # drop rows with NA / Inf (prevents glm NA/NaN/Inf errors)
        sample_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        sample_df.dropna(inplace=True)

        # require both genotype categories present and ≥10 samples
        if sample_df.shape[0] < 10 or sample_df['snp'].nunique() < 2:
            continue

        # optional outlier removal
        sample_r = (
            ro.r['removeoutlier'](pandas2ri.py2rpy(sample_df))
            if remove_outlier else pandas2ri.py2rpy(sample_df)
        )
        if int(base.nrow(sample_r)[0]) < 10:
            continue

        # fit full model
        try:
            m1 = pscl.zeroinfl(f_full, data=sample_r, dist=dist, EM=EM, control=ctrl)

            with localconverter(default_converter):          # NumPy/pandas disabled *inside* block
                coef_r = ro.r('coef')(m1)                    # rpy2 FloatVector, keeps .rx2
            beta = float(coef_r.rx2('count_snp')[0])         # NumPy/pandas disabled *inside* block
        except Exception:
            continue

        # fit null model and LRT
        pval_nan = False
        try:
            m0 = pscl.zeroinfl(f_null, data=sample_r, dist=dist, EM=EM, control=ctrl)
            ll1, ll0 = map(lambda m: float(stats.logLik(m)[0]), (m1, m0))
            pval = float(stats.pchisq(2 * (ll1 - ll0), df=1, lower_tail=False)[0])
        except Exception:
            pval = np.nan
            pval_nan = True

        if not pval_nan:
            donors = df_covar_donor.loc[cells, ['DonorID']]
            dfx = pd.concat([
                donors,
                expression.rename('expr'),
                snp.rename('dosage')
            ], axis=1)

            dfx = dfx[dfx['dosage'].isin([0,1,2])].copy()
            dfx['dosage'] = dfx['dosage'].astype(int)

            # cell-level Spearman
            rho_cell, _ = spearmanr(dfx['dosage'].to_numpy(), dfx['expr'].to_numpy())

            # donor-level medians within dosage
            g = (dfx.groupby(['dosage','DonorID'], as_index=False)
                    .agg(expr_med=('expr','median')))

            rho_donor, _ = spearmanr(g['dosage'].to_numpy(), g['expr_med'].to_numpy())

            # dosage-wise medians (fill missing with NaN)
            med_cell  = dfx.groupby('dosage')['expr'].median()
            med_donor = g.groupby('dosage')['expr_med'].median()

            med0_cell  = float(med_cell.get(0,  np.nan))
            med1_cell  = float(med_cell.get(1,  np.nan))
            med2_cell  = float(med_cell.get(2,  np.nan))
            med0_donor = float(med_donor.get(0, np.nan))
            med1_donor = float(med_donor.get(1, np.nan))
            med2_donor = float(med_donor.get(2, np.nan))

        else:
            rho_donor = np.nan
            rho_cell = np.nan
            med0_cell  = np.nan
            med1_cell  = np.nan
            med2_cell  = np.nan
            med0_donor  = np.nan
            med1_donor  = np.nan
            med2_donor  = np.nan

        buffer.append((gene, var, pval, beta, rho_donor, rho_cell, med0_cell, med1_cell, med2_cell, med0_donor, med1_donor, med2_donor))

        # flush buffer
        if len(buffer) >= flush_every:
            pd.DataFrame(buffer, columns=header).to_csv(out_path, mode='a', header=False, index=False)
            buffer.clear(); gc.collect()

    # final flush
    if buffer:
        pd.DataFrame(buffer, columns=header).to_csv(out_path, mode='a', header=False, index=False)

    return 'DONE'

### Chunk preparation & parallel run ###

genes   = extract_genes(chromosome)
chunks  = make_gene_chunks(genes, len(genes))
submit_df = pd.DataFrame({'chunk': range(len(chunks))})

pandarallel.initialize(progress_bar=True, nb_workers=threads)
submit_df['run'] = submit_df['chunk'].parallel_apply(lambda i: Run(i, chunks))
