import pandas as pd
import numpy as np
import vcf
import mygene
import time
import matplotlib.pyplot as plt
from statsmodels.stats import multitest
import random
import os, sys

# description: this script runs limma correlations for cell type specific eQTL variants using humanislets.com data.
# Allele dosage is compared to proteomic (cis-effects) and GSIS data

### GSIS ###

directory = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/GSIS_geno_v2'
GSIS_file = 'gsis.csv'
GSIS_file = os.path.join(directory, GSIS_file)

df_gsis = pd.read_csv(GSIS_file, index_col=0)

df_gsis = df_gsis.loc[:,['replicate','culturetime2', 'gluc_conc_group', 'stim_index', 'first_insulin_percent', 'second_insulin_percent']]

gsis_data = {}
donors = []

def outlier_removal(data_list):
    data = np.asarray(data_list, dtype=float)

    # mask invalid values
    valid = data > 0

    log10 = np.full_like(data, np.nan)
    log10[valid] = np.log10(data[valid])

    SD = np.nanstd(log10)
    mean = np.nanmean(log10)

    cutoff = mean + 4 * SD

    out = data.copy()
    out[log10 > cutoff] = np.nan

    return out

donors = list(set(df_gsis.index))
df1 = pd.DataFrame(index=donors)
for gluc_conc_group, dfg in df_gsis.groupby('gluc_conc_group'):
    gsis_data = []
    df1[gluc_conc_group] = np.zeros(len(donors))
    df1[f'{gluc_conc_group}_first_perc'] = np.zeros(len(donors))
    df1[f'{gluc_conc_group}_second_perc'] = np.zeros(len(donors))
    df1[gluc_conc_group] = np.zeros(len(donors))
    donor_order = []
    for donor, dfgd in dfg.groupby('record_id'):
        donor_order.append(donor)
        median = np.median(dfgd['stim_index'])
        gsis_data.append(median)
    gsis_data = outlier_removal(gsis_data)
    for donor, si in zip(donor_order, gsis_data):
        df1.loc[donor,gluc_conc_group] = si

    donor_order = []
    data_1st = []
    data_2nd = []
    for donor, dfgd in dfg.groupby('record_id'):
        donor_order.append(donor)
        median_1st = np.median(dfgd['first_insulin_percent'])
        median_2nd = np.median(dfgd['second_insulin_percent'])
        data_1st.append(median_1st)
        data_2nd.append(median_2nd)
    data_1st = outlier_removal(data_1st)
    data_2nd = outlier_removal(data_2nd)
    for donor, perc1, perc2 in zip(donor_order, data_1st, data_2nd):
        df1.loc[donor,f'{gluc_conc_group}_first_perc'] = perc1
        df1.loc[donor,f'{gluc_conc_group}_second_perc'] = perc2


### genotypes from Humanislet.com dataset ###
def get_geno_df(eQTLs):
    directory = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/GSIS_geno'
    vcf_file = 'eQTL_variants_extracted2_96samples.vcf.gz'
    vcf_file = os.path.join(directory, vcf_file)
    VCF = vcf.Reader(filename=vcf_file)

    # parse VCF file to extract genotypes
    genotypes_by_sample = {sample: [] for sample in VCF.samples}
    variant_ids = []
    for record in VCF:
        ref = record.REF
        alt = str(record.ALT).strip('[').strip(']')
        chrom = record.CHROM
        pos   = record.POS
        gene  = record.INFO.get("Gene", "NA")
        variant_id = f"{chrom}:{pos}_{ref}_{alt}_{gene}"
        if variant_id in eQTLs:
            variant_ids.append(variant_id)
            for sample in record.samples:
                gt = sample['GT']
                if gt.count('.') > 0:
                    dosage = np.nan
                else:
                    dosage = gt.count('1')
                genotypes_by_sample[sample.sample].append(dosage)

    df2 = pd.DataFrame(index = genotypes_by_sample.keys(), columns=variant_ids)
    for donor in genotypes_by_sample.keys():
        for variant, dosage in zip(variant_ids, genotypes_by_sample[donor]):
            df2.loc[donor, variant] = dosage
    return df2

### eQTL variants ###
def get_eQTL_df(celltype, group):
    directory = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/eQTL_v6'
    eQTL_file = 'merged/all-cells_allgroups.adjusted-pvalue.csv'
    eQTL_file = os.path.join(directory, eQTL_file)

    df_all = pd.read_csv(eQTL_file, index_col=0)
    df_sig = df_all[(df_all['adj_pvalue'] <= 0.05)&(df_all['Cell_group'].isin([f'{celltype}_{group}']))]

    return df_sig

### Proteomics ###

def get_proteomics_df(eQTLs):
    # read files
    directory = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/GSIS_geno/proteomic_data_humanislets'
    prot_file = 'proc_prot.csv'
    prot_file = os.path.join(directory, prot_file)

    df_prot = pd.read_csv(prot_file, index_col=0).T
    df_prot.index.names = ['record_id']
    df_prot.columns = [str(idx) for idx in df_prot.columns]
    genes = [variant.split('_')[-1] for variant in eQTLs]

    def get_conversion_dict(symbols):
        mg = mygene.MyGeneInfo()
        results = mg.querymany(symbols, scopes="symbol", fields="entrezgene", species="human")

        symbol2id = {}
        id2symbol = {}
        for r in results:
            if "notfound" not in r and "entrezgene" in r:
                symbol = r["query"]
                entrez = str(r["entrezgene"])
                symbol2id[symbol] = entrez
                id2symbol[entrez] = symbol
        return id2symbol, symbol2id

    id2symbol, symbol2id = get_conversion_dict(set(genes))

    protein_ids = list(set(symbol2id.values()).intersection(set(df_prot.columns)))
    df3 = df_prot.loc[:,protein_ids]
    df3.columns = [id2symbol.get(c) for c in df3.columns]

    return df3

### Perfusions ###

# Data from human islets, normalized by me using dynamic_gsis_norm.py

directory = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/GSIS_geno/'
summary_file = 'peri_gluc.norm_TdW.summary.csv'#'function_summary.csv'
summary_file = os.path.join(directory, summary_file)
df_perf = pd.read_csv(summary_file, index_col=0)
df7 = df_perf
for col in df7.columns:
    df7[col] = outlier_removal(list(df7[col]))

### Meta data ###
directory = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/GSIS_geno/meta_data/'
donor_file = 'donor.csv'
donor_file = os.path.join(directory, donor_file)
isolation_file = 'isolation.csv'
isolation_file = os.path.join(directory, isolation_file)
cellpro_file = 'cell_pro.csv'
cellpro_file = os.path.join(directory, cellpro_file)

df_donor = pd.read_csv(donor_file, index_col=0)
df_iso = pd.read_csv(isolation_file, index_col=0)
df_cellpro = pd.read_csv(cellpro_file, index_col=0)

df4 = df_donor
df5 = df_iso
df6 = df_cellpro

celltypes = ['alpha-cells', 'beta-cells']
groups = ['ND']


df1.columns = ['-'.join(['GSIS', c]) for c in df1.columns]
df4.columns = ['-'.join(['meta', c]) for c in df4.columns]
df5.columns = ['-'.join(['meta', c]) for c in df5.columns]
df6.columns = ['-'.join(['meta', c]) for c in df6.columns]
df7.columns = ['-'.join(['GSIS', c]) for c in df7.columns]

for celltype in celltypes:
    for group in groups:
        df_sig = get_eQTL_df(celltype, group)
        eQTLs = list(set(df_sig['Variant']))
        df2 = get_geno_df(eQTLs)
        df3 = get_proteomics_df(eQTLs)

        df2.columns = ['-'.join(['eQTL', c]) for c in df2.columns]
        df3.columns = ['-'.join(['Proteinexp', c]) for c in df3.columns]

        df_merge = pd.concat([df1, df2, df3, df4, df5, df6, df7], axis=1)

        print(df_merge.columns)
        print(df_merge.shape)

        if group == 'ND':
            df_merge = df_merge[(df_merge['meta-T2_diabetes'] == 0)&(df_merge['meta-T1_diabetes'] == 0)]
        elif group == 'T2D':
            df_merge = df_merge[(df_merge['meta-T2_diabetes'] == 1)&(df_merge['meta-T1_diabetes'] == 0)]

        # print(df_merge.shape)

        df_merge.to_csv(f'humanislets_data_{celltype}.{group}.csv')

        os.system(f'Rscript regression.R humanislets_data_{celltype}.{group}.csv')

        # FDR correction
        df_res = pd.read_csv(f'humanislets_data_{celltype}.{group}.regression.csv', index_col=0)
        columns = [col for col in df_res.columns if col != 'adj_P_within']
        df_res_prot = df_res.loc['Proteomics', columns]
        rejected, pvals_adjusted = multitest.fdrcorrection(list(df_res_prot['p_value']), alpha=0.05, method='indep', is_sorted=False)
        df_res_prot['adj_p-value'] = pvals_adjusted
        df_res_prot.to_csv(f'humanislets_data_{celltype}.{group}.regression.proteomics.csv')

        df_res_gsis = df_res[df_res['response'] == 'GSIS-1_16.7'].loc[:,columns]
        df_res_gsis = df_res_gsis.dropna(subset='p_value')
        rejected, pvals_adjusted = multitest.fdrcorrection(list(df_res_gsis['p_value']), alpha=0.05, method='indep', is_sorted=False)
        df_res_gsis['adj_p-value'] = pvals_adjusted
        df_res_gsis.to_csv(f'humanislets_data_{celltype}.{group}.regression.GSIS.csv')