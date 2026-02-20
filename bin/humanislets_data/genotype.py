import pandas as pd
import sys, os
import warnings
import math
import numpy as np
warnings.filterwarnings("ignore", category=FutureWarning)

# Description: Takes a CSV file created by genotypes extracted using samtools mpileup
# Filters, genotypes, and converts the files to VCF

### Parameters ###

read_depth = 10
eps = 0.01 # error-rate
EMPTY = '.'

# read CSV file

header_file = 'empty_header.txt'
file = sys.argv[1]
outdir=sys.argv[2]
# file = os.path.basename(file)
SAMPLE = os.path.basename(file).strip('.csv')

# directory = 'eQTL_variants_extracted2_96samples'
# file = os.path.join(directory, file)
if not os.path.exists(outdir):
    os.mkdir(outdir)

df = pd.read_csv(file)

# filter CSV

dff = df.loc[df.sum(axis=1) >= read_depth].copy()

dff['REF'] = dff['variant'].apply(lambda x:x.split('_')[1])
dff['ALT'] = dff['variant'].apply(lambda x:x.split('_')[2])
dff = dff.set_index('variant')

print(dff.head())

def _log_binom_pmf(k, n, p):
    if p == 0.0:
        return 0.0 if k == 0 else -math.inf
    if p == 1.0:
        return 0.0 if k == n else -math.inf
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + k * math.log(p)
        + (n - k) * math.log(1 - p)
    )

def genotype(ref_reads, alt_reads, eps=eps, min_depth=read_depth):
    """
    binomial-likelihood genotype caller.

    Returns: '0/0', '0/1', '1/1', or np.nan
    """
    n = ref_reads + alt_reads
    if n < min_depth:
        return np.nan

    k = alt_reads

    loglik = {
        "0/0": _log_binom_pmf(k, n, eps),
        "0/1": _log_binom_pmf(k, n, 0.5),
        "1/1": _log_binom_pmf(k, n, 1.0 - eps),
    }

    return n, alt_reads/n, max(loglik, key=loglik.get)

def get_dbSNP(CHROM, POS):
    pass
    return EMPTY

def VCF_line(var, df):
    REF = df.loc[var,'REF']
    ALT = df.loc[var,'ALT']
    if len(REF) == 1 and len(ALT) == 1:
        CHROM = var.split(':')[0]
        POS = var.split(':')[1].split('_')[0]
        Gene = var.split('_')[3]
        ID = get_dbSNP(CHROM, POS)
        QUAL = EMPTY
        FILTER = 'PASS'
        REF_counts = df.loc[var,REF]
        ALT_counts = df.loc[var,ALT]
        DP, AF, GT = genotype(REF_counts, ALT_counts)

        INFO = f'DP={DP};AF={AF};Gene={Gene}'
        FORMAT = f'GT:DP:AF'
        SAMPLE_FIELD = f'{GT}:{DP}:{AF}'
        END = '\n'

        line = '\t'.join([CHROM,
                          POS,
                          ID,
                          REF,
                          ALT,
                          QUAL,
                          FILTER,
                          INFO,
                          FORMAT,
                          SAMPLE_FIELD])+END

        return line

def write_VCF(df, outfile, verbose=False):
    with open(outfile, 'w') as vcf:
        with open(header_file, 'r') as header:
            for line in header.readlines():
                vcf.write(line.replace('<SAMPLE>', SAMPLE))
            vcf.write('\n')
            variants = list(df.index)
            for var in variants:
                try:
                    line = VCF_line(var, df)
                    if verbose:
                        print(line)
                    vcf.write(line)
                except:
                    continue

outfile_tmp = os.path.join(outdir, f'{SAMPLE}.vcf.tmp')
outfile = os.path.join(outdir, f'{SAMPLE}.vcf')
write_VCF(dff, outfile_tmp)
os.system(f'bcftools sort {outfile_tmp} > {outfile}')
os.system(f'bgzip -f {outfile}')
os.system(f'bcftools index {outfile}.gz')
os.remove(outfile_tmp)
print(f'{outfile}.gz made')