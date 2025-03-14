#!/bin/python

###########################
##    eQTL Figure 5B     ##
###########################

import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import random

filenameout = 'figure5b.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

# vcf files
vcf_files = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/2-variants/1-VCF/smartseqPancreas_3UTR_gene_annotated_chr*.vcf'

VCFfiles = glob.glob(vcf_files)
chromosomes = [i.rsplit("chr")[1].rsplit(".vcf")[0] for i in VCFfiles]
vcf_dict = {}
for i,chr in enumerate(chromosomes):
    vcf_dict[chr] = VCFfiles[i]

def get_genotypes(adata, vcf, variant, group='all'):
    """Extract the allele sequences for a specific genetic variant from a 
    variant call format file and make a boxplot from this variant """
    from vcf import Reader
    import scanpy
    import numpy as np
    from collections import Counter
    
    # import VCF file
    VCF = Reader(fsock=None, filename=vcf, prepend_chr=False, strict_whitespace=False, encoding='ascii')
    
    # filter adata for group
    
    if group == 'all':
        adatax = adata
    else:
        adatax = adata[(adata.obs['Group'] == group)]
        
    donor_list = list(set(adatax.obs['Donor'])) # list of donors to select which genotype to use in plot
   
    # read VCF file
    genotype_dict = {}
    for ii, record in enumerate(VCF):
        rec = str(record)
        varname = str(rec.rsplit("CHROM=")[1].rsplit(",")[0] + ":"
                  + rec.rsplit("POS=")[1].rsplit(",")[0]) # varname = chr:position
        if varname in variant:
            ref = str(rec.rsplit("REF=")[1].rsplit(",")[0]) # get the reference sequence
            alts = str(rec.rsplit('ALT=[')[1].rsplit("]")[0]) # get the alternative sequences
            alts = alts.rsplit(',') 
            allels = [ref] + [i.strip() for i in alts] # append the ref and alt sequences in 1 list
            for sample in record.samples: # loop over the genotype results for given variant
                s = str(sample)
                donor = str(s.rsplit("sample=")[1].rsplit(",")[0]) # get donor name
                #print(donor)
                if donor in donor_list:
                    GT = str(s.rsplit("GT=")[1].rsplit(",")[0]) # e.g. GT=1/1
                    GT = GT.replace('|','/')
                    try: # skip if genotype contains '.' (e.g. not called)
                        al1 = int(GT.rsplit('/')[0])
                        al2 = int(GT.rsplit('/')[1])
                    except:
                        continue

                    if al1 == 0 and al2 == 0:
                        gtype = '0-Homozygous REF'
                    elif al1 == 0 and al2 > 0:
                        gtype = '1-Heterozygous'
                    else:
                        gtype = '2-Homozygous ALT'

                    genotype_dict[donor] = []
                    genotype_dict[donor].append(str(allels[al1] + "/" + allels[al2])) # append variant sequences to dict
                    genotype_dict[donor].append(str(gtype + '-' + str(al2)))# append donor genotype to dict
            break
    return genotype_dict

def get_adata_sub(adata, geno, celltype, group):
    adata_sub = adata[adata.obs['Donor'].isin(list(geno.keys()))].copy()
    adata_sub = adata_sub[adata_sub.obs['Cell_type'] == celltype].copy()
    adata_sub = adata_sub[adata_sub.obs['Group'] == group].copy()
    adata_sub.obs['genotype'] = adata_sub.obs['Donor'].apply(lambda x: geno[x][1])
    adata_sub.obs['alleles'] = adata_sub.obs['Donor'].apply(lambda x: geno[x][0])
    
    return adata_sub

# get genotypes
geno = get_genotypes(adata, vcf_dict['10'], '10:87966988', group='all')
adata_sub = get_adata_sub(adata, geno, 'beta-cells', 'healthy')

def merge(allele):
    if "C" in allele:
        return 'C'
    else:
        return 'T'
    
# adata_sub.obs['merged_alleles'] = adata_sub.obs['alleles'].apply(lambda a: merge(a))

adata_sub = adata_sub[adata_sub.obs['alleles'].isin(['T/T', 'C/C'])]

# differential gene expression

sc.tl.rank_genes_groups(adata_sub, 'alleles', method='wilcoxon')

result = adata_sub.uns['rank_genes_groups']

group = 'T/T'

df = pd.DataFrame({'logfoldchanges': result['logfoldchanges'][group],
              'pvals_adj': result['pvals_adj'][group]}, index=result['names'][group])

fig, ax = plt.subplots(dpi=300, figsize=(10,5))

df_updown = df.copy()

dfx = df_updown[df_updown['logfoldchanges'] < 0].copy()

x = list(dfx['logfoldchanges'])
y = [-np.log10(p) for p in dfx['pvals_adj']]

ax.scatter(x,y, edgecolors='black', color='#B5C228')

dfx = df_updown[df_updown['logfoldchanges'] > 0].copy()

x = list(dfx['logfoldchanges'])
y = [-np.log10(p) for p in dfx['pvals_adj']]

ax.scatter(x,y, edgecolors='black', color='#9C3C8E')

ax.set_title('beta-cells')
ax.set_xlabel('log2(foldchanges)')
ax.set_ylabel('-log10(p-value)')

ax.axhline(-np.log10(0.05), color='black', linestyle='--')
ax.axvline(0, color='black', linestyle='--')

df_text_up = df[df['pvals_adj'] < 0.05].sort_values(by='logfoldchanges', ascending=True).head(6)
df_text_down = df[df['pvals_adj'] < 0.05].sort_values(by='logfoldchanges', ascending=False).head(6)

def annotate(x,y,s, xlim=(-25,25), ymax=3):
    if x < 0:
        if x < xlim[0]:
            xi = x + random.uniform(3, 4)
        else:
            xi = x - random.uniform(8, 9)
    else:
        if x > xlim[1]:
            xi = x - random.uniform(8, 9)
        else:
            xi = x + random.uniform(3, 4)
            
    yi = ymax+1
    while yi > 3:
        yi = y + random.uniform(-0.25, 0.25)
        
    ax.annotate(s, xy=(x,y), xytext=(xi,yi), arrowprops=dict(facecolor='black', arrowstyle='->'), size=10)
    


df_text = pd.concat([df_text_up, df_text_down])
for s in df_text.index:
    x = df_text.loc[s,'logfoldchanges']
    y = -np.log10(df_text.loc[s,'pvals_adj'])
    annotate(x,y,s)

ax.grid(False)

fig.savefig(filenameout, bbox_inches='tight')