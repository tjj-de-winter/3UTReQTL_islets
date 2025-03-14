#!/bin/python

###########################
##    eQTL Figure 5D     ##
###########################

import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import glob
import decoupler as dc
import statsmodels.stats.multitest as sm
import scipy.stats as stats

filenameout = 'figure5d.png'

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

def difference(acts, groupby='alleles', groups=['T/T', 'C/C']):
    group1 = groups[0]
    group2 = groups[1]
    
    acts1 = pd.DataFrame(acts[acts.obs[groupby] == group1].obsm['mlm_estimate'])
    acts2 = pd.DataFrame(acts[acts.obs[groupby] == group2].obsm['mlm_estimate'])
    
    
    
    index = acts1.columns
    df = pd.DataFrame(index=index, columns=[f'mean_{group1}', f'mean_{group2}',f'max_{group1}',f'max_{group2}', 'pval', 'adj_pval'])
    
    pvals = []
    for gene in index:
        data1, data2 = acts1[gene], acts2[gene]
        pval = stats.ttest_ind(data1, data2)
        pvals.append(pval[1])
        
    pvals_adj = sm.multipletests(pvals, method='bonferroni')
    
    df['pval'] = pvals
    df['adj_pval'] = pvals_adj[1] 
    df[f'mean_{group1}'] = list(acts1.loc[:,index].mean(axis=0))
    df[f'mean_{group2}'] = list(acts2.loc[:,index].mean(axis=0))
    df[f'max_{group1}'] = list(acts2.loc[:,index].max(axis=0))
    df[f'max_{group2}'] = list(acts2.loc[:,index].max(axis=0))
    df['fold_difference'] = df[f'mean_{group1}']/df[f'mean_{group2}']
    
    return df

progeny = dc.get_progeny(organism='human', top=500)

dc.run_mlm(
    mat=adata,
    net=progeny,
    source='source',
    target='target',
    weight='weight',
    verbose=True
)

acts_pat = dc.get_acts(adata, obsm_key='mlm_estimate')

geno = get_genotypes(adata, vcf_dict['10'], '10:87966988', group='all')
acts_pat_sub = get_adata_sub(acts_pat, geno, 'beta-cells', 'healthy')

df_diff = difference(acts_pat_sub, groupby='alleles', groups=['T/T', 'C/C'])
df_active = df_diff[(df_diff['max_T/T'] > 0.5) & (df_diff['max_C/C'] > 0.5)]
df_active = df_active.sort_values(by='adj_pval', ascending=True)

data = -np.log10(df_active['adj_pval'])
labels = df_active['adj_pval'].index

fig, ax = plt.subplots(dpi=300)

x = range(len(data))

ax.barh(x, data, color='#2102bd')
ax.axvline(-np.log10(0.05), ls='--', color='black')
ax.set_xlabel('adjusted -log10(p-value)')
ax.set_ylabel('signaling pathway activity')
ax.grid(False)
ax.set_title('T/T beta-cells vs C/C beta-cells')

ax.set_yticks(x, labels);

fig.savefig(filenameout, bbox_inches='tight')