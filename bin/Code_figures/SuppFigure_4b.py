#!/bin/python

###########################
##  eQTL SUPP Figure 4B  ##
###########################

import matplotlib.pyplot as plt
import scanpy as sc
import glob
import pandas as pd
import decoupler as dc

filenameout1 = 'supp_figure4b_1.png'
filenameout2 = 'supp_figure4b_2.png'

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

# vcf files
vcf_files = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/2-variants/1-VCF/smartseqPancreas_3UTR_gene_annotated_chr*.vcf'

VCFfiles = glob.glob(vcf_files)
chromosomes = [i.rsplit("chr")[1].rsplit(".vcf")[0] for i in VCFfiles]
vcf_dict = {}
for i,chr in enumerate(chromosomes):
    vcf_dict[chr] = VCFfiles[i]

def get_matrix(adata):
    matrix = pd.DataFrame(adata.raw.X, columns=adata.raw.var.index, index=adata.obs.index)
    return matrix

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

colors_cells = {'acinar-cells': '#1f77b4', 'ductal-cells':'#ff7f0e', 'endothelial-cells':'#279e68', 'PP-cells':'#d62728','stellate-cells': '#aa40fc', 'alpha-cells':'#8c564b', 'beta-cells':'#e377c2', 'delta-cells':'#b5bd61'}

def get_adata_sub(adata, geno, celltype, group):
    adata_sub = adata[adata.obs['Donor'].isin(list(geno.keys()))].copy()
    adata_sub = adata_sub[adata_sub.obs['Cell_type'] == celltype].copy()
    adata_sub = adata_sub[adata_sub.obs['Group'] == group].copy()
    adata_sub.obs['genotype'] = adata_sub.obs['Donor'].apply(lambda x: geno[x][1])
    adata_sub.obs['alleles'] = adata_sub.obs['Donor'].apply(lambda x: geno[x][0])
    
    return adata_sub

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
adata_sub = get_adata_sub(adata, geno, 'beta-cells', 'healthy')
acts_pat_sub = get_adata_sub(acts_pat, geno, 'beta-cells', 'healthy')

fig, ax = plt.subplots(dpi=300)
x = pd.DataFrame(acts_pat_sub.obsm['mlm_estimate'])['EGFR']
y = get_matrix(adata_sub)['PTEN']

ax.scatter(x,y, edgecolor='black')
ax.set_xlabel('EGFR activity')
ax.set_ylabel('PTEN')
ax.axvline(0, color='black')
ax.set_title('ND beta-cells')
ax.grid(False)

fig.savefig(filenameout1, bbox_inches='tight')

adata_sub = get_adata_sub(adata, geno, 'beta-cells', 'T2D')
acts_pat_sub = get_adata_sub(acts_pat, geno, 'beta-cells', 'T2D')

fig, ax = plt.subplots(dpi=300)
x = pd.DataFrame(acts_pat_sub.obsm['mlm_estimate'])['EGFR']
y = get_matrix(adata_sub)['PTEN']

ax.scatter(x,y, edgecolor='black')
ax.set_xlabel('EGFR activity')
ax.set_ylabel('PTEN')
ax.axvline(0, color='black')
ax.set_title('T2D beta-cells')
ax.grid(False)

fig.savefig(filenameout2, bbox_inches='tight')