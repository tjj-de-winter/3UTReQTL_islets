#!/bin/python

###########################
##    eQTL Figure 5G     ##
###########################

import scanpy as sc
import glob
import matplotlib.pyplot as plt
from vcf import Reader
import numpy as np
from collections import Counter
import os

filenameout = 'figure5g.png'

sc.set_figure_params(scanpy=True, fontsize=30, dpi=300, facecolor='white')
sc.settings.figdir = './'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

sc.set_figure_params(scanpy=True, fontsize=8, dpi=300, facecolor='white')

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

geno = get_genotypes(adata, vcf_dict['10'], '10:87966988', group='all')
adata_sub = get_adata_sub(adata, geno, 'beta-cells', 'healthy')

sc.pl.dotplot(adata_sub, ['FOXO3'], groupby='alleles', show=False, save=filenameout)
# plt.savefig(filenameout, dpi=300, bbox_inches='tight')

os.system('mv dotplot_figure5g.png figure5g.png')
