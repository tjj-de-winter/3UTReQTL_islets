#!/bin/python

###########################
##    eQTL Figure 1E     ##
###########################

import matplotlib.pyplot as plt
import pandas as pd
import glob
import scanpy as sc
import random
from vcf import Reader
import numpy as np
from collections import Counter

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv'
df = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df = df[df['pvalue'] < 0.05]
df['gene'] = df['gene'].apply(lambda x: x.rsplit('>')[0])

df = df.drop_duplicates()

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

def get_matrix(adata):
    matrix = pd.DataFrame(adata.raw.X, columns=adata.raw.var.index, index=adata.obs.index)
    return matrix

def get_data(genotypes, adata, gene):
    
    matrix = get_matrix(adata)
    
    data = []
    gt = sorted(set([genotypes[cell][1] for cell in genotypes]))
    for g in gt:
        cell_index = adata.obs[adata.obs['genotype'] == g].index
        matrix_list = list(matrix.loc[cell_index, gene])
        data.append(matrix_list)
        
    return data

def onecelltype_boxplot(data1, genotypes, color_rgb, title='rs0001', ylabel='gene X', fontsize=10):     
    ''' data1 = [REF, HET, ALT]
    genotypes = ['T/T', 'T/C', 'C/C']'''
    
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['figure.titlesize'] = fontsize
    plt.rcParams['font.size'] = fontsize

    fig, ax = plt.subplots(dpi=300, figsize=(2,3))

    x = range(len(genotypes))
    delta = 0
    position1 = [xi-delta for xi in x]

    def individual_points(data, position):
        xs = []
        ys = []
        for d, x in zip(data, position):
            scatter = [random.uniform(x-0.07, x+0.07) for i in range(len(d))]
            xs += scatter
            ys += list(d)
        return xs,ys
    
    def rgb_to_rgba(rgb, alpha=1.0):
        return tuple([channel / 255.0 for channel in rgb] + [alpha])
    
    x,y = individual_points(data1, position1)
    ax.scatter(x, y, s=1, color='black', zorder=10)

    ax.boxplot(data1, 
               positions = position1, 
               showcaps=False, 
               patch_artist=True, 
               boxprops=dict(facecolor=rgb_to_rgba(color_rgb, alpha=1.0), color='black'),
               medianprops=dict(color='black'),
               showfliers=False
              );

    ax.set_xticks(range(len(genotypes)))
    ax.set_xticklabels(genotypes, rotation=45)
    ax.set_title(title.replace('healthy', 'ND'))
    ax.set_ylabel(ylabel)
    ax.grid(False)
    
    return fig, ax
    
def hex_to_rgb(hexa):
    return tuple(int(hexa[i:i+2], 16)  for i in (0, 2, 4))
    
def plot(gene, chromosome, variant, celltype, group, color, set_ylim=False, ylim=(0,6), nameout='test', fontsize=10):
    genotypes = get_genotypes(adata, vcf_dict[chromosome], variant, group='all')
    adata_sub = get_adata_sub(adata, genotypes, celltype, group)
    alleles = {genotypes[cell][1]:genotypes[cell][0] for cell in genotypes}    

    data = get_data(genotypes, adata_sub, gene)
    
    alleles_list = [alleles[geno] for geno in sorted(alleles)]
        
    
    fig, ax = onecelltype_boxplot(data, alleles_list, color, title=f'{variant}\n{celltype} {group}', ylabel=f'{gene} expression', fontsize=fontsize)

    if set_ylim:
        ax.set_ylim(ylim)

    fig.savefig(nameout, bbox_inches='tight')

colors_cells = {'acinar-cells': '#1f77b4', 'ductal-cells':'#ff7f0e', 'endothelial-cells':'#279e68', 'PP-cells':'#d62728','stellate-cells': '#aa40fc', 'alpha-cells':'#8c564b', 'beta-cells':'#e377c2', 'delta-cells':'#b5bd61'}

plot_data = {
    '1': ['SLC30A8', '8:117172786', 'beta-cells', 'healthy'],
    '2': ['SLC30A8', '8:117173494', 'beta-cells', 'healthy'],
}

for idx in plot_data:
    gene = plot_data[idx][0]
    variant = plot_data[idx][1]
    chromosome = variant.rsplit(':')[0]
    celltype = plot_data[idx][2]
    group = plot_data[idx][3]
    color = colors_cells[celltype].strip('#')
    
    print(gene, chromosome, variant, celltype, group)

    plot(gene, chromosome, variant, celltype, group, hex_to_rgb(color), nameout=f'figure1e_{idx}.png',fontsize=14)