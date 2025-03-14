#!/bin/python

###########################
##  eQTL SUPP Figure 1J  ##
###########################

import PyComplexHeatmap as pch
import matplotlib.pyplot as plt
import glob
import pandas as pd
import json
import numpy as np

filenameout = 'supp_figure1j.png'

files = glob.glob('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/5-public_data/GSEA_genesets/*.json')

gene_sets = {}
for file in files:
    name = file.rsplit('/')[-1].rsplit('_')[2].lower()
    if name == 'mesenchymal':
        name = 'stellate'
    if name == 'gamma':
        name = 'PP'
    f = open(file)
    dataset = json.load(f)
    for i in dataset:
        gene_set = dataset[i]['geneSymbols']
        
    gene_sets[name] = gene_set

gwas_files =  glob.glob('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/5-public_data/GWAS/GWAS_SNPS/*.tsv')

gwas_dict = {}
for file in gwas_files:
    df = pd.read_csv(file, sep="\t")
    
    name = file.rsplit('/')[-1].rsplit('.EFO')[0].rsplit('.MONDO')[0].replace('Trait.', '').replace('.',' ').capitalize()
    df['Trait'] = [name for idx in df.index]
    
    gwas_dict[file.rsplit('.tsv')[0]] = df
gwas_dict

df_gwas = pd.concat([gwas_dict[key] for key in gwas_dict])
df_gwas.index = range(df_gwas.shape[0])

file="/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv"

df_eqtl_raw = pd.read_csv(file, sep=",", names=['gene', 'variant', 'pval', 'group'])
df_eqtl_raw['locations'] = df_eqtl_raw['variant'].apply(lambda v: v.rsplit('_')[0])
df_eqtl = pd.read_csv(file, sep=",", names=['gene', 'variant', 'pval', 'group'])
df_eqtl = df_eqtl.drop_duplicates()
df_eqtl = df_eqtl[df_eqtl['pval'] < 0.05]
df_eqtl['gene'] = df_eqtl['gene'].apply(lambda g: g.rsplit('>')[0])
df_eqtl['locations'] = df_eqtl['variant'].apply(lambda v: v.rsplit('_')[0])
df_eqtl['ND_T2D'] = df_eqtl['group'].apply(lambda x: x.rsplit('_')[1])


GWAS_SNPs = set(list(df_gwas['locations']))
eQTL_SNPs = set(list(df_eqtl['locations']))

shared_SNPs = list(GWAS_SNPs.intersection(eQTL_SNPs))

df_eqtl_selection = df_eqtl[df_eqtl['locations'].isin(shared_SNPs)]
df_eqtl_selection = df_eqtl_selection[df_eqtl_selection['ND_T2D'] == 'healthy']  

def get_average_pval_gwas(df, location):
    gwas_dict = {}
    
    df_sub = df[df['locations'] == location]
    for trait, dfi in df_sub.groupby('traitName'):
        pval = dfi.mean()['pValue']
        pval = -np.log10(pval)
        gwas_dict[trait] = pval
                
    return gwas_dict

def order_variants(variant_list):
    '''Order variants based on chromosome number and nucleotide position'''
    variant_dict = {v.rsplit('_')[0].replace('X','23').replace('Y','24'): v for v in variant_list}
    
    chromosomes = set(int(v.rsplit(':')[0]) for v in variant_dict.keys())
    
    variants_neworder = []
    for chrom in chromosomes:
        variants = []
        position = []
        for var in variant_dict.keys():
            if var.startswith(str(str(chrom)+':')):
                variants.append(var)
                position.append(int(var.rsplit(':')[1]))
                
        position, variants = zip(*sorted(zip(position, variants)))
        variants_neworder += variants
        
    return [variant_dict[var].replace('23:','X').replace('24:','Y') for var in variants_neworder]

variants_gwas = {}
all_variants = []
all_cts = []
all_traits = []
for ct in gene_sets:
    gs = gene_sets[ct]
    
    dfxx = df_eqtl_selection[(df_eqtl_selection['gene'].isin(gs))&(df_eqtl_selection['group'] == f'{ct}-cells_healthy')]
    variants = list(set(dfxx['variant']))
    variants = order_variants(variants)
    variants_gwas[ct] = variants
    all_variants += variants
    all_cts = all_cts+[ct]*len(variants)

df_data = pd.DataFrame(index=all_variants)
for variant in all_variants:
    location = variant.rsplit('_')[0]
    
    gwas_dict = get_average_pval_gwas(df_gwas, location)
    
    for trait in gwas_dict:
        if trait not in df_data.columns:
            df_data[trait] = np.zeros(df_data.shape[0])
            
        df_data.loc[variant,trait] = gwas_dict[trait]
        
df_data.index = [i.strip('_') for i in df_data.index]
       
df_data.columns = sorted(df_data.columns)
df_meta = pd.DataFrame(index=all_variants)
df_meta['cell_type'] = all_cts

df_meta.index = [i.strip('_') for i in df_meta.index]

fig = plt.figure(figsize=(6, 4), dpi=300)

colors = ['#ff7f0e','#e377c2','#aa40fc','#1f77b4']

row_ha = pch.HeatmapAnnotation(Celltype=pch.anno_simple(df_meta.cell_type, colors=colors),axis=0)

cm = pch.ClusterMapPlotter(
            data=df_data, left_annotation=row_ha,
            row_cluster=False, col_cluster=False,
            label='Avg -log10(p-val)\nGWAS',row_dendrogram=False,show_rownames=True,show_colnames=True,
            vmin=0,vmax=20,
            rasterized=True,
            tree_kws={'row_cmap': 'Set1'},verbose=0,legend_gap=10,
            #annot=True,linewidths=0.05,linecolor='gold',
            cmap='Purples',
            xlabel="GWAS traits")

fig.savefig(filenameout, bbox_inches='tight')