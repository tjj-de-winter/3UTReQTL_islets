#!/bin/python

###########################
##    eQTL Figure 3E     ##
###########################

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import PyComplexHeatmap as pch

filenameout = 'figure3e.png'

files = glob.glob('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/mirbase/miranda.out.filtered.*cells*.csv')
df_dict = {}
for file in files:
    split = file.rsplit('.')
    name = '_'.join([split[3],split[4]])
    df = pd.read_csv(file, index_col=0)
    df.index = range(df.shape[0])
    df['mirna'] = [mirna.strip('>') for mirna in df['mirna']]
    df = df[(df['slope_binding'] > 0)&(df['slope_binding_pval'] < 0.05)]
    df['variant-miRNA'] = ['-'.join([v,m]) for v, m in zip(df['variant'], df['mirna'])]
    df['cell_group'] = df['Target'].apply(lambda x: name)
    df['group'] = df['cell_group'].apply(lambda x: x.rsplit('_')[1])
    df_dict[name] = df

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv'
df_eqtl = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df_eqtl['gene'] = df_eqtl['gene'].apply(lambda x: x.rsplit('>')[0])

df_eqtl = df_eqtl.drop_duplicates()

df_eqtl_beta_ND = df_eqtl[df_eqtl['type'] == 'beta-cells_healthy']
df_eqtl_beta_T2D = df_eqtl[df_eqtl['type'] == 'beta-cells_T2D']

variants = set(list(df_eqtl_beta_ND['variant']) + list(df_eqtl_beta_T2D['variant']))

ND_T2D ={}
for var in variants:
#     print(df_eqtl_beta_ND[df_eqtl_beta_ND['variant'] == var]['pvalue'])
    try:
        ND = float(df_eqtl_beta_ND[df_eqtl_beta_ND['variant'] == var]['pvalue'])
        ND = -np.log10(ND)
    except:
        ND = -np.log10(1)

    
    try:
        T2D = float(df_eqtl_beta_T2D[df_eqtl_beta_T2D['variant'] == var]['pvalue'])
        T2D = -np.log10(T2D)
    except:
        T2D = -np.log10(1)
        
    
    if T2D >= -np.log10(0.05) or ND >= -np.log10(0.05):
        ND_T2D[var] = [ND, T2D]
        
ND_pvals = [ND_T2D[var][0] for var in ND_T2D.keys()]
T2D_pvals = [ND_T2D[var][1] for var in ND_T2D.keys()]
variants = [var for var in ND_T2D.keys()]

shared = []
T2D_unique = []
T2D_unique_p_T2D = []
T2D_unique_p_ND = []
ND_unique = []
ND_unique_p_T2D = []
ND_unique_p_ND = []

for nd, t2d, var in zip(ND_pvals, T2D_pvals, variants):
    if nd >= -np.log10(0.05) and t2d >= -np.log10(0.05):  
#         ax.scatter(nd, t2d, edgecolors='black', color='#e0a604')
        shared.append(var)
    if nd >= -np.log10(0.05) and t2d <= -np.log10(0.05):
#         ax.scatter(nd, t2d, edgecolors='black', color='#13b807')
        ND_unique.append(var)
        ND_unique_p_T2D.append(t2d)
        ND_unique_p_ND.append(nd)
    if nd <= -np.log10(0.05) and t2d >= -np.log10(0.05):
#         ax.scatter(nd, t2d, edgecolors='black', color='#6e50f2')
        T2D_unique.append(var)
        T2D_unique_p_T2D.append(t2d)
        T2D_unique_p_ND.append(nd)
        
shared = list(set(shared))
T2D_unique = list(set(T2D_unique))
ND_unique = list(set(ND_unique))

fontsize = 8
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize
plt.rcParams['figure.dpi'] = 300  # Set global DPI for the figure


def variant_conversion(var):
    split = var.rsplit('_')
    return '_'.join([split[2],split[0]])

def variant_label(var):
    split = var.rsplit('_')
    return ' '.join([split[2],split[0]])

def eQTL_miRNA(df_mirna, variants, mirnas):
    
    variant_intersection = set(variants).intersection(set(df_mirna['variant']))
    
    selection = df_mirna[(df_mirna['variant'].isin(variant_intersection)) & (df_mirna['mirna'].isin(mirnas))]
    
    print(selection.shape[0])
    
    return selection

def order_variants(variant_list):
    '''Order variants based on chromosome number and nucleotide position'''
    variant_dict = {v.rsplit('_')[1].replace('X:','23:'): v for v in variant_list}
    
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
        
    return [variant_dict[var] for var in variants_neworder]

# miRNAs from miRNA_seq_ND_T2D_AA
ND_mirnas = ['hsa-miR-1224-3p', 'hsa-miR-126-5p', 'hsa-miR-127-5p', 'hsa-miR-129-1-3p', 'hsa-miR-1306-5p', 'hsa-miR-132-3p', 'hsa-miR-1343-3p', 'hsa-miR-139-5p', 'hsa-miR-145-3p', 'hsa-miR-145-5p', 'hsa-miR-146a-5p', 'hsa-miR-148a-3p', 'hsa-miR-150-5p', 'hsa-miR-181c-5p', 'hsa-miR-187-3p', 'hsa-miR-18a-3p', 'hsa-miR-194-3p', 'hsa-miR-205-5p', 'hsa-miR-212-5p', 'hsa-miR-26b-5p', 'hsa-miR-27b-5p', 'hsa-miR-28-5p', 'hsa-miR-301a-3p', 'hsa-miR-3064-5p', 'hsa-miR-30a-3p', 'hsa-miR-30a-5p', 'hsa-miR-30c-2-3p', 'hsa-miR-31-5p', 'hsa-miR-369-3p', 'hsa-miR-370-3p', 'hsa-miR-378a-5p', 'hsa-miR-410-3p', 'hsa-miR-424-3p', 'hsa-miR-493-5p', 'hsa-miR-501-5p', 'hsa-miR-505-5p', 'hsa-miR-532-5p', 'hsa-miR-551b-3p', 'hsa-miR-625-5p', 'hsa-miR-6743-3p', 'hsa-miR-6858-3p', 'hsa-miR-769-3p', 'hsa-miR-769-5p', 'hsa-miR-7706', 'hsa-miR-96-5p', 'hsa-miR-98-5p', 'hsa-miR-99a-5p', 'hsa-miR-99b-3p', 'miRNA-2', 'miRNA-31', 'miRNA-40', 'miRNA-47', 'miRNA-51']


# select miRNA-eQTL pairs that are specific for ND and not for T2D

variants = [variant_conversion(var) for var in ND_unique]
variant_convert = {variant_conversion(var): var for var in ND_unique}

df_ND = eQTL_miRNA(df_dict['beta-cells_healthy'], variants, ND_mirnas)
df_ND['location'] = df_ND['variant'].apply(lambda var: var.rsplit('_')[1])
df_ND = df_ND.sort_values(by=['mirna'])

# order variants by miRNA and for each miRNA by variant position
dfs = []
for mir, dfi in df_ND.groupby('mirna'):
    dfi.index = dfi.variant
    
    index = order_variants(list(dfi.variant))
    dfi = dfi.loc[index,:]
    
    dfs.append(dfi)
    
df_ND_sorted = pd.concat(dfs)
df_ND_sorted.index = range(df_ND_sorted.shape[0])

# dataframes and matrix for plotting

index = []
for idx in df_ND_sorted['mirna']:
    if idx not in index:
        index.append(idx)
columns = []
for col in df_ND_sorted['variant']:
    if col not in columns:
        columns.append(col)
        
zeros = np.zeros((len(index), len(columns)))
matrix = pd.DataFrame(zeros, index=index, columns=columns)

for mir in matrix.index:
    for var in matrix.columns:
        try:
            var_mirna = '-'.join([var,mir])
            val = float(df_ND[df_ND['variant-miRNA'] == var_mirna]['Score'])
            matrix.loc[mir,var] = val
        except:
            continue
            
# eQTL pvalue annotation for ND and T2D beta cells
df_anno = pd.DataFrame(index=matrix.columns, columns=['ND_p', 'T2D_p'],dtype=float)
df_eqtl['location'] = df_eqtl['variant'].apply(lambda var: var.rsplit('_')[0])
for var in df_anno.index:
    loc = var.rsplit('_')[1]
    groups = ['beta-cells_healthy', 'beta-cells_T2D']
    group_col = ['ND_p', 'T2D_p']
    for group, col in zip(groups, group_col): 
        try:
            pval = float(df_eqtl[(df_eqtl['type']==group)&(df_eqtl['location']==loc)]['pvalue'])
        except:
            pval = 1.0
        pval = -np.log10(pval)
        if pval == -0.0:
            pval = 0
        df_anno.loc[var,col] = pval

# plot the data
fig = plt.figure(dpi=300, figsize=(6,2.8))

# change label notation
matrix.index = [idx.replace('hsa-','') for idx in matrix.index]
matrix.columns = [col.replace('_',' ') for col in matrix.columns]
df_anno.index = [idx.replace('_',' ') for idx in df_anno.index]

col_ha = pch.HeatmapAnnotation(ND=pch.anno_simple(df_anno.ND_p, 
                                                  add_text=False,
                                                  legend=True, 
                                                  cmap='BrBG',
                                                  vmin=0,
                                                  vmax=6),
                               T2D=pch.anno_simple(df_anno.T2D_p, 
                                                   add_text=False,
                                                   legend=False, 
                                                   cmap='BrBG',
                                                   vmin=0,
                                                   vmax=6),
                               label_kws={'color':'black'},
                               axis=1)

cm1 = pch.ClusterMapPlotter(data=matrix,top_annotation=col_ha,
                            col_cluster=False,row_cluster=False,
                            label='miRNA binding score',col_dendrogram=False,
                            show_rownames=True,show_colnames=True,
                            row_names_side='left', 
                            cmap='Greens',xticklabels_kws={'labelrotation':45,'labelcolor':'black'},
                            linewidths=0.5,           
                            linecolor='gray',
                            legend_gap=15)

fig.savefig(filenameout, bbox_inches='tight')