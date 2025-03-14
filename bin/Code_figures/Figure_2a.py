#!/bin/python

###########################
##    eQTL Figure 2A     ##
###########################

import matplotlib.pyplot as plt
import pandas as pd
import glob
from scipy.stats import pearsonr

filenameout = 'figure2a.png'

# ND variants

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv'
df = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df = df.drop_duplicates()

celltypes = [ct.rsplit('_')[0] for ct in set(df.type)]

# T2D variants
file_path = 'Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/downsampling/T2D_vs_ND/all_celltypes/*_downsampling_eQTL.adjusted-pvalue.csv'
files = glob.glob('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/downsampling/T2D_vs_ND/all_celltypes/*_downsampling_eQTL.adjusted-pvalue.csv'
)

groups_downsampled = {'acinar-cells': 'healthy', 'ductal-cells':'healthy', 'endothelial-cells':'healthy', 'PP-cells':'healthy','stellate-cells': 'healthy', 'alpha-cells':'healthy', 'beta-cells':'T2D', 'delta-cells':'healthy'}

def convert(group):
    if group == 'healthy':
        return 'T2D'
    else:
        return 'healthy'

def df_ct_specific(df, celltype):
    df_sub = df[df['type'] == celltype]
    return df_sub

def ND_vs_T2D(df1, df2, corr='jaccard'):
    variants = list(set(df1.variant).intersection(set(df2.variant)))
    
    df1.index = df1.variant
    df2.index = df2.variant
                
    g1_pvals = list(df1.loc[variants,'pvalue'])
    g2_pvals = list(df2.loc[variants,'pvalue'])
     
    if corr == 'jaccard':
        jaccard = jaccard_index(g1_pvals, g2_pvals)
        return jaccard
    elif corr == 'pearson':
        pearson = pearson_coeff(g1_pvals, g2_pvals)
        return pearson

def pearson_coeff(g1_pvals, g2_pvals):
    correlation, p_value = pearsonr(g1_pvals, g2_pvals)
    return correlation

eqtl_dict_dfs = {}
for ct in celltypes:
    group_down = groups_downsampled[ct]
    group = convert(group_down)
    
    nd_t2d = {}
    
    df1 = df_ct_specific(df, f'{ct}_{group}')
    nd_t2d[group] = df1
    for file in files:
        if ct in file:
            df2 = pd.read_csv(file, sep=',', names=['variant', 'gene','pvalue', 'type'])
            df2 = df2.drop_duplicates()
            nd_t2d[group_down] = df2
            break
    eqtl_dict_dfs[ct] = [nd_t2d['healthy'], nd_t2d['T2D']]

fontsize = 18
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize
    

data = []
labels = ['beta-cells', 'alpha-cells', 'delta-cells', 'PP-cells']
for ct in labels:
    df1 = eqtl_dict_dfs[ct][0]
    df2 = eqtl_dict_dfs[ct][1]
    pearson = ND_vs_T2D(df1, df2, corr='pearson')
    
    data.append(pearson)
    
fig, ax = plt.subplots(dpi = 300, figsize=(5,5))

colors_cells = {'acinar-cells': '#1f77b4', 'ductal-cells':'#ff7f0e', 'endothelial-cells':'#279e68', 'PP-cells':'#d62728','stellate-cells': '#aa40fc', 'alpha-cells':'#8c564b', 'beta-cells':'#e377c2', 'delta-cells':'#b5bd61'}
colors = [colors_cells[ct] for ct in labels]

x = range(len(labels))
ax.bar(x, data, color=colors, edgecolor='black')

ax.set_xticks(x,labels, rotation=90)
ax.set_ylabel('Pearson coefficient')
ax.set_title('ND vs T2D')
ax.set_ylim(0,1.1)

fig.savefig(filenameout, bbox_inches='tight')