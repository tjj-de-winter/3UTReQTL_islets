#!/bin/python

###########################
##    eQTL Figure 2G     ##
###########################

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import random

filenameout = 'figure2g.png'
filenameout_L = 'figure2g_left.png'
filenameout_R = 'figure2g_right.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

def shared_unique_eqtls_plot(df1, df2, group1, group2, xlim=(-1,17), ylim=(-1,17), fontsize = 15, plot=True, text=False):
    variants = set(list(df1['variant']) + list(df2['variant']))

    
    fontsize = fontsize
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['xtick.labelsize'] = fontsize
    plt.rcParams['ytick.labelsize'] = fontsize
    plt.rcParams['figure.titlesize'] = fontsize
    plt.rcParams['font.size'] = fontsize
    
    g1_g2 ={}
    for var in variants:
    #     print(df_eqtl_beta_ND[df_eqtl_beta_ND['variant'] == var]['pvalue'])
        try:
            g1 = float(df1[df1['variant'] == var]['pvalue'])
            g1 = -np.log10(g1)
        except:
            g1 = -np.log10(1)


        try:
            g2 = float(df2[df2['variant'] == var]['pvalue'])
            g2 = -np.log10(g2)
        except:
            g2 = -np.log10(1)


        if g2 >= -np.log10(0.05) or g1 >= -np.log10(0.05):
            g1_g2[var] = [g1, g2]

    g1_pvals = [g1_g2[var][0] for var in g1_g2.keys()]
    g2_pvals = [g1_g2[var][1] for var in g1_g2.keys()]
    variants = [var for var in g1_g2.keys()]
    
    shared = []
    g1_unique = []
    g2_unique = []
    
    if plot:
        fig, ax = plt.subplots(dpi=300, figsize=(4,4))
        # fig, ax = plt.subplots(dpi=300, figsize=(10,10))

        for g1_, g2_, var in zip(g1_pvals, g2_pvals, variants):
            if g1_ >= -np.log10(0.05) and g2_ >= -np.log10(0.05):  
                ax.scatter(g1_, g2_, edgecolors='black', color='#e0a604')
                if text:
                    ax.text(g1_,g2_,var.rsplit('_')[2], size=3)
                shared.append(var)
            if g1_ >= -np.log10(0.05) and g2_ <= -np.log10(0.05):
                ax.scatter(g1_, g2_, edgecolors='black', color='#13b807')
                g1_unique.append(var)
            if g1_ <= -np.log10(0.05) and g2_ >= -np.log10(0.05):
                ax.scatter(g1_, g2_, edgecolors='black', color='#6e50f2')
                g2_unique.append(var)
    else:
        for g1_, g2_, var in zip(g1_pvals, g2_pvals, variants):
            if g1_ >= -np.log10(0.05) and g2_ >= -np.log10(0.05):  
                shared.append(var)
            if g1_ >= -np.log10(0.05) and g2_ <= -np.log10(0.05):
                g1_unique.append(var)
            if g1_ <= -np.log10(0.05) and g2_ >= -np.log10(0.05):
                g2_unique.append(var)

    shared = list(set(shared))
    g2_unique = list(set(g2_unique))
    g1_unique = list(set(g1_unique))
        
    if not plot:
        return g1_unique, g2_unique, shared
    
    jaccard = jaccard_index(g1_pvals, g2_pvals)
    

    alleqtls = len(shared) + len(g2_unique) + len(g1_unique)
    print('shared', len(shared)/alleqtls * 100, f'jaccard index: {jaccard}')
    print(f'{group1}_unique', len(g1_unique)/alleqtls * 100)
    print(f'{group2}_unique', len(g2_unique)/alleqtls * 100)


    ax.set_xlabel(f'-log10(p-value) {group1} eQTL')
    ax.set_ylabel(f'-log10(p-value) {group2} eQTL')
    ax.axvline(-np.log10(0.05), color='red')
    ax.axhline(-np.log10(0.05), color='red')
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])

    ax.grid(False)
    
    return fig, ax

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv'
df_eqtl = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df_eqtl['gene'] = df_eqtl['gene'].apply(lambda x: x.rsplit('>')[0])

df_eqtl = df_eqtl.drop_duplicates()

df_eqtl_beta_ND = df_eqtl[df_eqtl['type'] == 'beta-cells_healthy']
df_eqtl_beta_T2D = df_eqtl[df_eqtl['type'] == 'beta-cells_T2D']

# Extract eGenes
ND_unique, T2D_unique, shared = shared_unique_eqtls_plot(df_eqtl_beta_ND, df_eqtl_beta_T2D, 'beta-cell ND', 'beta-cell T2D', plot=False)


#volcano plot with eGenes significant in beta cells

adata_sub = adata[adata.obs['Cell_type'] == 'beta-cells']

# genes = [var.rsplit('_')[2] for var in variants]

sc.tl.rank_genes_groups(adata_sub, 'Group', method='wilcoxon')
result = adata_sub.uns['rank_genes_groups']
group='T2D'
df = pd.DataFrame({str('logfoldchanges'): result['logfoldchanges'][group],
                   str('pvals_adj'): result['pvals_adj'][group]}, index=result['names'][group])

df.columns = ['log2(foldchanges)','pvals_adj_wilcoxon']

fig, ax = plt.subplots(dpi=300, figsize=(7,5))

total_genes = []

#ND unique
genes = [var.rsplit('_')[2] for var in ND_unique]
df_sub = df.loc[set(genes),:]

total_genes += genes

x = list(df_sub['log2(foldchanges)'])
y = [-np.log10(p) for p in df_sub['pvals_adj_wilcoxon']]

ax.scatter(x,y, edgecolors='black', color='#13b807')

#T2D unique
genes = [var.rsplit('_')[2] for var in T2D_unique]
df_sub = df.loc[set(genes),:]

total_genes += genes

x = list(df_sub['log2(foldchanges)'])
y = [-np.log10(p) for p in df_sub['pvals_adj_wilcoxon']]

ax.scatter(x,y, edgecolors='black', color='#6e50f2')

#shared
genes = [var.rsplit('_')[2] for var in shared]
df_sub = df.loc[set(genes),:]

total_genes += genes

x = list(df_sub['log2(foldchanges)'])
y = [-np.log10(p) for p in df_sub['pvals_adj_wilcoxon']]

ax.scatter(x,y, edgecolors='black', color='#e0a604')

ax.set_title('beta-cells eGenes')
ax.set_xlabel('log2(foldchanges)')
ax.set_ylabel('-log10(p-value)')

ax.axhline(-np.log10(0.05), color='black', linestyle='--')
ax.axvline(0, color='black', linestyle='--')

ax.set_xlim(-2.5,2.5)
ax.set_ylim(0,12)
ax.grid(False)

# annotation

def annotate(x,y,s, xlim=(-25,25), ymax=3, distance_ratio=1, hight_ratio=1):
    if x < 0:
        if x < xlim[0]:
            xi = x + random.uniform(xlim[1]/8.33*distance_ratio, xlim[1]/6.25*distance_ratio)
        else:
            xi = x - random.uniform(xlim[1]/3.125*distance_ratio, xlim[1]/2.78*distance_ratio)
    else:
        if x > xlim[1]:
            xi = x - random.uniform(xlim[1]/3.125*distance_ratio, xlim[1]/2.78*distance_ratio)
        else:
            xi = x + random.uniform(xlim[1]/8.33*distance_ratio, xlim[1]/6.25*distance_ratio)
            
    yi = ymax+1
    
    def ybetween(y):
        if y > ymax:
            return True
        elif y < 0:
            return True
    
    while ybetween(yi):
        yi = y + random.uniform(-ymax/12*hight_ratio, ymax/12*hight_ratio)
        
    ax.annotate(s, xy=(x,y), xytext=(xi,yi), arrowprops=dict(facecolor='black', arrowstyle='->'), size=10)
    
    
dfa = df.loc[set(total_genes)]
up = dfa[(dfa['pvals_adj_wilcoxon'] < 0.05)&(dfa['log2(foldchanges)'] > 0)].sort_values('log2(foldchanges)', ascending=False).head(5)
down = dfa[(dfa['pvals_adj_wilcoxon'] < 0.05)&(dfa['log2(foldchanges)'] < 0)].sort_values('log2(foldchanges)', ascending=True).head(5)
# interesting = dfa.loc[['PTEN']]

dftext = pd.concat([up, down])#, interesting])

for s in dftext.index:
    x = dftext.loc[s, 'log2(foldchanges)']
    y = -np.log10(dftext.loc[s, 'pvals_adj_wilcoxon'])

    annotate(x,y,s, xlim=(-2,2), ymax=100, distance_ratio=1, hight_ratio=0.1)
    
interesting = dfa.loc[['THAP5', 'PABPC1', 'ATP6V1C1', 'PTEN']]
dftext = interesting 
for s in dftext.index:
    x = dftext.loc[s, 'log2(foldchanges)']
    y = -np.log10(dftext.loc[s, 'pvals_adj_wilcoxon'])

    annotate(x,y,s, xlim=(-2,2), ymax=12, distance_ratio=2.5, hight_ratio=0.1)

fig.savefig(filenameout, bbox_inches='tight')

# pie chart left

# Extract eGenes
ND_unique, T2D_unique, shared = shared_unique_eqtls_plot(df_eqtl_beta_ND, df_eqtl_beta_T2D, 'beta-cell ND', 'beta-cell T2D', plot=False)


def diff_expr(variants, up_down='up', df=df):
    genes = [var.rsplit('_')[2] for var in variants]
    df_sub = df.loc[set(genes),:]
    
    if up_down == 'up':
        df_sub = df_sub[(df_sub['log2(foldchanges)'] > 0) & (df_sub['pvals_adj_wilcoxon'] < 0.05)]
        return df_sub.shape[0]
    
    else:
        df_sub = df_sub[(df_sub['log2(foldchanges)'] < 0) & (df_sub['pvals_adj_wilcoxon'] < 0.05)]
        return df_sub.shape[0]

    
# determine differential expressed cells

adata_sub = adata[adata.obs['Cell_type'] == 'beta-cells']


sc.tl.rank_genes_groups(adata_sub, 'Group', method='wilcoxon')
result = adata_sub.uns['rank_genes_groups']
group='T2D'
df = pd.DataFrame({str('logfoldchanges'): result['logfoldchanges'][group],
                   str('pvals_adj'): result['pvals_adj'][group]}, index=result['names'][group])

df.columns = ['log2(foldchanges)','pvals_adj_wilcoxon']

# count down regulated genes per eQTL type
up_down = 'down'

variant_lists = ND_unique, T2D_unique, shared

data = []
for variants in variant_lists:
    data.append(diff_expr(variants, up_down=up_down))

fig, ax = plt.subplots(dpi=300)

colors = ['#13b807', '#6e50f2', '#e0a604']

ax.pie(data, colors=colors, wedgeprops={'edgecolor': 'black', 'linewidth': 2});

# count up regulated genes per eQTL type
up_down = 'up'

variant_lists = ND_unique, T2D_unique, shared

data = []
for variants in variant_lists:
    data.append(diff_expr(variants, up_down=up_down))

fig.savefig(filenameout_L, bbox_inches='tight')

fig, ax = plt.subplots(dpi=300)

colors = ['#13b807', '#6e50f2', '#e0a604']

ax.pie(data, colors=colors, wedgeprops={'edgecolor': 'black', 'linewidth': 2});

fig.savefig(filenameout_R, bbox_inches='tight')

