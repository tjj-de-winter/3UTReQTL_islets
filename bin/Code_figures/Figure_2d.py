#!/bin/python

###########################
##    eQTL Figure 2D     ##
###########################

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

filenameout = 'figure2d.png'
filenameout_anno = 'figure2d_annotations.png'

def jaccard_index(g1_pvals, g2_pvals):
    from sklearn.metrics import jaccard_score
    
    def significant(pval):
        if pval >= -np.log10(0.05):
            return 1
        else:
            return 0
        
    g1 = [significant(pval) for pval in g1_pvals]
    g2 = [significant(pval) for pval in g2_pvals]

    jaccard = jaccard_score(g1, g2)
    return jaccard

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

df_eqtl_beta_T2D = df_eqtl[df_eqtl['type'] == 'beta-cells_T2D']

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/downsampling/CSV_alpha_T2D_n133/alpha-cells_T2D_downsampling_eQTL.adjusted-pvalue.csv'
df_eqtl_alpha_T2D = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df_eqtl_alpha_T2D['type'] = df_eqtl_alpha_T2D['type'].apply(lambda x: str(x+'_downsample'))

fig, ax = shared_unique_eqtls_plot(df_eqtl_beta_T2D, df_eqtl_alpha_T2D, '\nbeta-cell T2D', 'alpha-cell T2D\n(downsampled)', 
                         fontsize=15,xlim=(-1,20), ylim=(-1,18), text=False)

fig.savefig(filenameout, bbox_inches='tight')

fig, ax = shared_unique_eqtls_plot(df_eqtl_beta_T2D, df_eqtl_alpha_T2D, '\nbeta-cell T2D', 'alpha-cell T2D\n(downsampled)', 
                         fontsize=15,xlim=(-1,20), ylim=(-1,18), text=True)

fig.savefig(filenameout_anno, bbox_inches='tight')