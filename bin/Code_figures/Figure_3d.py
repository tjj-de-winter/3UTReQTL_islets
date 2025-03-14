#!/bin/python

###########################
##    eQTL Figure 3D     ##
###########################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as shc
import PyComplexHeatmap as pch

filenameout = 'figure3d.png'

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/deseq2_norm_data_AA.tsv'
df_ct = pd.read_csv(file, index_col=0, sep="\t")
df_ct.columns = [col.rsplit('.')[0] for col in df_ct]

df_ct_norm = np.log(df_ct+1)

log_df = df_ct_norm

nddf = log_df[[c for c in log_df.columns if 'ND' in c]]
t2ddf = log_df[[c for c in log_df.columns if 'T2D' in c]]

def bootstrap_variance(v):
    bootstrap_var = np.array([np.random.choice(v, size = len(v)).var() for i in range(100)])
    return pd.Series({'BSvar': bootstrap_var.mean(), 'BSvar_err': bootstrap_var.std()/np.sqrt(len(v))})

var_nndf = pd.DataFrame(log_df[nddf.columns].apply(lambda x: bootstrap_variance(x), axis = 1))
var_t2ddf = pd.DataFrame(log_df[t2ddf.columns].apply(lambda x: bootstrap_variance(x), axis = 1))
var_bs = var_nndf.merge(var_t2ddf, how = 'inner', left_index = True, right_index = True, suffixes = ['_ND','_T2D'])

fontsize = 8
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

# Hierarchical clustering of miRNAs

# 53 miRNAs as computed in ~/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/miRNA_seq_ND_T2D-AA.ipynb
ND_mirnas = ['hsa-miR-1224-3p', 'hsa-miR-126-5p', 'hsa-miR-127-5p', 'hsa-miR-129-1-3p', 'hsa-miR-1306-5p', 'hsa-miR-132-3p', 'hsa-miR-1343-3p', 'hsa-miR-139-5p', 'hsa-miR-145-3p', 'hsa-miR-145-5p', 'hsa-miR-146a-5p', 'hsa-miR-148a-3p', 'hsa-miR-150-5p', 'hsa-miR-181c-5p', 'hsa-miR-187-3p', 'hsa-miR-18a-3p', 'hsa-miR-194-3p', 'hsa-miR-205-5p', 'hsa-miR-212-5p', 'hsa-miR-26b-5p', 'hsa-miR-27b-5p', 'hsa-miR-28-5p', 'hsa-miR-301a-3p', 'hsa-miR-3064-5p', 'hsa-miR-30a-3p', 'hsa-miR-30a-5p', 'hsa-miR-30c-2-3p', 'hsa-miR-31-5p', 'hsa-miR-369-3p', 'hsa-miR-370-3p', 'hsa-miR-378a-5p', 'hsa-miR-410-3p', 'hsa-miR-424-3p', 'hsa-miR-493-5p', 'hsa-miR-501-5p', 'hsa-miR-505-5p', 'hsa-miR-532-5p', 'hsa-miR-551b-3p', 'hsa-miR-625-5p', 'hsa-miR-6743-3p', 'hsa-miR-6858-3p', 'hsa-miR-769-3p', 'hsa-miR-769-5p', 'hsa-miR-7706', 'hsa-miR-96-5p', 'hsa-miR-98-5p', 'hsa-miR-99a-5p', 'hsa-miR-99b-3p', 'miRNA-2', 'miRNA-31', 'miRNA-40', 'miRNA-47', 'miRNA-51']

df = df_ct_norm.loc[ND_mirnas,:]

# calculate z-scores
df = df.T
df = (df - df.mean()) / df.std()
df = df.T

heatmap_label = 'z-score'
data = df.to_numpy()

dend_miR = shc.dendrogram(shc.linkage(data, method='ward'), get_leaves=True, no_plot=True)
df_heatmap = df.iloc[dend_miR['leaves'],:].copy()

df_heatmap = df_heatmap.T

df_heatmap.index = [i.replace('ND','').replace('T2D','') for i in df_heatmap.index]

# generate additional heatmaps containing group annotation
df_annotation = pd.DataFrame(index=df_heatmap.index)
df_annotation['group'] = [i[4::] for i in df_ct_norm.columns]

# generate additional heatmaps containing logFC and variance in each group
df_heatmaps2 = pd.DataFrame(index=df.index)
df_heatmaps2['ND'] = df.loc[:,['R155ND', 'R167ND', 'R171ND']].mean(axis=1)
df_heatmaps2['T2D'] = df.loc[:,['R185T2D', 'R195T2D', 'R201T2D']].mean(axis=1)
df_heatmaps2['logFC'] = np.log2(df_heatmaps2['ND']/df_heatmaps2['T2D'])
df_heatmaps2['ND_var'] = var_bs['BSvar_ND']
df_heatmaps2['T2D_var'] = var_bs['BSvar_T2D']
df_heatmaps2['label'] = df_heatmaps2.index

# change index names
df_heatmaps2.index = [idx.replace('hsa-','') for idx in df_heatmaps2.index]
df_heatmap.columns = [col.replace('hsa-','') for col in df_heatmap.columns]

row_ha = pch.HeatmapAnnotation(group=pch.anno_simple(df_annotation.group, 
                                                     add_text=False,legend=True, 
                                                     colors=['#F3DD89', '#7D7F99']), 
                               label_kws={'color':'white'},
                               axis=0)

col_ha = pch.HeatmapAnnotation(ND=pch.anno_simple(df_heatmaps2.loc[:,'ND_var'], 
                                                  cmap='Spectral_r',
                                                  vmin=0,
                                                  vmax=1,
                                                  label='var(ND)',
                                                  legend=True),
                               T2D=pch.anno_simple(df_heatmaps2.loc[:,'T2D_var'], 
                                                  cmap='Spectral_r',
                                                  vmin=0,
                                                  vmax=1,
                                                  legend=False),
                               label_kws={'color':'black'},
                               label_side='left',
                               axis=1)

fig = plt.figure(dpi=400, figsize=(6,1.8))
cm1 = pch.ClusterMapPlotter(data=df_heatmap,top_annotation=col_ha,right_annotation=row_ha,
                            col_cluster=False,
                            row_cluster=False,
                            label=heatmap_label,
                            col_dendrogram=False,
                            show_rownames=True,
                            show_colnames=True,
                            row_names_side='left',
                            row_split=3,
                            row_split_gap=1,
                            vmin=-2,
                            vmax=2,
                            tree_kws={'col_cmap': 'Set1','colors':'blue'},
                            verbose=0,
                            xlabel='miRNA',
                            ylabel='Donor',
                            cmap='viridis',xticklabels_kws={'labelrotation':90,'labelcolor':'black'})

fig.savefig(filenameout, bbox_inches='tight')