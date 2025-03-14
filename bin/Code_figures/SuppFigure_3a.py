#!/bin/python

###########################
##  eQTL SUPP Figure 3A  ##
###########################

import scipy.cluster.hierarchy as shc
import PyComplexHeatmap as pch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

filenameout = 'supp_figure3a.png'

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/deseq2_norm_data_AA.tsv'
df_ct = pd.read_csv(file, index_col=0, sep="\t")
df_ct.columns = [col.rsplit('.')[0] for col in df_ct]

df_ct_norm = np.log(df_ct+1)

data = df_ct_norm.to_numpy()

dend_miR = shc.dendrogram(shc.linkage(data, method='ward'), get_leaves=True, no_plot=True)
df_heatmap = df_ct_norm.iloc[dend_miR['leaves']].copy()

df_heatmap.columns = [i.replace('ND','').replace('T2D','') for i in df_heatmap.columns]

df_annotation = pd.DataFrame(index=df_heatmap.columns)
df_annotation['group'] = [i[4::] for i in df_ct_norm.columns]

col_ha = pch.HeatmapAnnotation(group=pch.anno_simple(df_annotation.group, add_text=False,legend=True, colors=['#F3DD89', '#7D7F99']),
                               label_kws={'color':'white'},
                               axis=1)


fig = plt.figure(dpi=300, figsize=(3,10))
cm1 = pch.ClusterMapPlotter(data=df_heatmap,top_annotation=col_ha,
                            col_cluster=False,row_cluster=False,
                            label='Expression',col_dendrogram=False,
                            show_rownames=True,show_colnames=True,
                            vmin=0,
                            vmax=9,
                            row_names_side='right', 
                            tree_kws={'col_cmap': 'Set1','colors':'blue'},verbose=0,legend_gap=20,
                            cmap='viridis',xticklabels_kws={'labelrotation':-90,'labelcolor':'black'})

fig.savefig(filenameout, bbox_inches='tight')