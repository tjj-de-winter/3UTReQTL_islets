#!/bin/python

###########################
##  eQTL SUPP Figure 3D  ##
###########################

import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


filenameout = 'supp_figure3d.png'

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/deseq2_norm_data_AA.tsv'
df_ct = pd.read_csv(file, index_col=0, sep="\t")
df_ct.columns = [col.rsplit('.')[0] for col in df_ct]

df_ct_norm = np.log(df_ct+1)

ND_cols = [col for col in df_ct_norm.columns if "ND" in col]
ND = df_ct_norm.loc[:,ND_cols]

T2D_cols = [col for col in df_ct_norm.columns if "T2D" in col]
T2D = df_ct_norm.loc[:,T2D_cols]

ND.to_csv('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/group1.csv')
T2D.to_csv('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/group2.csv')

os.system('Rscript /Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/DESeQ2_DGE.R')

df_DESeq2 = pd.read_csv('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/DESeq2_results.csv', index_col=0)

fig, ax = plt.subplots(dpi=300, figsize=(4,3))

x = df_DESeq2.log2FoldChange
y = -np.log10(df_DESeq2.padj)
ax.scatter(x,y, edgecolor='black', color='#3d3dcc')
ax.axhline(-np.log10(0.05),linestyle='--', color='black')
ax.axvline(0, linestyle='--', color='black')

ax.set_ylim(-.2,4.2)
ax.set_xlabel('log2(foldchanges)')
ax.set_ylabel('-log10(p-value)')
ax.set_title('miRNAs')

fig.savefig(filenameout, bbox_inches='tight')