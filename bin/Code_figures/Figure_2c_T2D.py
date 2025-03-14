#!/bin/python

###########################
##    eQTL Figure 2C     ##
###########################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

filenameout = 'figure2c_T2D.png'

print(filenameout)

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv'
df = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df = df.drop_duplicates()

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/downsampling/CSV_alpha_T2D_n133/alpha-cells_T2D_downsampling_eQTL.adjusted-pvalue.csv'
df_eqtl_alpha_T2D = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df_eqtl_alpha_T2D['type'] = df_eqtl_alpha_T2D['type'].apply(lambda x: str(x+'_downsample'))

df = pd.concat([df, df_eqtl_alpha_T2D])
df.index = range(df.shape[0])

idxs = list(set(df['variant']))
cols = list(set(df['type']))

matrix = pd.DataFrame(np.zeros((len(idxs), len(cols))), index=idxs, columns=cols)

for idx in df.index:
    ct = df.iloc[idx, df.columns.to_list().index('type')]
    var = df.iloc[idx, df.columns.to_list().index('variant')]
    pval = df.iloc[idx, df.columns.to_list().index('pvalue')]
    pval = -np.log10(pval)
    matrix.loc[var,ct] = pval

correlation_matrix = matrix.corr(method='pearson')

cells = ['beta-cells_T2D', 'alpha-cells_T2D_downsample','delta-cells_T2D','PP-cells_T2D']

correlation_matrix = correlation_matrix.loc[cells, cells]

fontsize = 20
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

plt.rcParams['font.family'] = 'Arial'  # Set Arial globally for this session

fig, ax = plt.subplots(figsize=(5, 4),dpi=300)  # Adjust size as needed

# Create the heatmap
im = ax.imshow(correlation_matrix, interpolation='nearest', cmap='coolwarm', vmax=0.49)

# Add a colorbar
colorbar = plt.colorbar(im, ax=ax, label='Pearson correlation')
colorbar.set_label('Pearson correlation', fontsize=25)
colorbar.ax.tick_params(labelsize=25)


# Set ticks and labels
labels = [ct.replace('healthy', 'ND') for ct in cells]
labels = [ct.replace('_downsample', '\n(downsampled)') for ct in labels]
labels = [ct.replace('_', ' ') for ct in labels]

ax.set_xticks(np.arange(len(labels)))
ax.set_yticks(np.arange(len(labels)))
ax.set_xticklabels(labels, rotation=90, fontsize=25)
ax.set_yticklabels(labels,fontsize=25)

# for xi, idx in enumerate(correlation_matrix.index):
#     for yi, col in enumerate(correlation_matrix.columns):
#         r2 = correlation_matrix.loc[idx,col]
#         r2 = round(r2,1)
#         ax.text(xi, yi, r2, color='white', 
#                 horizontalalignment='center',
#                 verticalalignment='center',)

# Add a title
ax.set_title('eQTL correlation', fontsize=25)

# Show the plot
plt.grid(False)

fig.savefig(filenameout, bbox_inches='tight')
