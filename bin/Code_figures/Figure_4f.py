#!/bin/python

###########################
##    eQTL Figure 4F     ##
###########################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

filenameout1 = 'figure4f_1.png'
filenameout2 = 'figure4f_2.png'
filenameoutlegend = 'figure4f_legend.png'

csv_HEK = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-luciferase/luciferase_HEK.csv'

df = pd.read_csv(csv_HEK, index_col=0, names=['1','2','3','4','5','6','allele'])

fig, ax = plt.subplots(dpi=300, figsize=(4,6))
color_dict = {'C':'#b5c229', 'T':'#a8329d'}

for i, dfi in df.groupby('allele'):
    
    if i == 'C':
        continue
    
    dfi = dfi.iloc[:,:6]
    
    index = [f'Vector {i}', f'vector {i} + NT miRNA', f'vector {i} + miRNA']
    
    dfi = dfi.loc[index,:]
    
    yerr = list(dfi.std(axis=1)/np.sqrt(dfi.shape[1]))
    y = list(dfi.mean(axis=1))
    x = range(len(y))

    ax.bar(x[0],y[0],width=0.4, yerr=yerr[0], color=color_dict[i], edgecolor='black', label=f'control',hatch='///')
    ax.bar(x[1],y[1],width=0.4, yerr=yerr[1], color=color_dict[i], edgecolor='black', label=f'NT miR', hatch="\\\\\\")
    ax.bar(x[2:],y[2:],width=0.4, yerr=yerr[2:], color=color_dict[i], edgecolor='black', label=f'mimic')
    
    ys = dfi.values

    for xi, y in zip(x, ys):
        random = np.random.uniform(-0.1,0.1, len(y))
        xscatter = [xi + r for r in random]
        ax.scatter(xscatter,y, s=20,color='black', edgecolor='black')

# Remove x-axis labels
ax.set_xlabel("")  # Remove x-axis label
ax.set_xticklabels([])  # Remove tick labels

# Remove x-axis tick marks
ax.set_xticks([])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylabel('Normalized luciferase activity')
ax.set_title('rs701848-T')
ax.set_xlabel('hsa-miR-127-5p')

fig.savefig(filenameout1, bbox_inches='tight')

fig, ax = plt.subplots(dpi=300, figsize=(4,6))
color_dict = {'C':'#b5c229', 'T':'#a8329d'}

for i, dfi in df.groupby('allele'):
    
    if i == 'T':
        continue
    
    dfi = dfi.iloc[:,:6]
        
    index = [f'Vector {i}', f'vector {i} + NT miRNA', f'vector {i} + miRNA']
    
    dfi = dfi.loc[index,:]
    
    yerr = list(dfi.std(axis=1)/np.sqrt(dfi.shape[1]))
    y = list(dfi.mean(axis=1))
    x = range(len(y))

    ax.bar(x[0],y[0],width=0.4, yerr=yerr[0], color=color_dict[i], edgecolor='black', label=f'control',hatch='///')
    ax.bar(x[1],y[1],width=0.4, yerr=yerr[1], color=color_dict[i], edgecolor='black', label=f'NT miR', hatch="\\\\\\")
    ax.bar(x[2:],y[2:],width=0.4, yerr=yerr[2:], color=color_dict[i], edgecolor='black', label=f'mimic')
    
    ys = dfi.values

    for xi, y in zip(x, ys):
        random = np.random.uniform(-0.1,0.1, len(y))
        xscatter = [xi + r for r in random]
        ax.scatter(xscatter,y, s=20,color='black', edgecolor='black')

# Remove x-axis labels
ax.set_xlabel("")  # Remove x-axis label
ax.set_xticklabels([])  # Remove tick labels

# Remove x-axis tick marks
ax.set_xticks([])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylabel('Normalized luciferase activity')
ax.set_title('rs701848-C')
ax.set_xlabel('hsa-miR-127-5p')
ax.set_ylim(0,1.42)

fig.savefig(filenameout2, bbox_inches='tight')

fig, ax = plt.subplots(dpi=300, figsize=(4,6))
color_dict = {'C':'#b5c229', 'T':'white'}

for i, dfi in df.groupby('allele'):
    
    if i == 'C':
        continue
    
    dfi = dfi.iloc[:,:6]
    
    index = [f'Vector {i}', f'vector {i} + NT miRNA', f'vector {i} + miRNA']
    
    dfi = dfi.loc[index,:]
    
    yerr = list(dfi.std(axis=1)/np.sqrt(dfi.shape[1]))
    y = list(dfi.mean(axis=1))
    x = range(len(y))

    ax.bar(x[0],y[0],width=0.4, yerr=yerr[0], color=color_dict[i], edgecolor='black', label=f'control',hatch='///')
    ax.bar(x[1],y[1],width=0.4, yerr=yerr[1], color=color_dict[i], edgecolor='black', label=f'NT miR', hatch="\\\\\\")
    ax.bar(x[2:],y[2:],width=0.4, yerr=yerr[2:], color=color_dict[i], edgecolor='black', label=f'mimic')
    
    ys = dfi.values

    for xi, y in zip(x, ys):
        random = np.random.uniform(-0.1,0.1, len(y))
        xscatter = [xi + r for r in random]
        ax.scatter(xscatter,y, s=20,color='black', edgecolor='black')

# Remove x-axis labels
ax.set_xlabel("")  # Remove x-axis label
ax.set_xticklabels([])  # Remove tick labels

# Remove x-axis tick marks
ax.set_xticks([])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_ylabel('Normalized luciferase activity')
ax.set_title('rs701848-C')
ax.set_xlabel('hsa-miR-127-5p')

fig.savefig(filenameoutlegend, bbox_inches='tight')
