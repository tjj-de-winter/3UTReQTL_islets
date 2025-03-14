#!/bin/python

###########################
##    eQTL Figure 4B     ##
###########################

import matplotlib.pyplot as plt
import random
import pandas as pd
import numpy as np

filenameout = 'figure4b.png'

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/miRNA qPCR/miRNAqPCR_normalized.csv'

df = pd.read_csv(file, index_col=0)

fig, ax = plt.subplots(dpi=300, figsize=(2,3))

y = list(df.mean(axis=1))
x = range(len(y))
error_bar = df.std(axis=1)/np.sqrt(df.shape[1])


colors = ['#5A7EB7', '#D279AF']
for i, idx in enumerate(df.index):
    ax.bar([i],y[i], yerr=error_bar[i], color=colors[i], edgecolor='black')
#     ax.bar([i],y[i], yerr=error_bar[i], color='#D279AF', edgecolor='black')
    
    ys = df.loc[idx,:]
    xs = [i+random.uniform(-0.2, 0.2) for j in range(3)]
    ax.scatter(xs, ys, s=10,edgecolor='black', color='black')

ax.set_ylabel(f'miR-127-5p\nnormalized expression')
ax.set_xticks([0,1])
ax.set_xticklabels(['HEK293T', 'EndoC-βH1'])

fig.savefig(filenameout, bbox_inches='tight')