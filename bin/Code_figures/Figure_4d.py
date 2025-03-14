#!/bin/python

###########################
##    eQTL Figure 4D     ##
###########################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

filenameout = 'figure4d.png'

csv = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-luciferase/endoc_miRNA-inhibitor_qPCR.csv'

fontsize = 18
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

df = pd.read_csv(csv, index_col=0, sep=',')

df = df.dropna()

y = list(df.mean(axis=1))
x = range(len(y))

fig, ax = plt.subplots(dpi=300, figsize=(4,6))

yerr = list(df.std(axis=1)/np.sqrt(df.shape[1]))

ax.bar(x[0],y[0], yerr=yerr[0], color='#D279AF', edgecolor='black', label='control',hatch='///')
ax.bar(x[1],y[1], yerr=yerr[1], color='#D279AF', edgecolor='black', label='NT miR', hatch="\\\\\\")
ax.bar(x[2:],y[2:], yerr=yerr[2:], color='#D279AF', edgecolor='black', label='inhibitor')


ys = df.values

for xi, y in enumerate(ys):
    random = np.random.uniform(-0.25,0.25, 3)
    xscatter = [xi + r for r in random]
    ax.scatter(xscatter,y, s=20,color='black', edgecolor='black')


labels = ['', '', '10', '15', '20', '30']

ax.set_xticks(x,labels, rotation=90)

ax.set_ylabel('PTEN normalized expression')

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel('hsa-miR-127-5p (nM)')

fig.savefig(filenameout, bbox_inches='tight')