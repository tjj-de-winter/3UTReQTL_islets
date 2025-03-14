#!/bin/python

###########################
##    eQTL Figure 3G     ##
###########################

import matplotlib.pyplot as plt
import pandas as pd

filenameout = 'figure3g.png'

# from '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/mirbase/miranda.out.filtered.beta-cells.healthy.csv'
allele1 = [140]
allele2 = [108]

fontsize = 15
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

fig, ax = plt.subplots(dpi=300, figsize=(3,5))

ax.bar([0],allele1, color='#a8329d', edgecolor='black')
ax.bar([1],allele2, color='#b5c229', edgecolor='black')
ax.axhline(140, linestyle='--', color='red')
ax.set_ylim(20,180)
ax.set_xticks([0,1], ['T-allele','C-allele'], rotation=90)
ax.set_ylabel('miRanda binding score')
ax.set_title('hsa-miR-127-5p\n rs701848')

fig.savefig(filenameout, bbox_inches='tight')