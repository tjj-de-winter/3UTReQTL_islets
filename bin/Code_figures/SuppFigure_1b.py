#!/bin/python

###########################
##  eQTL SUPP Figure 1B  ##
###########################

import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np

filenameout = 'supp_figure1b.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

# cell type distribution per donor

donors = sorted(list(set([i for i in adata.obs['Donor']])))
celltypes = sorted(list(set([i for i in adata.obs['Cluster']])))
colors = adata.uns['Cluster_colors']
cellcount = {}
for ct in celltypes:
    cellcount[ct] = []

for d in donors:
    for ct in cellcount:
        df = adata.obs[(adata.obs['Donor'] == d)]
        x = len(df[(df['Cluster'] == ct)].index)
        t = len(df.index)
        p = (x/t)*100
        cellcount[ct].append(p)
        
fig, ax = plt.subplots(dpi=300, figsize=(3,6))

labels = {ct: ct_label for ct, ct_label in zip(adata.obs['Cluster'], adata.obs['Cell_type'])}

left = None
for i, (ct, color) in enumerate(zip(cellcount, colors)):
    if i == 0:
        ax.barh(donors, cellcount[ct], label=labels[ct], color=color)
        left = np.array(cellcount[ct])
    else: 
        ax.barh(donors, cellcount[ct], left=left, color=color, label=labels[ct])
        left = left + np.array(cellcount[ct])

labels = []
for d in donors:
    if d.startswith('WAN'):
        labels.append(str("WAN-" + d.rsplit('_')[1]))
    if d.startswith('SER'):
        labels.append(str("SEG-" + d.rsplit('_')[1]))
    if d.startswith('XIN'):
        labels.append(str("XIN-" + d.rsplit('_')[1]))
    if d.startswith('ENG'):
        labels.append(str("ENG-" + d.rsplit('_')[1]))

labels = [l.replace('donor','').replace('Donor','') for l in labels]
        
ax.grid(False)
ax.set_yticklabels(labels, rotation=0)
for label in ax.get_yticklabels():
    label.set_fontsize(10)
ax.legend()
ax.legend(bbox_to_anchor=(1, 0.75))
ax.set_xlabel('Percentage of \ncell types (%)')
ax.set_ylabel('Donors')

fig.savefig(filenameout, bbox_inches='tight')