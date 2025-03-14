#!/bin/python

###########################
##  eQTL SUPP Figure 1C  ##
###########################

import matplotlib.pyplot as plt
import scanpy as sc

filenameout = 'supp_figure1c.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

# total cell count per donor

donors = sorted(list(set(list(adata.obs['Donor']))))
colors = adata.uns['Donor_colors']

cellcount = []

for d in donors:
    df = adata.obs[(adata.obs['Donor'] == d)]
    t = len(df.index)
    cellcount.append(t)
        
fig, ax = plt.subplots(figsize=(10,5), dpi=300)

ax.bar(donors, cellcount,color='#3252a8')


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
        
labels = [l.replace('Donor','').replace('donor','') for l in labels]

ax.grid(False)
ax.set_xticklabels(labels, rotation=90, ha='center', size=15)
ax.set_ylabel('Total cell count\n(post filter)')
ax.set_xlabel('Donors')

fig.savefig(filenameout, bbox_inches='tight')