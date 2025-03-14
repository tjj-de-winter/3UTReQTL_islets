#!/bin/python

###########################
##  eQTL SUPP Figure 1I  ##
###########################

from sklearn.linear_model import LinearRegression
from collections import Counter
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np

filenameout = 'supp_figure1i.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv'
df_eqtl = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df_eqtl['gene'] = df_eqtl['gene'].apply(lambda x: x.rsplit('>')[0])

df_eqtl = df_eqtl.drop_duplicates()

def ExtractColor(adata,obsKey='louvain',keytype=int):
#   labels=sorted(adata.obs[obsKey].unique().to_list(),key=keytype)
    labels=adata.obs[obsKey].cat.categories
    colors=adata.uns[obsKey+'_colors']
    return dict(zip(labels,colors))

eQTL_counter = Counter(df_eqtl[df_eqtl['pvalue'] < 0.05]['type'])

adata.obs['cell-group'] = ['_'.join([c,g]) for c,g in zip(adata.obs['Cell_type'], adata.obs['Group'])]

celltype_counter = Counter(adata.obs['cell-group'])

colors = ExtractColor(adata,obsKey='Cluster',keytype=int)

labels = {ct: ct_new for ct, ct_new in zip(adata.obs['Cell_type'], adata.obs['Cluster'])}

fig, ax = plt.subplots(dpi=150, figsize=(5,5))

x = []
y = []
for ct in eQTL_counter:
    xi = eQTL_counter[ct]
    x.append(xi)
    yi = celltype_counter[ct]
    y.append(yi)
    
    label = labels[ct.rsplit('_')[0]]
    if 'T2D' in ct:
        ax.scatter(xi, yi, color=colors[label],edgecolor='black', label=ct)
    else:
        ax.scatter(xi, yi, color=colors[label],label=ct.replace('healthy','ND'))
            
num = len(eQTL_counter.values())
start = max(eQTL_counter.values())
stop = min(celltype_counter.values())
xseq = np.linspace(start, stop, num=num)

x = np.array(x).reshape((-1,1))

model = LinearRegression()
model.fit(x, y)
r2 = model.score(x,y)
coefficients = model.coef_
intercept = model.intercept_

ax.plot(xseq, intercept+coefficients[0]*xseq, color="black", lw=1.5)
ax.text(1000,1000, f'$r^2$ = {round(r2,3)}')
ax.set_xlabel('Number of eQTLs')
ax.set_ylabel('Number of cells')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
ax.grid(False)

fig.savefig(filenameout, bbox_inches='tight')