#!/bin/python

###########################
##    eQTL Figure 1B     ##
###########################

import matplotlib.pyplot as plt
import sys,os
import scanpy as sc

#
filenameout1 = 'figure1b_1.png'
filenameout2 = 'figure1b_2.png'
filenameout3 = 'figure1b_3.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

# plot the data part 1

sc.set_figure_params(scanpy=True, fontsize=30, dpi=300, facecolor='white')
# Create a figure with custom spacing
fig, axs = plt.subplots(1, 3, figsize=(17, 5))

# Adjust the space between subplots
fig.subplots_adjust(wspace=0.3)

# Plot UMAP with different colors in subplots
sc.pl.umap(adata, color='INS', ax=axs[0], show=False, scale_factor=300)
sc.pl.umap(adata, color='GCG', ax=axs[1], show=False, scale_factor=300)
sc.pl.umap(adata, color='SST', ax=axs[2], show=False, scale_factor=300)

# Adjust colorbar ticks size for each plot
for ax in axs:
    cbar = ax.collections[0].colorbar  
    cbar.ax.tick_params(labelsize=20)  

fig.savefig(filenameout1, bbox_inches='tight')

# plot the data part 2

sc.set_figure_params(scanpy=True, fontsize=30, dpi=300, facecolor='white')
# Create a figure with custom spacing
fig, axs = plt.subplots(1, 3, figsize=(17, 5))

# Adjust the space between subplots
fig.subplots_adjust(wspace=0.3)

# Plot UMAP with different colors in subplots
sc.pl.umap(adata, color='PPY', ax=axs[0], show=False, scale_factor=300)
sc.pl.umap(adata, color='CPA1', ax=axs[1], show=False, scale_factor=300)
sc.pl.umap(adata, color='KRT19', ax=axs[2], show=False, scale_factor=300)

# Adjust colorbar ticks size for each plot
for ax in axs:
    cbar = ax.collections[0].colorbar  
    cbar.ax.tick_params(labelsize=20) 

fig.savefig(filenameout2, bbox_inches='tight')

# plot the data part 3

sc.set_figure_params(scanpy=True, fontsize=30, dpi=300, facecolor='white')
# Create a figure with custom spacing
fig, axs = plt.subplots(1, 2, figsize=(11.3, 5))

# Adjust the space between subplots
fig.subplots_adjust(wspace=0.3)

# Plot UMAP with different colors in subplots
sc.pl.umap(adata, color='COL1A1', ax=axs[0], show=False, scale_factor=300)
sc.pl.umap(adata, color='PECAM1', ax=axs[1], show=False, scale_factor=300)

# Adjust colorbar ticks size for each plot
for ax in axs:
    cbar = ax.collections[0].colorbar  
    cbar.ax.tick_params(labelsize=20)

fig.savefig(filenameout3, bbox_inches='tight')
