#!/bin/python

###########################
##  eQTL SUPP Figure 1E  ##
###########################

import scanpy as sc
import os

filenameout = 'supp_figure1e.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

color_all = '#0B3652'
color_T2D = '#7E7F9A'
color_healthy = '#F3DE8A'

adata.uns['Group_colors'] = [color_T2D, color_healthy]

sc.settings.figdir = './'

sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white')
sc.pl.umap(adata, color=['Group'],ncols=3,scale_factor=400, show=False, save=filenameout)

os.system(f'mv umap{filenameout} {filenameout}')