#!/bin/python

###########################
##  eQTL SUPP Figure 1F  ##
###########################

import scanpy as sc
import os

filenameout = 'supp_figure1f.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

sc.settings.figdir = './'

sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white')
sc.pl.umap(adata, color=['GHRL'],ncols=3,scale_factor=400, color_map='Reds', size=100, show=False, save=filenameout)

os.system(f'mv umap{filenameout} {filenameout}')