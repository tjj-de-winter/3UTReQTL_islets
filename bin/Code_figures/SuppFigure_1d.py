#!/bin/python

###########################
##  eQTL SUPP Figure 1D  ##
###########################

import scanpy as sc
import os

filenameout = 'supp_figure1d.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

sc.settings.figdir = './'

sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white')

dataset_dict = {'SER':'Segerstolpe et al. 2016', 'WAN':'Wang et al. 2016', 'XIN16':'Xin et al. 2016',
               'ENG':'Enge et al. 2017'}
adata.obs['Dataset_label'] = adata.obs['Dataset'].apply(lambda x: dataset_dict[x])

sc.pl.umap(adata, color=['Dataset_label'],ncols=3,scale_factor=400, title='Dataset', show=False, save=filenameout)

os.system(f'mv umap{filenameout} {filenameout}')