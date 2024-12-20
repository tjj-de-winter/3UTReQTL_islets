### Description ###

# Code to make a filtering matrix (PASS or FAIL) for each gene and each cell type
# Matrix is used to filter out eQTLs in cell types with low expression of the gene

### Import packages ###
import os
import pandas as pd
import numpy as np
import scanpy as sc
import argparse

### Input variables ### 
parser = argparse.ArgumentParser(description='Generate a gene filter matrix based on expression thresholds per cell type')
parser.add_argument('--h5ad',help='a scanpy generated h5ad object with cell type annotation')
parser.add_argument('--group_header', default="Group",help='Name of the column in adata.obs that contains group or condition information')
parser.add_argument('--celltype_header', default="Cell_type",help='Name of the column in adata.obs that contains cell type information')
parser.add_argument('--cell_and_group', action="store_true",help='Enable to compute gene filter for each cell type and group seperately, if not used only cell type is used')
parser.add_argument('--outprefix', default='', help='outprefix for output files')
parser.add_argument('--outdir', default='.', help='output path')

args = parser.parse_args()

h5ad = args.h5ad #'/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/1-scRNA-seq/3-h5ad/3UTRpancreas.h5ad'
group_col = args.group_header 
celltype_col = args.celltype_header
cell_and_group = args.cell_and_group
output_path = args.outdir
outprefix = args.outprefix

adata = sc.read_h5ad(h5ad)

### Parameters ###

median_filter = 0.25 # median gene expression
percentage_filter = 20 # percentage of cells expressing the gene

### Functions ###

def gene_expression_filtering(adata, celltype, gene):
    '''
    Based on expression parameters generates PASS or FAIL value indicating if gene is inluded or exluded for analysis, respectively.
    Function looks at median expression levels per cell type and percentage of cells that express the gene. 

    Args:
        adata: adata object
        celltype (str): cell type name
        gene (str): gene name

    Returns:
        str: PASS or FAIL
        '''
    adata_sub = adata[adata.obs[column_obs] == celltype]
    
    matrix = pd.DataFrame(adata_sub.raw.X, columns=adata_sub.raw.var.index, index=adata_sub.obs.index)
    
    matrix_gene = list(matrix[gene])
    
    percentage_expressed = np.count_nonzero(matrix_gene)/len(matrix_gene)*100
    mean_expressed = np.mean(matrix_gene)
    median_expressed = np.median(matrix_gene)
    
    if median_expressed < median_filter: # filter out genes with a certain median expression value
        return 'FAIL'
    elif percentage_expressed < percentage_filter: # filter out genes that are expressed in a certain percentage of cells
        return 'FAIL'
    else:
        return "PASS"

### Code ###

if not os.path.exists(output_path):
        os.makedirs(output_path)

# Combine cell and group annotation as a single string
if cell_and_group:
    column_obs = 'cell_group'
    adata.obs[column_obs] = ['_'.join([c,g]) for c, g in zip(adata.obs[celltype_col], adata.obs[group_col])]
else:
    column_obs = celltype_col

# perform the filtering and generate a CSV file
df_filter = pd.DataFrame(index=adata.raw.var.index, columns=list(set(adata.obs[column_obs])))

for gene in df_filter.index:
    for ct in df_filter.columns:
        df_filter.loc[gene, ct] = gene_expression_filtering(adata, ct, gene)

outfile = f'{output_path}/{outprefix}_celltype_gene_filter.csv'
df_filter.to_csv(outfile)

