import sys
import pandas as pd
import scanpy as sc

file = '../DEseq_matrix/eQTL_revision_filtered_count_matrix_normalized.csv'
cell_type = sys.argv[1] #'beta-cells'
group = sys.argv[2] #'ND'

h5ad = '../notebook/eQTL_revision.h5ad'
adata = sc.read_h5ad(h5ad)

DEseq_matrix = pd.read_csv(file, index_col=0)
genes = [gene.split('>')[0] for gene in DEseq_matrix.index]
cells = [cell.replace('.','-') for cell in DEseq_matrix.columns]

DEseq_matrix.index = genes
DEseq_matrix.columns = cells

print('dataframe loaded')

adata_sub = adata[(adata.obs['Cell_type'] == cell_type) & (adata.obs['Group'] == group)]
cells = list(adata_sub.obs.index)

DEseq_matrix_subset = DEseq_matrix.loc[:,cells]

DEseq_matrix_subset.index = genes

DEseq_matrix_subset = DEseq_matrix_subset[DEseq_matrix_subset.sum(axis=1)>0]

cell_type = cell_type.replace('/','-')

DEseq_matrix_subset.to_csv('../bin/{cell_type}_{group}_filtered_count_matrix_normalized.csv')