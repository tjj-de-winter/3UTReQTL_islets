import sys
import pandas as pd
import scanpy as sc
import numpy as np
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=Warning)

cell_type = sys.argv[1] # 'beta-cells'
group = sys.argv[2] # 'ND'
celltype = '_'.join([cell_type, group])
VCFfile = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/VCF/15.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf'
h5ad_file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/notebook/eQTL_revision.h5ad'
n_PCs_gene = 10
n_PCs_geno = 10

# h5ad
adata = sc.read_h5ad(h5ad_file)

print('adata loaded')

meta = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/eQTL_pipeline/metadata/3UTRpancreas_donor_metadata.v2.csv'
metadata = pd.read_csv(meta, sep=",")
metadata['donor_ID'] = ['_'.join([lib, donor]) for lib, donor in zip(metadata['library_ID'], metadata['donor_ID'])]
metadata.index = metadata['donor_ID']
adata.obs['Sex'] = ['Male' if metadata['Sex'][i] in 'Male' else 'Female' for i in adata.obs['Donor']]
adata.obs['Sex'] = adata.obs['Sex'].astype('category')

# prepare adata objects
adata = adata.raw.to_adata() # to get all genes, use adata.raw

adata_sub = adata[(adata.obs['Cell_type'] == cell_type)& (adata.obs['Group'] == group)]
adata_sub.obs['log_total_counts'] = np.log(adata_sub.obs['total_counts'])
donors = list(set(adata_sub.obs['DonorID']))
matrix_norm = pd.DataFrame(adata_sub.X, index=adata_sub.obs.index, columns=adata_sub.var.index)

# PCs gene expression
PCA_name='X_pca_harmony'

pc_df = pd.DataFrame(adata_sub.obsm[PCA_name], index=adata_sub.obs.index)
pc_df.columns = [f'genePC{i+1}' for i in pc_df]

columns = [f'genePC{i+1}' for i in range(n_PCs_gene)]
pc_df = pc_df.loc[:,columns]

index = adata_sub.obs.index
for pcn in columns:
    adata_sub.obs[pcn] = list(pc_df.loc[index][pcn])

print('PCs of gene expression added')

cell_type = cell_type.replace('/','-')
adata_sub.obs.to_csv(f'{cell_type}_{group}.covariates.csv')