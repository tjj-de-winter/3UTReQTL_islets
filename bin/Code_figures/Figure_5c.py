#!/bin/python

###########################
##    eQTL Figure 5C     ##
###########################

import matplotlib.pyplot as plt
import scanpy as sc
import gseapy
import numpy as np
import os

filenameout = 'figure5c.png'

# load scRNAseq data

h5ad = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/RESULTS/ANALYSIS/3UTRpancreas.h5ad'
adata = sc.read_h5ad(h5ad)

def GO(gene_list, background='all', geneset='Reactome_2022', title='', top=10, figsize=(10,5), color='red'):
    
    if background == 'all':
        background = list(adata.raw.var.index)
    
    enr = gseapy.enrichr(gene_list=gene_list,
                         organism='Human',
                         background=background ,
                         gene_sets=geneset,
                         description='pathway',
                         cutoff = 0.1)
    enr_df = enr.results
    
    df_sig = enr_df[enr_df['Adjusted P-value']< 0.05]
    
    df_sig = df_sig.iloc[:top,:]
    
    # plot
    data = list(-np.log10(df_sig['Adjusted P-value']))
    labels = list(df_sig['Term'])

    data.reverse()
    labels.reverse()
    
    fig, ax = plt.subplots(dpi=300, figsize=figsize)
    ax.barh(range(len(df_sig.index)), data, align='center', color=color, edgecolor='black')
    ax.set_yticks(range(len(df_sig.index)))
    ax.set_yticklabels(labels)
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('-log10(Adjusted P-value)')

    ax.set_title(f'{title}\n{geneset}', fontweight="bold");
    
    return fig

up = ['MTDH', 'TMEM212', 'ATP5F1E', 'CSNK1E', 'KCNJ3', 'MT1X', 'GIGYF1', 'POGZ', 'SPTBN1', 'TOMM7', 'AHCYL1', 'INSR', 'FXYD2', 'EIF3A', 'COX7B', 'NDUFV3', 'H4C3', 'PRELID3B', 'LRPPRC', 'NDUFB4', 'RPL38', 'ATRX', 'FAM76A', 'USP2', 'RBBP5', 'RPS19', 'NKX3-1', 'AC004805.1', 'ZNF516', 'ZBTB20', 'ZFX', 'LARP1', 'ARF6', 'RPS8', 'SYT14', 'PCBP2', 'RAP2B', 'ACTR10', 'MIDN', 'CRIPT', 'CYP4V2', 'ZMYND8', 'SMIM7', 'LIN7C', 'NDUFA1', 'WDR43', 'LSM12', 'MT-TS2', 'SLC29A4', 'AC119674.2', 'OGT', 'ANKRD40', 'CAMK2N1', 'PAN3', 'RPL34', 'NPTN', 'IRF2BPL', 'EEF2', 'TLK1', 'EAF1', 'PTBP2', 'COX6C', 'EXPH5', 'CACNA1A', 'RNASEK', 'ATP5ME', 'IFI6', 'VEGFA', 'RPS14', 'MT-TG', 'CDC42', 'DLG5', 'EIF4G3', 'TBC1D7', 'ATXN7', 'RASGEF1B', 'NDUFB3', 'DDX17', 'BRD4', 'SKIL', 'ANP32E', 'RPL18', 'TSPYL1', 'TIAL1', 'CDKN1C', 'TULP4', 'KIAA1109', 'YWHAE', 'PABPC1', 'CTTN', 'FXR1', 'INS', 'RRAGD', 'THAP5', 'MAP1B', 'RPL36AL', 'B2M', 'ALYREF', 'PCMTD1', 'EFNA5', 'PFDN5', 'MT-TL1', 'SST', 'CTBP1', 'CLK1', 'BTBD9', 'DNMT3A', 'MZT2A', 'RPS27A', 'MLXIPL', 'TBL1XR1', 'LUC7L3', 'CSNK1A1', 'RPL39', 'LDB1', 'RSRC1', 'RUFY3', 'RPL7A', 'BAX', 'PDCD4', 'EIF1AD', 'PER2', 'RPLP1', 'CCDC58', 'EGR1', 'PPP3CA', 'CBX6', 'XIAP', 'SRPRB', 'MT-TE', 'RPL36', 'SSR4', 'CNOT9', 'WSB1', 'MPP6', 'TNRC6B', 'RPS21', 'C6orf62', 'PRPF39', 'PHPT1', 'CALR', 'RPS11', 'SAR1A', 'ERO1A', 'MORF4L1', 'FNIP1', 'STAT2', 'MARCKS', 'FOXP1', 'UBE2B', 'MZT2B', 'HSP90AA1', 'THRAP3', 'STRN3', 'SRRM2', 'GOLIM4', 'TMED4', 'ADCYAP1', 'CAPZA1', 'PDGFA', 'TMEM33', 'TRA2A', 'JUND', 'ARID1B', 'ZDHHC2', 'RAB12', 'MIDEAS', 'SET', 'DNAJC3', 'ARGLU1', 'EXOSC6', 'WTAP', 'PTEN', 'MARCHF6', 'EIF2S3', 'MDM4', 'CTNNB1', 'PPP1CB', 'HNRNPL', 'HNRNPH1', 'PDS5B']

fig = GO(up, background='all', geneset='KEGG_2021_Human', title='Genes upregulated in C/C beta cells', top=20, figsize=(7,5), color='#B5C228')

fig.savefig(filenameout, bbox_inches='tight')

os.system('rm -r Enrichr')