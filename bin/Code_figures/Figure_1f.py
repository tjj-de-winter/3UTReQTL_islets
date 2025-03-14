#!/bin/python

###########################
##    eQTL Figure 1F     ##
###########################

import glob
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt

filenameout = 'figure1f.png'

gwas_files =  glob.glob('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/5-public_data/GWAS/GWAS_SNPS/*.tsv')

gwas_dict = {}
for file in gwas_files:
    df = pd.read_csv(file, sep="\t")
    
    name = file.rsplit('/')[-1].rsplit('.EFO')[0].rsplit('.MONDO')[0].replace('Trait.', '').replace('.',' ').capitalize()
    df['Trait'] = [name for idx in df.index]
    
    gwas_dict[file.rsplit('.tsv')[0]] = df
gwas_dict

df_gwas = pd.concat([gwas_dict[key] for key in gwas_dict])
df_gwas.index = range(df_gwas.shape[0])

file="/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv"

df_eqtl_raw = pd.read_csv(file, sep=",", names=['gene', 'variant', 'pval', 'group'])
df_eqtl_raw['locations'] = df_eqtl_raw['variant'].apply(lambda v: v.rsplit('_')[0])
df_eqtl = pd.read_csv(file, sep=",", names=['gene', 'variant', 'pval', 'group'])
df_eqtl = df_eqtl.drop_duplicates()
df_eqtl = df_eqtl[df_eqtl['pval'] < 0.05]
df_eqtl['gene'] = df_eqtl['gene'].apply(lambda g: g.rsplit('>')[0])
df_eqtl['locations'] = df_eqtl['variant'].apply(lambda v: v.rsplit('_')[0])

group = 'beta-cells_healthy'
gene_chr = 'SLC30A8>8'
gene = 'SLC30A8'

dfx = df_eqtl_raw[(df_eqtl_raw['gene'] == gene_chr) & (df_eqtl_raw['group'] == group)]

pvals_eqtl = {}
pvals_gwas = {}
for l, p in zip(dfx.locations, dfx.pval):
    if p == 'None':
        p = 1
    p = float(p)
    
    if math.isnan(p):
        p = 1
    
    pvals_eqtl[l] = -np.log10(p)
    
    dfg = df_gwas[df_gwas['locations'] == l]
    
    pvals_gwas[l] = [],[]
    
    if len(dfg.index) == 0:
        pvals_gwas[l][0].append(-np.log10(1))
        pvals_gwas[l][1].append('None')
    
    for pg, t in zip(dfg.pValue, dfg.Trait):
        pg = -np.log10(pg)
        
        pvals_gwas[l][0].append(pg)
        pvals_gwas[l][1].append(t)

# plot

fontsize = 30
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

fig, (ax1, ax2) = plt.subplots(2, sharex=True, dpi=300, figsize=(10,10))
# fig, (ax1, ax2) = plt.subplots(2, sharex=True, dpi=300, figsize=(10,15))


locations = [int(l.rsplit(':')[1]) for l in dfx.locations]
xs = [xn-locations[0] for xn in locations]

label1 = False
label2 = False
dotsize = 80
for x, l in zip(xs, pvals_gwas): 
    lists = pvals_gwas[l]
    for p, t in zip(lists[0], lists[1]):
            if p > 1.3:
                if not label1:
                    ax1.scatter(x, p, color='red',s=dotsize , edgecolors='black',label='p <= 0.05')
                    label1 = True
                else:
                    ax1.scatter(x, p, color='red',s=dotsize , edgecolors='black')
            else:
                if not label2:
                    ax1.scatter(x, p, color='gray',edgecolors='black',s=dotsize , label='p >= 0.05')
                    label2 = True
                else:
                    ax1.scatter(x, p, color='gray',edgecolors='black',s=dotsize )
                
ax1.set_ylim(-2, 50)
ax1.text(-180, 30, 'rs3802177', size=20)
ax1.text(400, 40, 'rs11558471', size=20)
                
ax1.set_title("ND beta cells GWAS-eQTL colocalization\nSLC30A8 3'UTR")
ax1.legend(loc='upper right', fontsize=30)

label1 = False
label2 = False
for x, l in zip(xs, pvals_gwas):
    p = pvals_eqtl[l]
    if p > 1.3: 
        if not label1:
            ax2.scatter(x, p, color='#D279AF',s=dotsize , edgecolors='black',label='p <= 0.05')
            label1 = True
        else:
            ax2.scatter(x, p, color='#D279AF',s=dotsize , edgecolors='black')
    else:
        if not label2:
            ax2.scatter(x, p, color='gray',s=dotsize , edgecolors='black', label='p >= 0.05')
            label2 = True
        else:
            ax2.scatter(x, p,s=dotsize , color='gray', edgecolors='black')

ax1.set_ylabel('-log10(p-value) \nGWAS')
ax2.set_ylabel('-log10(p-value) \nbeta-cells ND eQTL')
ax2.legend(loc='upper right', fontsize=30)
ax1.grid(False)
ax2.grid(False)

ax2.set_xticks(xs)
ax2.set_xticklabels(['']*len(xs), rotation=90);

# ax2.set_xticklabels(list(dfx.locations), rotation=90);

print(list(dfx.locations))
plt.subplots_adjust(hspace=0.03)
ax1.tick_params(axis='x', which='both', length=7, width=2)
ax2.tick_params(axis='x', which='both', length=7, width=2)

fig.savefig(filenameout, bbox_inches='tight')
