#!/bin/python

###########################
##  eQTL SUPP Figure 1K  ##
###########################

import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats


filenameout = 'supp_figure1k.png'

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv'
df_eqtl = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df_eqtl['gene'] = df_eqtl['gene'].apply(lambda x: x.rsplit('>')[0])

df_eqtl = df_eqtl.drop_duplicates()

df_ND = df_eqtl[df_eqtl['type'].isin(['beta-cells_healthy'])]
df_T2D = df_eqtl[df_eqtl['type'].isin(['beta-cells_T2D'])]

def count(df):
    n_sig = df[df['pvalue']<=0.05].shape[0]
    n_nonsig = df[df['pvalue']>=0.05].shape[0]
    total = df.shape[0]
    
    return n_sig, n_nonsig, total

fontsize = 14

plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

ND = count(df_ND)

T2D = count(df_T2D)

data_sig = [ND[0]/ND[2]*100, T2D[0]/T2D[2]*100]
data_nonsig = [ND[1]/ND[2]*100, T2D[1]/T2D[2]*100]

# Fisher's exact test

contingency_table = [[ND[0], ND[1]], 
                     [T2D[0], T2D[1]]]

odds_ratio, p_value = stats.fisher_exact(contingency_table)

print(odds_ratio, p_value)


fig, ax = plt.subplots(dpi=300, figsize=(5,5))

x = [0,1]
ax.bar(x, data_sig,align='edge', width=1,label='eQTLs', edgecolor='black', color='#61bac9')
ax.bar(x, data_nonsig, align='edge', width=1,bottom=data_sig, label='variants', edgecolor='black', color='#8a465a')

ax.set_ylabel('%')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fontsize)
ax.set_xticks([0.5,1.5])
ax.set_xticklabels(['ND beta-cells', 'T2D beta-cells'])

ax.grid(False)

ax.set_ylim(0,115)

fig.savefig(filenameout, bbox_inches='tight')