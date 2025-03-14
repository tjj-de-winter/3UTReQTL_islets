#!/bin/python

###########################
##  eQTL SUPP Figure 3C  ##
###########################

import glob
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import numpy as np

filenameout = 'supp_figure3c.png'

files = glob.glob('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/mirbase/miranda.out.filtered.*cells*.csv')

df_dict = {}
for file in files:
    split = file.rsplit('.')
    name = '_'.join([split[3],split[4]])
    df = pd.read_csv(file, index_col=0)
    df.index = range(df.shape[0])
    df['mirna'] = [mirna.strip('>') for mirna in df['mirna']]
    df = df[(df['slope_binding'] > 0)&(df['slope_binding_pval'] < 0.05)]
    df['variant-miRNA'] = ['-'.join([v,m]) for v, m in zip(df['variant'], df['mirna'])]
    df['cell_group'] = df['Target'].apply(lambda x: name)
    df['group'] = df['cell_group'].apply(lambda x: x.rsplit('_')[1])
    df_dict[name] = df

df_merge = pd.concat(df_dict.values())

mirna_counts = {}
mirnas = []

for mirna in set(df_merge['mirna']):
    df_mirna = df_merge[df_merge['mirna'] == mirna]
    mirna_counts[mirna] = len(set(df_mirna['variant-miRNA']))
    mirnas += [mirna]*len(set(df_mirna['variant-miRNA']))
    
order = dict(Counter(mirnas).most_common()).keys()
mirna_counts = {mirna: mirna_counts[mirna] for mirna in order}

def determine_group(df):
    groups = list(set(df['group']))
    if len(groups) == 2:
        return 'shared'
    elif len(groups) == 1 and 'T2D' in groups:
        return 'T2D'
    else:
        return 'ND'


mirnas = []
for mirna in set(df_merge['mirna']):
    df_mirna = df_merge[df_merge['mirna'] == mirna]
    mirna_counts[mirna] = len(set(df_mirna['variant-miRNA']))
    mirnas += [mirna]*len(set(df_mirna['variant-miRNA']))    

order = dict(Counter(mirnas).most_common()).keys()
    
ND = np.zeros(len(set(df_merge['mirna'])))
T2D = np.zeros(len(set(df_merge['mirna'])))
shared = np.zeros(len(set(df_merge['mirna'])))

for idx, mirna in enumerate(list(order)):
    df_mirna = df_merge[df_merge['mirna'] == mirna]
    for v, dfv in df_mirna.groupby('variant'):
        group = determine_group(dfv)
        if group == 'ND':
            ND[idx] += 1
        elif group == 'T2D':
            T2D[idx] += 1
        elif group == 'shared':
            shared[idx] += 1
            
fig, ax = plt.subplots(dpi=300, figsize=(25,4))

x = range(len(order))

ax.bar(x, ND, color='#3BAA35', edgecolor='black', label='ND unique')
ax.bar(x, T2D, color='#55529F', edgecolor='black', bottom=ND, label='T2D unique')
ax.bar(x, shared, color='#E0A604', edgecolor='black', bottom=ND+T2D, label='shared')

ax.legend(fontsize=20)
ax.set_yticks([0,10,20,30,40,50,60],[0,10,20,30,40,50,60],fontsize=20)
ax.set_ylabel('number of eQTLs', fontsize=20)

def get_xticks(mirna_list, order):
    x = []
    label = []
    for idx, mirna in enumerate(order):
        if mirna in mirna_list:
            x.append(idx)
            label.append(mirna)
    return x, label

mirna_list = ['hsa-miR-143-5p', 'hsa-miR-127-5p', 'hsa-miR-654-5p', 'hsa-miR-186-5p', 'hsa-miR-543',
             'hsa-miR-299-5p', 'hsa-miR-150-5p', 'hsa-miR-4443', 'hsa-miR-151a-3p', 'hsa-miR-323b-3p',
             'hsa-miR-1468-5p', 'hsa-miR-30c-2-3p', 'hsa-miR-224-5p', 'hsa-miR-149-5p', 'hsa-miR-23a-5p',
             'hsa-miR-29b-1-5p', 'hsa-miR-30b-3p']
xticks, xlabels = get_xticks(mirna_list, order)

ax.set_xticks(xticks, xlabels, rotation=90);
ax.set_xlim(-1,len(order))

fig.savefig(filenameout, bbox_inches='tight')