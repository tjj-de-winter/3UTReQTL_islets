#!/bin/python

###########################
##  eQTL SUPP Figure 1A  ##
###########################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob

filenameout = 'supp_figure1a.png'

filter_file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/1-scRNA-seq/1-count_tables/filter_criteria_3utrcount_smart-seq.csv'

count_files = glob.glob('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/1-scRNA-seq/1-count_tables/*.tsv')

# Make dataframe from count table and append to dictionary
dfs = {}

for file in count_files:
    name = file.rsplit('/')[-1].rsplit('.')[0]
    df = pd.read_csv(file, sep='\t', index_col=0)
    dfs[name] = df
    
# Make dataframe from filtering criteria

df_filter = pd.read_csv(filter_file, sep=';', index_col=0)
df_filter.index = [i.upper() for i in df_filter.index]

def read_counts_per_cell(df):
    read_counts = list(df.sum())
    return np.array(read_counts)

def mito_ratio_per_cell(df):
    MT_index = [i for i in df.index if i.startswith('MT-')]
    df_mito = df.loc[MT_index,:]
    mito_counts = list(df_mito.sum())
    all_counts = list(df.sum())
    mito_ratio = np.array(mito_counts)/np.array(all_counts)
    return mito_ratio

def get_treshold(df, name, col):
    name = name.upper()
    return df.loc[name,col]

fig, (ax1,ax2) = plt.subplots(2, sharex=True, dpi=300)

### top plot

data_reads = [] 
labels = []
names = sorted([name for name in dfs.keys()])
for name in names:
    df = dfs[name]
    labels.append(name)
    data_reads.append(read_counts_per_cell(df))
    
labels = [l.rsplit('_count')[0].replace('donor','').replace('Donor','') for l in labels]    
    
x = range(len(labels))
violin = ax1.violinplot(data_reads, 
                       positions=x, 
                       widths=1, 
                       showmeans=False,
                       showextrema=False)

# Customize the appearance
for pc in violin['bodies']:
    pc.set_facecolor('#3252a8')  # Change fill color
    pc.set_edgecolor('black')      # Change edge color
    pc.set_alpha(1)             # Set transparency

for xi, name in zip(x,names):
    name = name.rsplit('_count-matrix')[0]
    name = name.upper()
    value = get_treshold(df_filter, name, 'min_reads')
    if xi == 0:
        ax1.hlines(y=value, xmin=xi-0.4,xmax=xi+0.4,color='red', zorder=2, linestyles='solid', label='threshold')
    else:
        ax1.hlines(y=value, xmin=xi-0.4,xmax=xi+0.4,color='red', zorder=2, linestyles='solid')

ax1.set_yscale('log')
ax1.legend(loc='lower right')
ax1.set_ylabel('reads per cell\n(prefilter)')

### bottom plot

data_reads = [] 
labels = []
names = sorted([name for name in dfs.keys()])
for name in names:
    df = dfs[name]
    labels.append(name)
    data_reads.append(mito_ratio_per_cell(df))
    
labels = [l.rsplit('_count')[0].replace('donor','').replace('Donor','') for l in labels]    
    
x = range(len(labels))
violin = ax2.violinplot(data_reads, 
                       positions=x, 
                       widths=1, 
                       showmeans=False,
                       showextrema=False)

# Customize the appearance
for pc in violin['bodies']:
    pc.set_facecolor('#3252a8')  # Change fill color
    pc.set_edgecolor('black')      # Change edge color
    pc.set_alpha(1)             # Set transparency

for xi, name in zip(x,names):
    name = name.rsplit('_count-matrix')[0]
    name = name.upper()
    value = get_treshold(df_filter, name, 'mito_frac')
    if xi == 0:
        ax2.hlines(y=value, xmin=xi-0.4,xmax=xi+0.4,color='red', zorder=2, linestyles='solid', label='threshold')
    else:
        ax2.hlines(y=value, xmin=xi-0.4,xmax=xi+0.4,color='red', zorder=2, linestyles='solid')
    
ax2.set_xticks(range(len(labels)), labels, rotation=90, zorder=0, size=8)
# ax.set_yscale('log')
ax2.legend()
ax2.set_ylabel('Mito gene fraction per cell\n(prefilter)')

fig.savefig(filenameout, bbox_inches='tight')

