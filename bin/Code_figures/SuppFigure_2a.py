#!/bin/python

###########################
##  eQTL SUPP Figure 2A  ##
###########################

import pandas as pd
import matplotlib.pyplot as plt

filenameout1 = 'supp_figure2a_1.png'
filenameout2 = 'supp_figure2a_2.png'

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_healthy.adjusted-pvalue.csv'
df = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df = df[df['pvalue'] < 0.05].drop_duplicates()


neQTLs = {}
for cell_type in set(df.type):
    df_ct = df[df['type'] == cell_type]
    df_rest = df[df['type'] != cell_type]
    
    ct = set(df_ct['variant'])
    rest = set(df_rest['variant'])
    
    unique = len(ct.difference(rest))
    shared = len(ct.intersection(rest))
    
    neQTLs[cell_type.rsplit('_')[0]] = [unique, shared]

colors_cells = {'acinar-cells': '#1f77b4', 'ductal-cells':'#ff7f0e', 'endothelial-cells':'#279e68', 'PP-cells':'#d62728','stellate-cells': '#aa40fc', 'alpha-cells':'#8c564b', 'beta-cells':'#e377c2', 'delta-cells':'#b5bd61'}
colors_cells_light = {'acinar-cells': '#B5CDE4', 'ductal-cells':'#F5D3AC', 'endothelial-cells':'#B9DAC8', 'PP-cells':'#E8B4B4','stellate-cells': '#CFBBDB', 'alpha-cells':'#D3C4C0', 'beta-cells':'#ECD0E4', 'delta-cells':'#E4E7C7'}


fig, ax = plt.subplots(dpi=300, figsize=(5,5))

cell_order = ['beta-cells', 'alpha-cells', 'delta-cells', 'PP-cells', 'acinar-cells', 'ductal-cells', 'stellate-cells', 'endothelial-cells']

for x, key in enumerate(cell_order):
    d = neQTLs[key]
    total = sum(d)
    percentage = [n/total*100 for n in d]
    print(key, percentage)
    
    ax.bar(x, percentage[0], color=colors_cells[key], edgecolor='black')
    ax.bar(x, percentage[1], bottom=percentage[0], color=colors_cells_light[key], edgecolor='black')

# make a plot with black and gray for the labels
ax.bar(x+2, percentage[0], color='black', alpha=.3, label='Shared eQTLs')
ax.bar(x+2, percentage[1], bottom=percentage[0],color='black', label='Unique eQTLs')
    
ax.grid(False)
labels = [key for key in cell_order]
ax.set_xlim(-1,len(labels))
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=90)
ax.set_ylabel('eQTL percentage (%)')
ax.set_title('ND')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

fig.savefig(filenameout1, bbox_inches='tight')

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_T2D.adjusted-pvalue.csv'
df = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df = df[df['pvalue'] < 0.05].drop_duplicates()

neQTLs = {}
for cell_type in set(df.type):
    df_ct = df[df['type'] == cell_type]
    df_rest = df[df['type'] != cell_type]
    
    ct = set(df_ct['variant'])
    rest = set(df_rest['variant'])
    
    unique = len(ct.difference(rest))
    shared = len(ct.intersection(rest))
    
    neQTLs[cell_type.rsplit('_')[0]] = [unique, shared]

colors_cells = {'acinar-cells': '#1f77b4', 'ductal-cells':'#ff7f0e', 'endothelial-cells':'#279e68', 'PP-cells':'#d62728','stellate-cells': '#aa40fc', 'alpha-cells':'#8c564b', 'beta-cells':'#e377c2', 'delta-cells':'#b5bd61'}
colors_cells_light = {'acinar-cells': '#B5CDE4', 'ductal-cells':'#F5D3AC', 'endothelial-cells':'#B9DAC8', 'PP-cells':'#E8B4B4','stellate-cells': '#CFBBDB', 'alpha-cells':'#D3C4C0', 'beta-cells':'#ECD0E4', 'delta-cells':'#E4E7C7'}



fig, ax = plt.subplots(dpi=300, figsize=(5,5))

cell_order = ['beta-cells', 'alpha-cells', 'delta-cells', 'PP-cells', 'acinar-cells', 'ductal-cells', 'stellate-cells', 'endothelial-cells']


for x, key in enumerate(cell_order):
    d = neQTLs[key]
    total = sum(d)
    percentage = [n/total*100 for n in d]
    
    print(key, percentage)
    
    ax.bar(x, percentage[0], color=colors_cells[key], edgecolor='black')
    ax.bar(x, percentage[1], bottom=percentage[0], edgecolor='black', color=colors_cells_light[key])
    

# make a plot with black and gray for the labels
ax.bar(x+2, percentage[0], color='black', edgecolor='black',alpha=.3, label='Shared eQTLs')
ax.bar(x+2, percentage[1], bottom=percentage[0],color='black', edgecolor='black', label='Unique eQTLs')
    
ax.grid(False)
labels = [key for key in cell_order]
ax.set_xlim(-1,len(labels))
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=90)
ax.set_ylabel('eQTL percentage (%)')
ax.set_title('T2D')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

fig.savefig(filenameout2, bbox_inches='tight')