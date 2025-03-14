#!/bin/python

###########################
##  eQTL SUPP Figure 1H  ##
###########################

import matplotlib.pyplot as plt
from upsetplot import plot
from upsetplot import from_memberships
import pandas as pd
import numpy as np
import itertools

filenameout = 'supp_figure1h.png'

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/3-eQTL/2-output_CSV/all-cells_allgroups.adjusted-pvalue.csv'
df_eqtl = pd.read_csv(file, sep=',', names=['gene', 'variant','pvalue', 'type'])
df_eqtl['gene'] = df_eqtl['gene'].apply(lambda x: x.rsplit('>')[0])

df_eqtl = df_eqtl.drop_duplicates()

df_eqtl_sig = df_eqtl[df_eqtl['pvalue'] < 0.05]

cell_types = set(df_eqtl_sig['type'])

ct_T2D = [ct for ct in cell_types if 'T2D' in ct]

def get_combinations(cell_types):
    combinations = []
    for n in range(len(cell_types)):
        n = n+1
        combinations += list(itertools.combinations(cell_types, n))
    return combinations

def count_variants(df, combination):
    
    matching_variants = []
    if len(combination) > 1:
        ct1 = set(df[df['type'] == combination[0]]['variant'])
        ct2 = set(df[df['type'] == combination[1]]['variant'])
        matching_variants = list(ct1.intersection(ct2))
    
    if len(combination) > 2:
        for ct in combination[2:]:
            ct = set(df[df['type'] == ct]['variant'])
            matching_variants = list(set(matching_variants).intersection(ct))
            
    if len(combination) == 1:
        matching_variants = list(set(df[df['type'] == combination[0]]['variant']))
        
    anti_combination = set(df['type']).difference(set(combination))
    
    df_rest = df[df['type'].isin(anti_combination)]
    
    variants = list(set(matching_variants).difference(set(df_rest['variant'])))
            
    return len(variants)
    
combinations = get_combinations(ct_T2D)

data = []
combinations_list = []
for com in combinations:
    count = count_variants(df_eqtl_sig, list(com))
    if count > 0:
        combinations_list.append(com)
        data.append(count)


plotdata = from_memberships(combinations_list,data=data)

plot(plotdata, facecolor='#3252a8', shading_color='lightgray')
plt.ylabel('eQTL frequency')
plt.yticks(np.linspace(0, round(max(data), -3), num=5))

plt.savefig(filenameout, bbox_inches='tight')


