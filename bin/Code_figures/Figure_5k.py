#!/bin/python

###########################
##    eQTL Figure 5K     ##
###########################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import random

filenameout = 'figure5k.png'

df = pd.read_csv('/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/UBC-stanford/static_gsis_ND.csv')

days = [4,5,6]
concentration = '1_16.7'

TT = df[(df['Culture days'].isin(days)) & (df['geno'] == 'TT')][concentration]
TC = df[(df['Culture days'].isin(days)) & (df['geno'] == 'TC')][concentration]
CC = df[(df['Culture days'].isin(days)) & (df['geno'] == 'CC')][concentration]

fig, ax = plt.subplots(dpi=300, figsize=(2,3))

x = [1,2,3]
data = [TT, TC, CC]

print('n=',sum([len(TT), len(TC), len(CC)]))

y = [np.mean(d) for d in data]
yerr = [np.std(d)/np.sqrt(len(d)) for d in data]

ax.bar(x, y, yerr=yerr, color='#F3DD88', edgecolor='black')

for i, geno in enumerate([TT, TC, CC]):
    x = [i+1+random.uniform(-0.1,0.1) for j in geno]
    ax.scatter(x, geno, color='black',s=10, edgecolor='black',zorder=3)

ax.set_xticks([1,2,3], ['T/T', 'T/C', 'C/C'])

#     ax.set_ylim(0,6)
ax.set_xlim(0.5,3.5)
ax.set_ylabel('Stimulation index')

ax.set_title(f'Insulin secretion\nHumanIslets Dataset') #\nculture time: {days} days')
ax.grid(False)

x = []
y = []
for i,d in enumerate(data):
    x += len(d)*[i+1]
    y += list(d)
    
print(stats.linregress(list(x),list(y)))

slope, intercept, r_value, p_value, std_err = stats.linregress(list(x),list(y))

print(f'r2 = {r_value}')
print(f'p = {p_value}')

ys = [intercept+slope*xs for xs in [1,2,3]]
ax.plot([1,2,3], ys, color="black", lw=1.5)

print('TT vs TC', stats.ttest_ind(TT, TC))
print('TT vs CC',stats.ttest_ind(TT, CC))
print('TC vs CC',stats.ttest_ind(TC, CC))

fig.savefig(filenameout, bbox_inches='tight')
