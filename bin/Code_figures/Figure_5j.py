#!/bin/python

###########################
##    eQTL Figure 5J     ##
###########################

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import random

filenameout = 'figure5j.png'

days = 4
    
TT = [5.23, 3.56, 2.04]
TC = [2.35, 1.91, 2.19, 1.50, 0.89]
CC = [0.46, 0.77]

fig, ax = plt.subplots(dpi=300, figsize=(2,3))

x = [1,2,3]
data = [TT, TC, CC]
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

ax.set_title(f'Insulin secretion')
ax.grid(False)
    
x = []
y = []
for i,d in enumerate(data):
    x += len(d)*[i+1]
    y += d
    
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

