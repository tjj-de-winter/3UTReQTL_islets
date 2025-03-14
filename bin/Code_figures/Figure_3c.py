#!/bin/python

###########################
##    eQTL Figure 3C     ##
###########################

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

filenameout = 'figure3C.png'

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/deseq2_norm_data_AA.tsv'
df_ct = pd.read_csv(file, index_col=0, sep="\t")
df_ct.columns = [col.rsplit('.')[0] for col in df_ct]

df_ct_norm = np.log(df_ct+1)
log_df = df_ct_norm

nddf = log_df[[c for c in log_df.columns if 'ND' in c]]
t2ddf = log_df[[c for c in log_df.columns if 'T2D' in c]]

def bootstrap_variance(v):
    bootstrap_var = np.array([np.random.choice(v, size = len(v)).var() for i in range(100)])
    return pd.Series({'BSvar': bootstrap_var.mean(), 'BSvar_err': bootstrap_var.std()/np.sqrt(len(v))})

var_nndf = pd.DataFrame(log_df[nddf.columns].apply(lambda x: bootstrap_variance(x), axis = 1))
var_t2ddf = pd.DataFrame(log_df[t2ddf.columns].apply(lambda x: bootstrap_variance(x), axis = 1))
var_bs = var_nndf.merge(var_t2ddf, how = 'inner', left_index = True, right_index = True, suffixes = ['_ND','_T2D'])

fontsize = 15
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

fig, ax = plt.subplots(dpi=300)

#ax.errorbar(var_nndf.loc['hsa-miR-127-5p','BSvar'], var_t2ddf.loc['hsa-miR-127-5p','BSvar'], xerr = var_nndf.loc['hsa-miR-127-5p','BSvar_err'], yerr = var_t2ddf.loc['hsa-miR-127-5p','BSvar_err'], fmt = 'o', c = 'orange')
xra = np.linspace(0,1,100)
ax.plot(xra, xra, ls = '--', c = 'k', lw = 1)
ax.set_xlabel("Var(ND miRNA expression)"); ax.set_ylabel("Var(T2D miRNA expression)")
# ax.plot(xra, (xra+0.05)*1.1, ls = '--', c = 'k', lw = 1)
# ax.plot(xra, (xra-0.05)*0.9, ls = '--', c = 'k', lw = 1)
ax.fill_between(xra, (xra-0.05)*0.9,  (xra+0.05)*1.1, color = 'gray', alpha = 0.25, zorder = 0)


top_t2d_variable = var_bs[(var_bs['BSvar_T2D']-var_bs['BSvar_err_T2D'])>1.1*(var_bs['BSvar_ND']+0.05)]
top_t2d_variable = top_t2d_variable[nddf.mean(axis=1)>0] # make sure miRNA is expressed in ND
top_nd_variable = var_bs[(var_bs['BSvar_T2D']+var_bs['BSvar_err_T2D'])<0.9*(var_bs['BSvar_ND']-0.05)]

non_variable = var_bs.loc[~var_bs.index.isin(list(top_t2d_variable.index) + list(top_nd_variable.index)),:]

# plot T2D variable
color='#7D7F99'
ax.scatter(top_t2d_variable['BSvar_ND'], top_t2d_variable['BSvar_T2D'],c=color, edgecolor='black', 
           zorder = 3, alpha = 0.7, s = 20, label='T2D variable')
ax.errorbar(top_t2d_variable['BSvar_ND'], top_t2d_variable['BSvar_T2D'], 
            xerr = top_t2d_variable['BSvar_err_ND'], yerr = top_t2d_variable['BSvar_err_T2D'], 
            lw = 0.5, c = 'silver', zorder = 1, fmt = 'o', ms = 0)

# plot ND variable
color='#F3DD89'
ax.scatter(top_nd_variable['BSvar_ND'], top_nd_variable['BSvar_T2D'],c=color, edgecolor='black', 
           zorder = 3, alpha = 0.7, s = 20, label = 'ND variable')
ax.errorbar(top_nd_variable['BSvar_ND'], top_nd_variable['BSvar_T2D'], 
            xerr = top_nd_variable['BSvar_err_ND'], yerr = top_nd_variable['BSvar_err_T2D'], 
            lw = 0.5, c = 'silver', zorder = 1, fmt = 'o', ms = 0)

# plot non variable
color='#3d3dcc'
ax.scatter(non_variable['BSvar_ND'], non_variable['BSvar_T2D'],c=color, edgecolor='black', 
           zorder = 2, alpha = 0.7, s = 20, label = 'stable/mixed')
ax.errorbar(non_variable['BSvar_ND'], non_variable['BSvar_T2D'], 
            xerr = non_variable['BSvar_err_ND'], yerr = non_variable['BSvar_err_T2D'], 
            lw = 0.5, c = 'silver', zorder = 1, fmt = 'o', ms = 0)


plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

fig.savefig(filenameout, bbox_inches='tight')