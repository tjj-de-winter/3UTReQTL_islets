#!/bin/python

###########################
##    eQTL Figure 3B     ##
###########################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

filenameout = 'figure3b.png'

def get_pca(X, nc = 55):
    pca = PCA(n_components = nc).fit(X.T)
    rpca = []
    
    explained_variance_percentages = pca.explained_variance_ratio_ * 100  # Convert to percentages
    for i in range(20):
        RX = X.copy(); RX = np.array(RX)
        for x in RX:
            np.random.shuffle(x)
        RX = pd.DataFrame(RX)
        rpca.append(PCA(n_components = nc).fit(RX.T))
    rand_var_expl = pd.DataFrame([x.explained_variance_ratio_ for x in rpca])
    return pca, rand_var_expl, explained_variance_percentages

file = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/8-miRNA_seq/count_matrix/deseq2_norm_data_AA.tsv'
df_ct = pd.read_csv(file, index_col=0, sep="\t")
df_ct.columns = [col.rsplit('.')[0] for col in df_ct]

df_ct_norm = np.log(df_ct+1)

norm_pca, norm_rnd_pcavar, explained_perc = get_pca(df_ct_norm, nc = 5)

df = pd.DataFrame(norm_pca.fit_transform(df_ct_norm.T), index = df_ct_norm.columns).T # pca components per sample

fontsize = 20
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

fig, ax = plt.subplots(dpi=300, figsize=(5,5))
for label,c in zip(['ND','T2D'],['#F3DD89','#7D7F99']):
    cols = [c for c in df.columns if label in c]
    ax.scatter(df.loc[0,cols], df.loc[1,cols],s = 50, c = c, edgecolor='black',label = label)

pc1 =round(explained_perc[0],1); pc2=round(explained_perc[1],1)  
ax.set_xlabel(f'PC 1 ({pc1}%)'); ax.set_ylabel(f'PC 2 ({pc2}%)')
for c in df.columns:
    name = c.replace('ND','').replace('T2D','')
    ax.text(df.loc[0,c], df.loc[1,c], '  '+name, va = 'center', ha = 'left', fontsize = 15)
ax.set_xlim(-5, 13)
ax.legend()

fig.savefig(filenameout, bbox_inches='tight')