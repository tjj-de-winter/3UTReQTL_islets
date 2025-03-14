#!/bin/python

###########################
##  eQTL SUPP Figure 4C  ##
###########################

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

filenameout = 'supp_figure4c.png'

# Load the data
genotype_csv = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/functional validation/genotyping/islets_genotype.csv'
df = pd.read_csv(genotype_csv, sep=",", index_col=0)
df.sort_values(by='VIC-T', inplace=True)

# Define markers and assign them to the dataframe
markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
df['marker'] = markers[:df.shape[0]]

# Set up the plot
fig, ax = plt.subplots(dpi=300, figsize=(5,5), constrained_layout=True)

colors = ['r', 'g', 'b', 'y']

# Store the legend handles
color_legend_elements = []
marker_legend_elements = []

# Plot the data
for ii, (geno, df_sub) in enumerate(df.groupby('Genotype')):
    color = colors[ii % len(colors)]  # Cycle through colors
    color_legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=geno))  # Add to genotype legend
    
    for i, idx in enumerate(df_sub.index):
        donor = df_sub.iloc[i, 0]
        marker = df_sub.iloc[i, 4]
        xi = df_sub.iloc[i, 1]
        yi = df_sub.iloc[i, 2]
        
        # Plot the points with different markers and colors, adding edgecolor
        ax.scatter(xi, yi, marker=marker, color=color, edgecolor='black', zorder=10)
        
        # Add unique donors to the marker legend
        if donor not in [h.get_label() for h in marker_legend_elements]:
            marker_legend_elements.append(Line2D([0], [0], marker=marker, color='black', markersize=10, linestyle='None', label=donor))

# Set labels, title, and axis limits
ax.set_ylabel('Normalized FAM fluorescence (rs701848-C)')
ax.set_xlabel('Normalized VIC fluorescence (rs701848-T)')
ax.set_title('PTEN rs701848')
ax.set_ylim(-150, 1150)
ax.set_xlim(-120, 2400)

# Add the genotype (color) legend
legend1 = ax.legend(handles=color_legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), title="Genotype")

# Add the donor (marker) legend
legend2 = ax.legend(handles=marker_legend_elements, loc='lower left', bbox_to_anchor=(1.05, 0), title="Donor")

# Make sure both legends are added properly without overwriting
ax.add_artist(legend1)

fig.savefig(filenameout, bbox_inches='tight')
