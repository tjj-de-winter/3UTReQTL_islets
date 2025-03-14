#!/bin/python

###########################
##    eQTL Figure 1C     ##
###########################
import matplotlib.pyplot as plt
import sys,os

#
filenameout = 'figure1c.png'

# set fontsize
fontsize = 30
plt.rcParams['axes.titlesize'] = fontsize
plt.rcParams['axes.labelsize'] = fontsize
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['figure.titlesize'] = fontsize
plt.rcParams['font.size'] = fontsize

# count SNPs and indels for variants and eQTLs
vcf_vars = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/2-variants/1-VCF/smartseqPancreas_3UTR.vcf'
vcf_eqtls = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/2-variants/1-VCF/smartseqPancreas_3UTR_gene_annotated.eQTLs.vcf'

data_file = 'data.txt' 
snps_vars = os.system(f'bcftools view -v snps {vcf_vars} | grep -v "^#" | wc -l > {data_file}')
all_vars = os.system(f'bcftools view {vcf_vars} | grep -v "^#" | wc -l >> {data_file}')
snps_eqtls = os.system(f'bcftools view -v snps {vcf_eqtls} | grep -v "^#" | wc -l >> {data_file}')
all_eqtls = os.system(f'bcftools view {vcf_eqtls} | grep -v "^#" | wc -l >> {data_file}')

f = open(data_file, 'r')
data = []
for line in f.readlines():
	data.append(int(line.strip('\n')))

# idx 0 = SNP variants
# idx 1 = All variants
# idx 3 = SNP eQTLs
# idx 4 = All eQTLs

variants = [data[1]-data[0], data[0]]
eqtls = [data[3]-data[2], data[2]]

print('variants', variants)
print('eQTLs', eqtls)

os.system(f'rm {data_file}')

fig, ax = plt.subplots(dpi=300, figsize=(12,4))

ax.barh([1,0], [variants[0], eqtls[0]], color=["#333680", "#386633"], edgecolor='black')
ax.barh([3,3], [variants[0], eqtls[0]], color='gray', edgecolor='black', label='Indel variants')

ax.barh([1,0], [variants[1], eqtls[1]],left=[variants[0], eqtls[0]], hatch='//', color=["#6E6EB1", "#619A41"], edgecolor='black')
ax.barh([3,3], [variants[1], eqtls[1]],left=[variants[0], eqtls[0]], hatch='//', color='white', edgecolor='black', label='SNP variants')


ax.grid(False)
ax.set_yticks([1,0], ["3'UTR variants", "3'UTR eQTLs"])
ax.set_xlabel('Count')
ax.set_ylim(-0.5,1.5)

# ax.set_xticks([100000,200000,300000,400000,500000], [100000,200000,300000,400000,500000]);

ax.legend(fontsize=fontsize)

ax.tick_params(axis='x', which='major', length=10, width=1)  # Length and width of major ticks

fig.savefig(filenameout, bbox_inches='tight')