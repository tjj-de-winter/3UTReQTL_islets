#!/bin/python

## Description ###

# Downloading VCF files from 1000genomes 

import os

chromosomes = [i+1 for i in range(22)]
chromosomes.append('X')
chromosomes.append('Y')

for chromosome in chromosomes:
	file = f"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr{chromosome}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
	filename = os.path.basename(file)
	print(filename)
	if not os.path.isfile(filename):
		print(f'start download VCF chromosome {chromosome}')
		os.system(f'wget {file}')
	else:
		print(f'VCF chromosome {chromosome} already exists')
	if not os.path.isfile(filename+'.csi'):
		os.system(f'bcftools index {filename}')