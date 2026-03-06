#!/bin/bash

sra_txt=$1 # SRR_Acc_list.txt file
outdir=$2
sra_dir=$3

for srr in $(awk '{print $0}' $sra_txt);
do
	fastq-dump --outdir $outdir --split-files --gzip ${sra_dir}/${srr}.sra;
done





