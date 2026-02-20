#!/bin/bash

# generate fastq files from the downloaded .sra files from patch-seq datasets processed with smartseq2

sra_dir=/exports/ana-scarlab/tdwinter/SRA_downloads/sra

Accesion=( GSE83139, GSE81547, GSE81608, GSE86469, GSE73727, GSE124742, GSE164875, GSE270484 )
for ID in "${Accesion[@]}"
do
	echo "$ID"
	outdir=${ID}_fastq
	mkdir -p $outdir

	sbatch --export=All -J ${ID#GSE} -e ${ID}.err -o ${ID}.out -t 7-24:00:00 --mem=30G --wrap="./sra2fastq.sh SRR_Acc_List_${ID}.txt $outdir $sra_dir"
done





