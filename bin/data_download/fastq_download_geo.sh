#!/bin/bash
# download fastq files from patch-seq datasets processed with smartseq2

Accesion=( GSE83139, GSE81547, GSE81608, GSE86469, GSE73727, GSE124742, GSE164875, GSE270484 )
for ID in "${Accesion[@]}"
do
	echo "$ID"
	sbatch --export=All -J ${ID#GSE} -e ${ID}.err -o ${ID}.out -t 7-24:00:00 --mem=30G --wrap="prefetch -X 100G -p --option-file SRR_Acc_List_${ID}.txt"
done