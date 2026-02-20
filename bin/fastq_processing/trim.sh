#!/bin/bash

fastq=$1
outdir=trimmed_fastq

mkdir -p $outdir

sbatch -t 5-05:00:00 --mem=25G -J trim -e ${fastq}_trim.err -o ${fastq}_trim.out --wrap="trim_galore -o ${outdir} $fastq"

echo "trimming $fastq done"
