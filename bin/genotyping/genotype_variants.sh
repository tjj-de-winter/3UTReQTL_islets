#!/bin/bash

#SBATCH --job-name=genotype
#SBATCH --output=genotype.R1.R2.out
#SBATCH --error=genotype.R1.R2.err
#SBATCH --time=3-24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1

### Description ###

# This script takes multiple GVCF files and outputs a filtered VCF file for a specific chromosome including SNPs and Indels.
# The code is based on the best Best Practices Workflows from GATK

### Input variables ###

if [ $# -ne 5 ]
then
    echo "Usage: $0 <chromosome> <cohort> <Genome FASTA file> <dbSNP VCF> <output prefix>"
    echo "1) chromosome to genotype"
    echo "2) cohort name"
    echo "3) reference genome FASTA file, unzipped"
    echo "4) dbSNP feature file (.vcf/.vcf.gz)"
    echo "5) output path"
    exit
fi

chromosome=$1
cohort_name=$2
reference=$3
# for reference use /exports/ana-scarlab/group_references/ensembl/human/102/Homo_sapiens.GRCh38.dna.primary_assembly_ERCC.fa
dbsnp=$4
outpath=$5
outpath=${outpath%/}

### Paths to utilities and scripts ###

GVCFspath=${GVCFspath%/}
outpath=${outpath%/}
p2gatk=$(which gatk)
p2java=$(which java)

### Code ###

# check java version, needs to be JDK version 8 to work with GATK 4.1.9.0

v="$(${p2java} -version 2>&1)"
version=$(echo ${v} | awk -F "." '/version/ {print $2}')
if [ ${version} -ne "8" ]
then
    echo "Java JDK needs to be version 8"
    exit
fi

# genotype GVCF files
${p2gatk} --java-options "-Xmx80g" GenotypeGVCFs \
    --reference $reference \
    --intervals $chromosome \
    --dbsnp $dbsnp \
    --variant ${outpath}/${cohort_name}.g.vcf.gz \
    --output ${outpath}/${chromosome}.${cohort_name}.vcf.gz

exit
