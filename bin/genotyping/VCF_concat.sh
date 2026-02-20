#!/bin/bash

#SBATCH --job-name=concat
#SBATCH --output=concat.out
#SBATCH --error=concat.err
#SBATCH --time=5-24:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1

VCF_path=/exports/ana-scarlab/tdwinter/3UTRpancreas/eQTL_revision/GVCFs_R1_R2/out

bcftools concat \
        ${VCF_path}/1.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/2.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/3.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/4.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/5.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/6.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/7.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/8.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/9.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/10.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/11.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/12.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/13.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/14.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/15.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/16.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/17.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/18.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/19.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/20.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/21.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/22.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        ${VCF_path}/X.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz \
        -Oz -o ${VCF_path}/smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf.gz