#!/bin/bash

for csv in $(ls eQTL_variants_extracted2_96samples/*.csv)
do
	echo $csv
	python genotype.py $csv VCF
done

bcftools merge VCF/*.gz -o eQTL_variants_extracted2_96samples.vcf.gz