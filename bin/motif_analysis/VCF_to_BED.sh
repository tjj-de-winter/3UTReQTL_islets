#!/bin/bash

VCF=$1
BED=${VCF%.vcf*}.bed

bcftools query -f '%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT\t%INFO\n' $VCF \
| awk '{OFS="\t"; print $1, $2-1, $3, $4, $5, $6, $7}' > $BED