snp=/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/VCF/combined/smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.bed
utr=/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/6-support_files/BED/Homo_sapiens.GRCh38.102.minusINS-IGF2_uniq3utr_fromstop.three_prime_utr.bed

#snp is made by bin/VCF_to_BED.sh

bedtools intersect -a $utr -b $snp -u > 3UTR.bed

bedtools intersect -a 3UTR.bed -b $snp > SNP_3UTR.bed

for motif_bed in $(ls ../AREsite2/*.chrmod.bed.gz)
do
        name=$(basename $motif_bed)
        outname=${name%.bed*}.3UTR.bed

        bedtools intersect -a $motif_bed -b 3UTR.bed -u > $outname
done

motif_bed=../polyasite/atlas.clusters.3.0.GRCh38.GENCODE_42.chrmod.PAS.window.bed.gz
name=$(basename motif_bed)
outname=${name%.bed*}.3UTR.bed
bedtools intersect -a $motif_bed -b 3UTR.bed -u > $outname
