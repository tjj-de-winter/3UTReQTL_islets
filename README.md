# 3UTReQTL_islets

This repository contains the code for: 

Cell type-specific eQTL detection from single-cell RNA-seq reveals post-transcriptional regulatory mechanisms in human islets

bioRxiv 2025.01.21.633530; doi: https://doi.org/10.1101/2025.01.21.633530

### This repository provides the codes to:
* Retrieving and Processing SMART-seq single cell RNA-seq (scRNA-seq) FASTQ files into BAM alignment files 
- `/bin/data_download`
- `/bin/fastq_preprocessing`
- `/bin/fastq_processing`
- `/bin/mapping`
* Genotype 3'UTR variants from BAM
- `/bin/genotyping`
* Transcriptome analysis from BAM
- `/bin/transcriptomics`
* Cell type specific-eQTL computation
- `/bin/eQTLcomputation`
* eQTL regulatory motif analysis
- `/bin/motif_analysis`
* Correlation with omics data from HumanIslets.com
- `/bin/humanislets_data`
* miRNA binding prediction using eQTL sites
- `/bin/mirna_prediction`

### Public data files used in this repository:
* GRCh38 human reference genome: https://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/
* GRCh38 human GTF file: https://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/
* Clinvar variant VCF file: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
* 1000 genomes VCF files: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/
* dbSNP variant VCF file: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/
* Mature miRNA sequences: https://www.mirbase.org/download/mature.fa
* AREsite2 motif BED files: http://rna.tbi.univie.ac.at/AREsite2/bulk
* HumanIslets.com omics data: https://www.humanislets.com/#/download

### Requiered software packages:
* Python v3.8.5
* R v4.2.1
* STAR v2.7.7a (https://github.com/alexdobin/STAR)
* samtools v1.16.1 (https://github.com/samtools/samtools)
* picard v2.25.0 (https://github.com/broadinstitute/picard)
* OpenJDK v8 (https://github.com/openjdk/jdk8u)
* GATK v4.1.9.0 (https://github.com/broadinstitute/gatk)
* featureCounts v2.0.1 (https://github.com/ShiLab-Bioinformatics/subread)
* bedtools v2.30.0 (https://github.com/arq5x/bedtools2)
* bcftools v1.17 (https://github.com/samtools/bcftools)
* vcftools v0.1.13 (https://github.com/vcftools/vcftools)
* miRanda v3.3a (https://github.com/hacktrackgnulinux/miranda)
* PyVCF v0.6.8 (https://github.com/jamescasbon/PyVCF/)
* pandas v1.5.1 (https://github.com/pandas-dev/pandas)
* numpy v1.23.3 (https://github.com/numpy/numpy)
* scanpy v1.9.3 (https://github.com/scverse/scanpy)
* matplotlib v3.7.1 (https://github.com/matplotlib/matplotlib)
* scipy v1.10.1 (https://github.com/scipy/scipy)
* scikit-learn v1.3.2 (https://github.com/scikit-learn/scikit-learn)
* statsmodels v0.14.1 (https://github.com/statsmodels/statsmodels)
* pandarallel v1.6.1 (https://github.com/nalepae/pandarallel)
* func_timeout v4.3.5 (https://github.com/kata198/func_timeout)
* rpy2 v3.5.4 (https://github.com/rpy2/rpy2)
* biopython v1.79 (https://github.com/biopython/biopython)
* SCeQTL v0.2.0 (https://github.com/XuegongLab/SCeQTL)
* pscl v1.5.9 (https://github.com/atahk/pscl)
* limma v3.54.2 (https://git.bioconductor.org/packages/limma)
* mygene v3.2.2 (https://github.com/biothings/mygene.py)