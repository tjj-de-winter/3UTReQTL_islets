# 3UTReQTL_islets

Code and Jupyter notebooks related to ###reference to paper###.

### This repository provides the codes to:
* Process SMART-seq single cell RNA-seq (scRNA-seq) FASTQ files into BAM alignment files
* Genotype 3'UTR variants from BAM
* Transcriptome analysis from BAM
* Cell type specific-eQTL computation
* miRNA binding prediction using eQTL sites

### Usage Instructions
See usage.md for detailed examples on how to use the scripts in this repository.

### Requiered public data files:
* GRCh38 human reference genome: https://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/
* GRCh38 human GTF file: https://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens/
* Clinvar variant VCF file: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
* dbSNP variant VCF file: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/
* Mature miRNA sequences: https://www.mirbase.org/download/mature.fa

### Requiered software packages:
* Python v3.8.5
* R v4.2.1
* STAR v2.7.7a (https://github.com/alexdobin/STAR)
* samtools v1.16.1 (https://github.com/samtools/samtools)
* picard v2.25.0 (https://github.com/broadinstitute/picard)
* GATK v4.1.9.0 (https://github.com/broadinstitute/gatk)
* featureCounts v2.0.1 (https://github.com/ShiLab-Bioinformatics/subread)
* bedtools v2.30.0 (https://github.com/arq5x/bedtools2)
* bcftools v1.17 (https://github.com/samtools/bcftools)
* vcftools v0.1.13 (https://github.com/vcftools/vcftools)
* miRanda v3.3a (https://github.com/hacktrackgnulinux/miranda)
* pysam v0.19.1 (https://github.com/pysam-developers/pysam)
* PyVCF v0.6.8 (https://github.com/jamescasbon/PyVCF/)
* pandas v1.5.1 (https://github.com/pandas-dev/pandas)
* numpy v1.23.3 (https://github.com/numpy/numpy)
* scanpy v1.9.3 (https://github.com/scverse/scanpy)
* decoupler v1.6.0 (https://github.com/saezlab/decoupler-py)
* matplotlib v3.7.1 (https://github.com/matplotlib/matplotlib)
* scipy v1.10.1 (https://github.com/scipy/scipy)
* PyComplexHeatmap v1.7.7 (https://github.com/DingWB/PyComplexHeatmap)
* scikit-learn v1.3.2 (https://github.com/scikit-learn/scikit-learn)
* UpSetPlot v0.9.0 (https://github.com/jnothman/UpSetPlot)
* statsmodels v0.14.1 (https://github.com/statsmodels/statsmodels)
* pandarallel v1.6.1 (https://github.com/nalepae/pandarallel)
* func_timeout v4.3.5 (https://github.com/kata198/func_timeout)
* rpy2 v3.5.4 (https://github.com/rpy2/rpy2)
* SCeQTL v0.2.0 (https://github.com/XuegongLab/SCeQTL)
* DESeq2 v1.38.3 (https://github.com/thelovelab/DESeq2)
