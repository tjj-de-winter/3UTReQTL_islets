### FASTQ file processing

Retrieve Smart-Seq FASTQ files, the script assumes that each file contains the transcriptome of a single cell.

Add cell IDs to the read identifier line in the FASTQ files.
run:
```
python add_cell_ID.py --fq fastqs.txt --cellid cell_ID.txt --outdir /path/to/out_folder --outprefix <output prefix>
```
**_EXAMPLE: fastqs.txt:_**
```
sampleA_1_R2.fastq.gz
sampleA_2_R2.fastq.gz
sampleA_3_R2.fastq.gz
...
```
**_EXAMPLE: cell_ID.txt:_**
```
sampleA_1_R2.fastq.gz,cell_1
sampleA_2_R2.fastq.gz,cell_2
sampleA_3_R2.fastq.gz,cell_3
...
```

merge the FASTQ files containing a unique cell identifier into a single FASTQ file.
```
./merge_fastq_files.sh fastq_files.txt <output fastq file> -f
```
**_EXAMPLE: fastq_files.txt:_**
```
sampleA_1_R2_cbc.fastq.gz
sampleA_2_R2_cbc.fastq.gz
sampleA_3_R2_cbc.fastq.gz
...
```

Trim the fastq files
```
trim_galore <input fastq file>
```

align the FASTQ files to the reference genome using STAR
```
./genome_alignment.sh <input fastq file> <output prefix> <sampleID> /path/to/indexed_genome
```

### Transcriptome analysis
Split the sample BAM file into smaller cell BAM files (each file is one cell)
```
./split_bam.sh <sample BAM file> /path/to/out_folder
```

Count the genes using featurecounts to produce a tab-separated count table
```
./count_genes.sh /path/to/single_BAM_files <GTF file> <feature to count> <output path> <output prefix>
```

Make quality control plots using the raw gene count tables
```
count_matrix_QC.py <count TSV file> <sampleID> /path/to/out_folder
```

Perform filtering and processing of the count tables of the different samples and define cell types based on gene expression markers.
***Make sure that the adata.obs object contains a cell-type column which contains the cell type name, a sampleID column which contains the sampleID and/or a group/condition column 
that contains the group name.***

Save the final anndata (adata) object in Python to a .h5ad file. 
```
Notebook: count_processing_celltype_calling
```

### 3'UTR variant genotyping
Perform genotyping according to the GATK pipeline to produce GVCFs
```
./find_variants.sh <BAM file> <Genome FASTA file> <dbSNP VCF> <clinvar VCF> /path/to/out_folder
```
Genotype the sample-specific GVCF files and generate a VCF
```
./genotype_variants.sh /path/to/GVCF <cohort name> <Genome FASTA file> <dbSNP VCF> <output prefix>
```
**_EXAMPLE:_**
```
├── /path/to/GVCF_folder/
│   ├── sampleA.g.vcf.gz
│   ├── sampleB.g.vcf.gz
│   ├── sampleC.g.vcf.gz
...
```

Select called variants from a VCF that map to the 3'UTR and annotate the variants with SNP rsID and gene name
```
./VCF_filtering_annotation.sh <VCF file> <BED 3'UTR annotation file> <dbSNP VCF> <output refix> /path/to/out_folder
```

### Compute eQTLs

Generate a gene filter matrix based on expression thresholds per cell type that is used to filter out eQTLs in lowly expressed genes.
```
python generate_gene_filter.py --h5ad <H5AD file> --celltype_header <celltype header>  --group_header <group header> --cell_and_group --outprefix <output prefix> --outdir /path/to/out_folder
```

Generate a matrix containing genotypes (encoded as integers) of gene-annotated variants per sample from the gene and SNPID annotated 3'UTR VCF.
Output is a comma-separated file.
```
python variant_matrix.py --VCF <VCF file> --outprefix <output prefix> --outdir /path/to/out_folder
```

Generate a matrix with the raw gene expression for each sample present and normalize the raw counts using DESeq
Output is a comma-separated file.
```
python DEseq_normalization.py --count_folder /path/to/TSV_folder/ --h5ad <H5AD file> <output prefix> --outdir /path/to/out_folder
```
**_EXAMPLE:_**
```
├── /path/to/TSV_folder/
│   ├── sampleA_count-matrix.tsv
│   ├── sampleB_count-matrix.tsv
│   ├── sampleC_count-matrix.tsv
...
```

Calculate eQTLs and output the results in an unfiltered comma-separated file
```
python eQTL_computation.py --geno <variant CSV> --genes <gene expression CSV> --h5ad <H5AD file> --cell_type <cell type>  --group <group> --sampleID_header <sampleID header> --celltype_header <celltype header>  --group_header <group header> --cell_and_group --outprefix <output prefix> --outdir /path/to/out_folder   
```

Filter eQTLs and adjust the p-values using the Benjamini-Hochberg FDR p-values, output as a new comma-separated file
```
eQTL-adjusted-pvalue.py --eqtl <eQTL CSV> --filter_file <gene filter CSV> --filter
```

visualize the eQTLs in boxplots
```
Notebook: eQTL_visualization
```

**_Optional:_**
Calculate eQTLs by downsampling the cell type of interest to a certain cell number, output is stored as an unfiltered comma-separated file for each iteration.
Note: The new number of cells has to be lower than the total cells of that cell type in the dataset.
To perform the downsampling 100 times run the following code:
```
for i in {1..100}; do
  python eQTL_computation_downsampled.py --cell_number <number of cells> --iteration <i> --geno <variant CSV> --genes <gene expression CSV> --h5ad <H5AD file> --cell_type <cell type>  --group <group> --sampleID_header <sampleID header> --celltype_header <celltype header>  --group_header <group header> --cell_and_group <output prefix> --outdir /path/to/out_folder
done
```

Average downsampled eQTL p-values, filter and adjust the p-values using the Benjamini-Hochberg FDR p-values, output as a new comma-separated file
```
python eQTL_downsampled_average.py --input_folder /path/to/CSV_folder --filter
```
**_EXAMPLE:_**
```
├── /path/to/CSV_folder/
│   ├── celltypeA_control_1_eQTL.csv
│   ├── celltypeA_control_2_eQTL.csv
│   ├── celltypeA_control_3_eQTL.csv
...
```

### Compute eQTL-miRNA binding site creation or destruction
Given a list of miRNAs, the following scripts compute if miRNA binding sites are created or destroyed by eQTL variant sequences.

Generates two FASTA files (1: reference alleles, 2: alternative alleles) with sequences flanking eQTLs for each eQTL in the CSV file
```
./eQTL_to_fasta.sh <filtered eQTL CSV> <variant VCF> <GTF BED> <Genome FASTA file> /path/to/out_folder <output prefix>
```

Next, compute miRNA binding using miRanda software. 
Mature miRNA sequences can be downloaded from miRBase. The code will use all miRNAs in the FASTA file, so a smaller FASTA can be made if only a subset of miRNA needs to be assessed.
unfiltered miRanda results are generated as a comma-separated file.
```
eQTL_miRanda.sh <reference alleles FASTA> <alternative alleles FASTA> <mature miRNA FASTA> /path/to/out_folder <output prefix>
```

Filter the output to only select binding sites that are altered by eQTL variants.
Note: VCF files need to split into a file specific for each chromosome to speed up the computation.
```
python eQTL_miRanda_filter.py --miranda_file <unfiltered miRanda CSV> --fasta_file <merged alternative and reference alleles FASTA>
```

Perform a second filtering step by selecting miRNA-eQTL pairs where miRNA binding leads to decreased gene expression (expected behavior of miRNA targeting).
Output is a single filtered comma-separated file for each cell type.
```
eQTL_miRanda_slope.py --miranda_file <miRanda CSV> --eQTL_file <filtered eQTL CSV> --h5ad <H5AD file> --celltype_header <celltype header> --group_header <group header> --sampleID_header <sampleID header> --cell_and_group --VCF_folder /path/to/VCF_folder
```

To split VCFs per chromosome run:
```
vcf=<VCF file>
for chr in $(bcftools query -f '%CHROM\n' ${vcf} | sort -u); do
    bcftools view -r ${chr} ${vcf} -o ${vcf%.vcf})_chr${chr}.vcf
done
```
**_EXAMPLE:_**
```
├── /path/to/VCF_folder/
│   ├── cohortA_3UTR_gene_annotated_chrX.vcf
│   ├── cohortA_3UTR_gene_annotated_chrY.vcf
│   ├── cohortA_3UTR_gene_annotated_chr1.vcf
│   ├── cohortA_3UTR_gene_annotated_chr2.vcf
...
```
