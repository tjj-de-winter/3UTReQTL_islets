# make variant matrix, encoding the SNP as the gene dosage of the alternative allele

import argparse as argp
import vcf
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
import glob
import os
warnings.simplefilter(action='ignore', category=Warning)

### Input variables ###

parser = argp.ArgumentParser(description='Generates a variant matrix, encoding the SNP as the gene dosage of the alternative allele')
parser.add_argument('--vcf', help = 'VCF file path (supports wildcards).', type = str, default = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/VCF/*.smartseqPancreasv2_3UTR_gene_annotated_MAF0.07.vcf')
parser.add_argument('--h5ad', help = 'Input .h5ad file path.', type = str, default = '/Users/twanw/Documents/LUMC/projects/3UTRpancreas/Data_3UTR_project/revisions/notebook/eQTL_revision.h5ad')
parser.add_argument('--celltype', help = 'Cell type to subset (e.g., beta-cells).', type = str, default = 'beta-cells')
parser.add_argument('--group', help = 'Group label to subset (e.g., ND).', type = str, default = 'ND')
parser.add_argument('--outdir', default='.',help='output path')
parser.add_argument('--prefix', default='.',help='filename prefix')

args = parser.parse_args()

VCFfile = args.vcf
h5ad_file = args.h5ad
cell_type = args.celltype
group = args.group
outdir = args.outdir
prefix = args.prefix

if not os.path.exists(outdir):
        os.makedirs(outdir)

adata = sc.read_h5ad(h5ad_file)

adata_sub = adata[(adata.obs['Cell_type'] == cell_type)& (adata.obs['Group'] == group)]

donors = list(set(adata_sub.obs['DonorID']))
minimum_alt_carriers = 5


def check_gene_dosage_criteria(genotype_dict, donors):
    REF = set()
    ALT = set()
    REF_allele = ''
    ALT_allele = ''
    alternative_alleles = 0
    reference_alleles = 0
    for donor in genotype_dict:
        if donor in donors:
            genotypes = genotype_dict[donor]
            geno_bases = str(genotypes[0]).replace('|','/').split('/')
            dosage = genotypes[1]
            if dosage == 0 or dosage == 1:
                REF.add(geno_bases[0])
                REF_allele = geno_bases[0]
                reference_alleles += 1
            elif dosage == 1:
                ALT.add(geno_bases[1])
                alternative_alleles += 1
                ALT_allele = geno_bases[1]
            elif dosage == 2:
                ALT.add(geno_bases[0])
                alternative_alleles += 1
                ALT_allele = geno_bases[0]
    if len(ALT) == 1 and "*" not in ALT:
        return reference_alleles, alternative_alleles, REF_allele, ALT_allele
    else:
        return 0, 0, 0, 0 # more than 1 alternative allele or not defined (e.g. *)

def allele_dosage(genotype_dict, donors):
    allele_dosage_dict = {}
    for donor in genotype_dict:
        if donor in donors:
            genotypes = genotype_dict[donor]
            dosage = genotypes[1]
            if dosage != None:
                allele_dosage_dict[donor] = dosage
    return allele_dosage_dict

cells = list(adata_sub.obs.index)
donors_cells = list(adata_sub.obs.DonorID)

dummy_df = pd.DataFrame(np.array(donors_cells), index=cells, columns=['Donor'])

def dict_to_df(dict, variant, dummy_df=dummy_df):
	alleles = []
	for donor in dummy_df['Donor']:
		if donor in dict.keys():
			dosage = dict[donor]
			alleles.append(int(dosage))
		else:
			alleles.append(np.nan)

	dummy_df[variant] = alleles
	dummy_df = dummy_df.loc[:,[variant]]
	return dummy_df

variant_matrix = pd.DataFrame(index=donors_cells)

VCFfiles = glob.glob(VCFfile)

first_variant = True
for VCFi in VCFfiles:
	VCF = vcf.Reader(fsock=None, filename=VCFi, prepend_chr=False, strict_whitespace=False, encoding='ascii')

	basename = os.path.basename(VCFi)
	chromosome_number = basename.split('.')[0]
	print(f'start reading chromosome {chromosome_number}')

	for ii, record in enumerate(VCF):
	    if 'Gene' in record.INFO:
		    for gene in record.INFO['Gene']:

		        genotype_dict = {}
		        for sample in record:
		            genotype_dict[sample.sample] = [sample.gt_bases, sample.gt_type]

		        reference_alleles, alternative_alleles, REF, ALT = check_gene_dosage_criteria(genotype_dict, donors)
		        if  reference_alleles >= minimum_alt_carriers and alternative_alleles >= minimum_alt_carriers:
		            variant = '_'.join([f'{record.CHROM}:{record.POS}',REF,ALT,gene])
		            allele_dosage_dict = allele_dosage(genotype_dict, donors)
		            if first_variant:
		            	variant_matrix = dict_to_df(allele_dosage_dict, variant)
		            	first_variant = False
		            else:
		            	variant_matrix = pd.concat([variant_matrix, dict_to_df(allele_dosage_dict, variant)], axis=1)

outfile = f'{cell_type}_{group}_variant_matrix.csv'

if prefix != None:
    outfile = prefix + "_" + outfile

if outdir != None:
    outfile = os.path.join(outdir, outfile)

variant_matrix.to_csv(outfile)
