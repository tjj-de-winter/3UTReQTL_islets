# determine MAF from 1000genomes VCF file, convert to csv file

import vcf
import argparse as argp

### Input variables ###
parser = argp.ArgumentParser(description = 'Retrieve minor allele frequency from VCF file')
parser.add_argument('--vcf', help = '1000genomes VCF file')

args = parser.parse_args()

vcf_file = args.vcf

vcf_reader = vcf.Reader(filename=vcf_file)

def get_MAF(alt_freq):
	'''record.aaf calculates the frequency of the alternative allele
	   to get the MAF check if the alt or ref is the minor allele'''
	MAF = min([alt_freq, 1-alt_freq])
	return MAF

with open(vcf_file.replace('.vcf.gz','.MAF.csv'),'w') as csv_out:
	for record in vcf_reader:
		chrom = record.CHROM
		pos = record.POS
		ref = record.REF
		alt = record.ALT
		MAF =  get_MAF(record.aaf[0])

		csv_out.write(f'{chrom}\t{pos}\t{ref}\t{alt}\t{MAF}\n')

print(vcf_file.replace('.vcf.gz','.MAF.csv'), ' made')

csv_out.close()