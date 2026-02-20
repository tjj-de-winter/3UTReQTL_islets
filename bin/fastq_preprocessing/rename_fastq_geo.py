#!/usr/bin/python3

### Description ###

#update fastq file sample name using the updated sraRunTable sample title column
#please run update_SraRunTable.py

### Import packages ###
import argparse as argp
import pandas as pd
import os

### Functions ###

### Input variables ###
parser = argp.ArgumentParser(description = 'rename fastq file using sraRunTable information')
parser.add_argument('--fq', help = 'fastq file name in SRRxxx_x.fastq.gz format')
parser.add_argument('--sratable', help = ' SraRunTable_update_GSExxxx.csv')
args = parser.parse_args()

fq = args.fq
sratable = args.sratable

### Code ###

dirname = os.path.dirname(fq)
fq = os.path.basename(fq)

SRR = fq.split('_')[0]

df = pd.read_csv(sratable)
sample_title = list(df[df['Run'] == SRR]['sample_title'])

if len(sample_title) == 1:
    fq_out = '-'.join([sample_title[0], fq])

    # rename the files
    fq = '/'.join([dirname, fq])
    fq_out = '/'.join([dirname, fq_out])
    print(fq, fq_out)
    os.rename(fq, fq_out)
