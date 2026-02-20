#!/usr/bin/python3

### Description ###

#update sraRunTable with sample title, using the corresponding soft.gz meta data

### Import packages ###
import pandas as pd
import gzip

### Functions ###

def combine_metadata(soft_file, sraRunTable):
    soft = gzip.open(soft_file)
    sample_name_dict = {}
    for line in soft.readlines():
        line = line.decode("utf-8").strip('\n')
        if '^SAMPLE = ' in line:
            sample = line.split(' = ')[-1]
        if '!Sample_title = ' in line:
            title = line.split(' = ')[-1]
            sample_name_dict[sample] = title


    df = pd.read_csv(sraRunTable)
    df['sample_title'] = df['Sample Name'].apply(lambda GSM: sample_name_dict[GSM])
    return df

### Code ###

# GSE124742
ID = 'GSE124742'
soft_file = f'{ID}_family.soft.gz'
sraRunTable = f'SraRunTable_{ID}.csv'
df = combine_metadata(soft_file, sraRunTable)
df['sample_title'] = [s.split(',')[0] for s in df['sample_title']]
filenameout = f'SraRunTable_update_{ID}.csv'
df.to_csv(filenameout, index=False)

# GSE164875
ID = 'GSE164875'
soft_file = f'{ID}_family.soft.gz'
sraRunTable = f'SraRunTable_{ID}.csv'
df = combine_metadata(soft_file, sraRunTable)
df['sample_title'] = [s.split('[')[-1].split(']')[0] for s in df['sample_title']]
filenameout = f'SraRunTable_update_{ID}.csv'
df.to_csv(filenameout, index=False)

# GSE270484
ID = 'GSE270484'
soft_file = f'{ID}_family.soft.gz'
sraRunTable = f'SraRunTable_{ID}.csv'
df = combine_metadata(soft_file, sraRunTable)
filenameout = f'SraRunTable_update_{ID}.csv'
df.to_csv(filenameout, index=False)