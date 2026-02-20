#!/usr/bin/python3

# description based on the average FASTQ file read length retrieve which STAR index is most suitable
# run first: ./read_lengt.sh

import pandas as pd
import os

df = pd.read_csv('read_length.csv', names=['file','read_length'])
df['donor'] = df['file'].apply(lambda x: os.path.basename(x).split('_R1_')[0])

def get_index(rl):
        if rl <55:
                index = 49
        elif rl >= 55 and rl <= 75:
                index = 74
        elif rl > 75:
                index = 99
        return index

df['STAR_index'] = df['read_length'].apply(lambda x: get_index(x))

df = df.loc[:,['donor','read_length','STAR_index']]

df.to_csv('STARindex.csv')
print(df)
