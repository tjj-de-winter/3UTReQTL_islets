#!/bin/bash

fq=$1

avg_readlength=$(zcat $fq | head -n 1000000 | awk 'NR % 4 == 2 { total += length($0); count++ } END { if (count > 0) print total / count; else print 0 }')

echo "${fq},${avg_readlength}" >> read_length.csv

echo $fq average read length $avg_readlength