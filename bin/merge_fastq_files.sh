#!/bin/bash

### Description ###

# merge a list of files (e.g. fastq files) into a single output file, automatic gzipping if the original file input was gzipped. 

### Input variables ###
if [ $# -ne 3 ]
then
    echo "Usage: $0 <file list txt> <output name> -f"
    echo "  1) Headerless text file with a list of files to merge"
    echo "  2) output file name"
    echo "  3) overwrite existing merged file? use '-f' or leave empty if not"
    exit 1
fi

file_list=$1 
output=$2 
overwrite=$3

### Code ###

# check if output file already exists, stop if not forced to overwrite
if [ -f $output ]
then 
	if [ -z $overwrite ]
	then
		echo "file '${output}' already exists"
		exit
	else
		rm $output
	fi
fi

# if output file exist but forced to overwrite check if unzipped file exist and remove

if [ -f ${output%.gz} ]
then
        if [ -z $overwrite ]
        then
                echo "file '${output%.gz}' already exists"
                exit
        else
                rm ${output%.gz}
        fi
fi

# merge files, check if unput exists and whether file is zipped (cat vs zcat) 
for file in $(awk '{print $NF}' $file_list)
do 
	if [ -f $file  ]
	then
		echo ""
	else
		echo "input file '${file}' does not exist"
		rm ${output%.gz}
		exit
	fi

	gz=${file: -3}
	if [ ${gz} != ".gz" ]
	then
		echo "start merging '${file}'"
		cat $file >> ${output%.gz}
	else
		echo "start merging '${file}'"
		zcat $file >> ${output%.gz}
	fi
done

# if output file is written as a zipped file, zip the newly made merged file

gzip=${output: -3}
if [ $gzip = ".gz" ]
then
	echo "start gzip of '${output%.gz}'"
	gzip ${output%.gz}
fi

exit
