#!/bin/bash

## basename_1.fastq.gz -> basename.fastq.1.gz
## basename_2.fastq.gz -> basename.fastq.2.gz

for fullfilename in $(find ${raw_fastq_dir} -name '*.fastq.gz')
do
	echo "Processing ${fullfilename}"
	basefilename=$(basename ${fullfilename})
	echo "basefilename: ${basefilename}"
	newbasename=$(echo ${basefilename} | sed -re 's|_(.)\.fastq\.gz|.fastq.\1.gz|')
	echo "newbasename: ${newbasename}" 
done
