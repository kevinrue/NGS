#!/bin/bash

# general form
# fastq-dump [options] <accession>

# Useful options
#  -O	|	--outdir <path>	Output directory, default is current working directory ('.').

# SRR to download
srr_list=(SRR2038194 SRR2038195 SRR2038196 SRR203819)

for srr_id in ${srr_list[@]}
do
        echo "srr_id: ${srr_id}"
        cmd="mkdir -v ${srr_id}"
        echo "cmd: ${cmd}"
        eval $cmd                                                                                                                                                                           
        cmd="fastq-dump -O ${srr_id} --bzip2 ${srr_id}"
        echo "cmd: ${cmd}"
        eval $cmd
done

echo "Complete"
