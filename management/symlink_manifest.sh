#!/bin/bash

#############
# Objective #
#############
# Set up symlinks to raw data
# Manifest should have two columns:
# 1: basename for fastq files
# 2: basename for sample

# Notes
# - Currently hard-coded for paired data
# - Currently hard-coded to disambiguate fastq files from samples sequenced on multiple lanes
#   (_i.e._, that need to be later concatenated)


# settings ----

sourceDir="/gfs/archive/sansom/rawdata/mtec/tmp_adam"
manifestFile="manifest"

destinationDir="/gfs/work/kralbrecht/04_lincRNAs/data"


# process inputs ----

baseFastqs=($(cut -f1 ${sourceDir}/${manifestFile} | xargs))
baseSamples=($(cut -f2 ${sourceDir}/${manifestFile} | xargs))

# echo ${baseFastqs[*]} # debug
# echo ${baseSamples[*]} # debug

for rowIndex in $(seq 1 ${#baseFastqs[*]})
do
	let itemIndex=${rowIndex}-1
	echo "itemIndex: ${itemIndex}"

	baseFastq=${baseFastqs[${itemIndex}]}
	echo "baseFastq: ${baseFastq}"

	baseSample=${baseSamples[${itemIndex}]}
	echo "baseSample: ${baseSample}"

    fastq_1="${sourceDir}/${baseFastq}_1.fastq.gz"
    fastq_2="${sourceDir}/${baseFastq}_2.fastq.gz"

    sampleFile_1="${destinationDir}/${baseSample}-${rowIndex}_1.fastq.gz"
    sampleFile_2="${destinationDir}/${baseSample}-${rowIndex}_2.fastq.gz"

	cmd_1="ln -s ${fastq_1} ${sampleFile_1}"
	cmd_2="ln -s ${fastq_2} ${sampleFile_2}"
    echo "${cmd_1}" # debug
    echo "${cmd_2}" # debug

    eval ${cmd_1}
    eval ${cmd_2}
done
