#!/usr/bin/env nextflow

// params
//   choices:
//     - species: { Mus_musculus, Homo_sapiens, ... }
//     - source: { UCSC, Ensembl }
//     - release: { UCSC, Ensembl }

log.info " === PARAMS ==="
params.species			= 'Homo_sapiens'
params.source			= 'UCSC'
params.release			= 'hg19'

log.info " === ENV ==="
log.info "(env) progsDir :\t\t${progsDir}"
log.info "(env) igenomesDir :\t\t${igenomesDir}"

// Post-process inputs //

reference_in = "${igenomesDir}/${params.species}/${params.source}/${params.release}/Sequence/Chromosomes/*.fa"
ht2_base = "${igenomesDir}/${params.species}/${params.source}/${params.release}/Sequence/Hisat2Index"


process hisat2build {
	tag 'hisat2-build'
	echo true

	script:
	"""
	echo hisat2-build \
		[options]*\
		${reference_in} \
		<ht2_base>
	"""
}