#!/usr/bin/env nextflow

// params

params.fastqs			= "${baseDir}/raw/*.fastq.{1,2}.gz" // TODO: error if no file?
params.outDir			= "${baseDir}/qc"

params.fastq_version	= '0.11.5'
params.python_version	= '2.7.11'
params.multiqc_version	= '1.0.dev0'

// Inform user //

log.info " === PARAMS ==="
log.info "(params) fastqs :\t\t${params.fastqs}"
log.info "(params) outDir :\t\t${params.outDir}"
log.info " =--- versions ---="
log.info "(params) fastq_version :\t${params.fastq_version}"
log.info "(params) python_version :\t${params.python_version}"
log.info "(params) multiqc_version :\t${params.multiqc_version}"

log.info " === ENV ==="
log.info "(env) progsDir :\t\t${progsDir}"

// Post-process inputs //

fastq_exe = "${progsDir}/fastqc/${params.fastq_version}/FastQC/fastqc" // NOTE: Error if fastq_exe does not exist ?
python_exe = "\$(which python)" // NOTE: Error if fastq_exe does not exist ?
multiqc_exe = "${progsDir}/multiqc/${params.multiqc_version}/multiqc" // NOTE: Error if multiqc_exe does not exist ?

fastqFiles = Channel.fromPath( params.fastqs )


process fastqc {
	tag { fastqFile }
	echo true
	publishDir "${params.outDir}/fastqc", mode: 'copy'

	input:
	file fastqFile from fastqFiles

	output:
	file "*_fastqc.{zip,html}" into fastqc_results

	"""
	echo ${fastq_exe}
	${fastq_exe} --format fastq --threads 10 ${fastqFile}
	"""

}


process multiqc {
	tag 'multiqc'
	echo true
	publishDir "${params.outDir}/multiqc", mode: 'copy'

	// input: used to link all relevant input files for MultiQC in the workdir
	// collect() used to run once this step no matter how many inputs
	input:
	file ('qc/fastqc/*') from fastqc_results.collect()

	output:
	file '*multiqc_report.html'
    file '*multiqc_data'

    """
    module load python/${params.python_version}
    echo ${multiqc_exe}
    $multiqc_exe -f .
    """

}


