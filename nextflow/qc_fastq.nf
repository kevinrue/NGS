#!/usr/bin/env nextflow

// params

params.project_root		= '${baseDir}' // in/out-put files relative to this path
params.fastqs			= 'testdata/*.fastq.{1,2}.gz' // relative to project_root
params.outDir			= 'qc'
params.logDir			= 'log'

params.fastq_version	= '0.11.5'
params.python_version	= '2.7.11'
params.multiqc_version	= '1.0.dev0'

// Inform user //

log.info " === PARAMS ==="
log.info "(params) project_root :\t\t${params.project_root}"
log.info "(params) fastqs :\t\t${params.fastqs}"
log.info "(params) outDir :\t\t${params.outDir}"
log.info "(params) logDir :\t\t${params.logDir}"
log.info " =--- versions ---="
log.info "(params) fastq_version :\t${params.fastq_version}"
log.info "(params) python_version :\t${params.python_version}"
log.info "(params) multiqc_version :\t${params.multiqc_version}"

log.info " === ENV ==="
log.info "(env) progsDir :\t\t${progsDir}"

// Post-process inputs //

fastq_exe = "${progsDir}/fastqc/${params.fastq_version}/FastQC/fastqc"
python_exe = "\$(which python)"
multiqc_exe = "${progsDir}/multiqc/${params.multiqc_version}/multiqc"

fastqFiles = Channel.fromPath( "${params.project_root}/${params.fastqs}" )

process fastqc {
	tag { fastqFile }
	echo true
	publishDir "${params.project_root}/${params.outDir}/fastqc", mode: 'copy'
	afterScript "cat .command.sh | grep -v '#' >> ${params.project_root}/${params.logDir}/cmd.sh"

	input:
	file fastqFile from fastqFiles

	output:
	file "*_fastqc.{zip,html}" into fastqc_results

	"""
	${fastq_exe} --format fastq --threads 10 ${fastqFile}
	"""

}


process multiqc {
	tag 'multiqc'
	echo true
	publishDir "${params.project_root}/${params.outDir}/multiqc", mode: 'copy'
	afterScript "cat .command.sh | grep -v '#' >> ${params.project_root}/${params.logDir}/cmd.sh"

	// input: used to link all relevant input files for MultiQC in the workdir
	// collect() used to run once this step no matter how many inputs
	input:
	file ('qc/fastqc/*') from fastqc_results.collect()

	output:
	file 'multiqc_report.html'
    file 'multiqc_data'

    """
    module load python/${params.python_version}
    $multiqc_exe --force .
    """

}


