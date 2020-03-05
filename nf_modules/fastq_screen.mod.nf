nextflow.preview.dsl=2
params.bisulfite = ''

process FASTQ_SCREEN {
	label 'bigMem'
	label 'multiCore'

    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastq_screen_args)
		val (verbose)

	output:
	    path "*png",  emit: png
	    path "*html", emit: html
		path "*txt",  emit: report

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:

		// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
		// This is only a temporary workaround until Paolo has fixed the Nextflow bug.
		// https://github.com/nextflow-io/nextflow/issues/1519
		fastq_screen_args = fastq_screen_args.replaceAll(/'/,"")

		if (verbose){
			println ("[MODULE] FASTQ SCREEN ARGS: "+ fastq_screen_args)
		}

		if (reads instanceof List) {
			reads = reads[0]
		}

	"""
	module load fastq_screen
	fastq_screen $params.bisulfite $fastq_screen_args $reads
	"""

}