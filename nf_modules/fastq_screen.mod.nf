nextflow.preview.dsl=2

params.fastq_screen_args = ''

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
fastq_screen_args = params.fastq_screen_args.replaceAll(/'/,"")
// println ("ARGS ARE [FASTQ SCREEN MODULE, replaced]: "+ fastq_screen_args + "\n")

process FASTQ_SCREEN {	
    input:
	    tuple val(name), path(reads)

	output:
	    path "*png", emit: png
	    path "*html", emit: html

    script:
	if (reads instanceof List) {
		reads = reads[0]
	}

	"""
	module load fastq_screen
	fastq_screen $fastq_screen_args	$reads
	"""

}