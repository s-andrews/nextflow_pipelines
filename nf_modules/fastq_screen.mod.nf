nextflow.preview.dsl=2

params.fastq_screen_args = ''
params.verbose = false

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
fastq_screen_args = params.fastq_screen_args.replaceAll(/'/,"")
if (params.verbose){
	println ("[MODULE] FASTQ SCREEN ARGS: "+ fastq_screen_args)
}

process FASTQ_SCREEN {	
    input:
	    tuple val(name), path(reads)

	output:
	    path "*png",  emit: png
	    path "*html", emit: html
		path "*txt",  emit: report

    script:
	if (reads instanceof List) {
		reads = reads[0]
	}

	"""
	module load fastq_screen
	fastq_screen $fastq_screen_args	$reads
	"""

}