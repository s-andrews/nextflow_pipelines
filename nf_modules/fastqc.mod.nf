nextflow.preview.dsl=2

process FASTQC {	
	input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastqc_args)
		val (verbose)

	output:
	    path "*fastqc*"
	
	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:

		// We need to replace single quotes in the arguments so that they are not taken passed in as a single string
		// This will likely be fixed soon
		fastqc_args = fastqc_args.replaceAll(/'/,"")
		if (verbose){
			println ("[MODULE] FASTQC ARGS: "+ fastqc_args)
		}


	"""
	module load fastqc
	fastqc $fastqc_args -q -t 2 ${reads}
	"""

}