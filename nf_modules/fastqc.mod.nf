nextflow.preview.dsl=2


params.fastqc_args = ''
// We need to replace single quotes in the arguments so that they are not taken passed in as a single string'
fastqc_args = params.fastqc_args.replaceAll(/'/,"")
// println ("ARGS ARE [FASTQC MODULE, replaced]: "+ fastqc_args + "\n")

process FASTQC {	
    input:
	    tuple val(name), path(reads)

	output:
	    path "*html" 
	
	script:
	"""
	module load fastqc
	fastqc $fastqc_args -q -t 2 ${reads}
	"""

}