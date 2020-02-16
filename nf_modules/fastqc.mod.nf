nextflow.preview.dsl=2

params.fastqc_args = ''
params.verbose = false

// We need to replace single quotes in the arguments so that they are not taken passed in as a single string'
fastqc_args = params.fastqc_args.replaceAll(/'/,"")
if (params.verbose){
	println ("[MODULE] FASTQC ARGS: "+ fastqc_args)
}

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