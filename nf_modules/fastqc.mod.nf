nextflow.preview.dsl=2


params.fastqc_args = ''
//println ("ARGS ARE [FASTQC MODULE]: "+ params.fastqc_args + "\n")
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