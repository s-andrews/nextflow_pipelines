nextflow.preview.dsl=2

process FASTQC {	
    input:
	    tuple val(name), path(reads)


	output:
	    path "*html" 

    script:
	"""
	module load fastqc
	fastqc -q -t 2 ${reads}
	"""

}