nextflow.preview.dsl=2

process FASTQC {	
    input:
	    file query


	output:
	    file "*html" 

    script:
	"""
	module load fastqc
	fastqc -q ${query}
	"""

}