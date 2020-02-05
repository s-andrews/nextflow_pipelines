nextflow.preview.dsl=2

process FASTQ_SCREEN {	
    input:
	    path query


	output:
	    path "*png", emit: png
	    path "*html", emit: html

    script:
	"""
	module load fastq_screen
	fastq_screen ${query}
	"""

}