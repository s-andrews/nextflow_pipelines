nextflow.preview.dsl=2

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
	fastq_screen $reads
	"""

}