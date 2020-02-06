nextflow.preview.dsl=2

process TRIM_GALORE {	
    input:
	    tuple val(name), path(reads)


	output:
	    tuple val(name), path ("*fq.gz"), emit: reads

    script:

	pairedString = ""
	if (reads instanceof List) {
		pairedString = "--paired"
	}

	"""
	module load trim_galore
	trim_galore ${pairedString} ${reads}
	"""

}