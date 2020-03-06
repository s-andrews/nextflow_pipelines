nextflow.preview.dsl=2

process MULTIQC {
	
	// label 'bigMem'
	// label 'multiCore'

    input:
	    file (file)
		val (outputdir)
		val (multiqc_args)
		val (verbose)

	output:
	    path "*html",       emit: html
		// path "*stats.txt", emit: stats 

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
		// This is only a temporary workaround until Paolo has fixed the Nextflow bug.
		// https://github.com/nextflow-io/nextflow/issues/1519
		multiqc_args = multiqc_args.replaceAll(/'/,"")
		if (verbose){
			println ("[MODULE] MULTIQC ARGS: " + multiqc_args)
		}

		"""
		module load multiqc
		multiqc . -x work
		"""

}