nextflow.enable.dsl=2

process MULTIQC {
	
	// label 'bigMem'
	// label 'multiCore'

    input:
	    path (file)
		val (outputdir)
		val (multiqc_args)
		val (verbose)

	output:
	    path "*html",       emit: html
		// path "*stats.txt", emit: stats 

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		
		if (verbose){
			println ("[MODULE] MULTIQC ARGS: " + multiqc_args)
		}

		"""
		module load multiqc
		multiqc $multiqc_args -x work .
		"""

}