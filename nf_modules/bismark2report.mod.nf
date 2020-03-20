nextflow.preview.dsl=2

process BISMARK2REPORT {
	
    input:
	    file (file)
		val (outputdir)
		val (bismark2report_args)
		val (verbose)

	output:
	    path "*html",       emit: html
		
	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		if (verbose){
			println ("[MODULE] BISMARK2REPORT ARGS: " + bismark2report_args)
		}

		"""
		module load bismark
		bismark2report
		"""

}