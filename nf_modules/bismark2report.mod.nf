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
		// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
		// This is only a temporary workaround until Paolo has fixed the Nextflow bug.
		// https://github.com/nextflow-io/nextflow/issues/1519
		
		bismark2report_args = bismark2report_args.replaceAll(/'/,"")
		if (verbose){
			println ("[MODULE] BISMARK2REPORT ARGS: " + bismark2report_args)
		}

		"""
		module load bismark
		bismark2report
		"""

}