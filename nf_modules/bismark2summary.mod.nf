nextflow.enable.dsl=2

process BISMARK2SUMMARY {
    
	input:
	    file (file)
		val (outputdir)
		val (bismark2summary_args)
		val (verbose)

	output:
	    path "*html",       emit: html
		path "*txt",        emit: report 

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
		// This is only a temporary workaround until Paolo has fixed the Nextflow bug.
		// https://github.com/nextflow-io/nextflow/issues/1519
		bismark2summary_args = bismark2summary_args.replaceAll(/'/,"")
		if (verbose){
			println ("[MODULE] BISMARK2SUMMARY ARGS: " + bismark2summary_args)
		}

		"""
		module load bismark
		bismark2summary
		"""

}