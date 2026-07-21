nextflow.enable.dsl=2

process BISMARK2REPORT {

	publishDir { outputdir },
		mode: "link", overwrite: true

    input:
	    file (file)
		val (outputdir)
		val (bismark2report_args)
		val (verbose)

	output:
	    path "*html",       emit: html

    script:
		if (verbose){
			println ("[MODULE] BISMARK2REPORT ARGS: " + bismark2report_args)
		}

		"""
		module load bismark
		bismark2report
		"""

}