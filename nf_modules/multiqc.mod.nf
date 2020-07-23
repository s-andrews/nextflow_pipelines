nextflow.enable.dsl=2

process MULTIQC {
	
	label 'quadCore'

	// dynamic directive
	memory { 20.GB * task.attempt }  
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
	maxRetries 3

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