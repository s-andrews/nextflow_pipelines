nextflow.enable.dsl=2

// let's use an empty default prefix which can be set for each of the pipelines invoking MultiQC
// using --prefix, e.g. "--prefix lane8075_L001_" so that we can copy files using 'copy_back_files'. 
params.prefix = "" 

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
		
	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:

		
		if (verbose){
			println ("[MODULE] MULTIQC ARGS: " + multiqc_args)
		}
	
		"""
		module load multiqc
		multiqc $multiqc_args -x work --filename ${params.prefix}multiqc_report.html .
		"""
}