nextflow.enable.dsl=2
params.bisulfite = ''
params.single_end = false

process FASTQ_SCREEN {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	// label 'hugeMem'
	
	label 'multiCore'
	
	memory { 30.GB * task.attempt }  
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
  	maxRetries 3
  	
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastq_screen_args)
		val (verbose)

	output:
	    path "*png",  emit: png
	    path "*html", emit: html
		path "*txt",  emit: report

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:

		if (verbose){
			println ("[MODULE] FASTQ SCREEN ARGS: "+ fastq_screen_args)
		}

		if (params.single_end){
			// TODO: Add single-end parameter
		}
		else{
			// for paired-end files we only use Read 1 (as Read 2 tends to show the exact same thing)
			if (reads instanceof List) {
				reads = reads[0]
			}
		}
		if (params.bisulfite){
			// println("Setting --bisulfite")
			fastq_screen_args += " --bisulfite "
			// println (fastq_screen_args)
		}	

	"""
	module load fastq_screen
	fastq_screen $fastq_screen_args $reads
	"""

}