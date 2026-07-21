nextflow.enable.dsl=2
params.single_end = false
params.no_output = false

process FASTQ_SCREEN {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	// label 'hugeMem'

	label 'multiCore'

	memory { 30.GB * task.attempt }
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
  	maxRetries 3

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastq_screen_args)
		val (verbose)
		val (bisulfite)

	output:
	    path "*png",  emit: png
	    path "*html", emit: html
		path "*txt",  emit: report

    script:

		if (verbose){
			println ("[MODULE] FASTQ SCREEN ARGS: "+ fastq_screen_args)
		}

		def screen_reads = reads
		if (params.single_end){
			// TODO: Add single-end parameter
		}
		else{
			// for paired-end files we only use Read 1 (as Read 2 tends to show the exact same thing)
			if (reads instanceof List) {
				screen_reads = reads[0]
			}
		}

		def fastq_screen_opts = fastq_screen_args
		if (bisulfite){
			// println("Setting --bisulfite")
			fastq_screen_opts += " --bisulfite "
			// println (fastq_screen_opts)
		}

	"""
	module load fastq_screen
	fastq_screen $fastq_screen_opts $screen_reads
	"""

}
