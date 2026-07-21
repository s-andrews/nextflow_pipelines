nextflow.enable.dsl=2

process BISMARK_DEDUPLICATION {

	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	// dynamic directive to increase memory as required
	cpus 1
	memory { 20.GB * task.attempt }
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
  	maxRetries 5

	publishDir { outputdir },
		mode: "link", overwrite: true

    input:
	    tuple val(name), path(bam)
		val (outputdir)
		val (deduplicate_bismark_args)
		val (verbose)

	output:
		path "*report.txt", emit: report
		tuple val(name), path ("*bam"),        emit: bam

    script:
		if (verbose){
			println ("[MODULE] BISMARK DEDUPLICATION ARGS: " + deduplicate_bismark_args)
		}

		// Options we add are
		def deduplication_options = deduplicate_bismark_args
		deduplication_options += " --bam "

		"""
		module load bismark
		deduplicate_bismark ${deduplication_options} ${bam}
		"""

}
