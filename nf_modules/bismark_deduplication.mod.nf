nextflow.preview.dsl=2

process BISMARK_DEDUPLICATION {
	label 'hugeMem'
	// consider dynamic directive to increase memory
	// memory { 2.GB * task.attempt }
    // time { 1.hour * task.attempt }
    // errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    // maxRetries 3

    input:
	    path(bam)
		val (outputdir)
		val (deduplicate_bismark_args)
		val (verbose)

	output:
		path "*report.txt", emit: report
		path "*bam",        emit: bam

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		if (verbose){
			println ("[MODULE] BISMARK DEDUPLICATION ARGS: " + deduplicate_bismark_args)
		}

		// Options we add are
		deduplication_options = deduplicate_bismark_args
		deduplication_options += " --bam "

		"""
		module load bismark
		deduplicate_bismark ${deduplication_options} ${bam}
		"""

}