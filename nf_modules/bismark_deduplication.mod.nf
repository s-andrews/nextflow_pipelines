nextflow.preview.dsl=2

process BISMARK_DEDUPLICATION {
	label 'bigMem'
		
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

		// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
		deduplicate_bismark_args = deduplicate_bismark_args.replaceAll(/'/,"")
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