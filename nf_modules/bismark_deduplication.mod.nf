nextflow.preview.dsl=2

params.deduplicate_bismark_args = ''
params.verbose = false

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
deduplicate_bismark_args = params.deduplicate_bismark_args.replaceAll(/'/,"")
if (params.verbose){
	println ("[MODULE] BISMARK DEDUPLICATION ARGS: " + deduplicate_bismark_args)
}

process BISMARK_DEDUPLICATION {
	label 'bigMem'
		
    input:
	    path(bam)

	output:
		path "*report.txt", emit: report
		path "*bam",        emit: bam

    script:
	// readString = ""

	// Options we add are
	deduplication_options = deduplicate_bismark_args
	deduplication_options = deduplication_options + " --bam "

	"""
	module load bismark
	deduplicate_bismark ${deduplication_options} ${bam}
	"""

}