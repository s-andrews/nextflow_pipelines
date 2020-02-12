nextflow.preview.dsl=2

params.deduplicate_bismark_args = ''

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
deduplicate_bismark_args = params.deduplicate_bismark_args.replaceAll(/'/,"")
println ("[BISMARK DEDUPLICATION MODULE, replaced] ARGS ARE: " + deduplicate_bismark_args)

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
	// --bedGraph --buffer 10G
	
	// if (reads instanceof List) {
	// 	readString = "-1 "+reads[0]+" -2 "+reads[1]
	// }
	// else {
	// 	readString = reads
	// }

	// index = " --genome " + params.genome["bismark"]
	
	"""
	module load bismark
	deduplicate_bismark ${deduplication_options} ${bam}
	"""

}