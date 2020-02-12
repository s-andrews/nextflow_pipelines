nextflow.preview.dsl=2

params.bismark_methylation_extractor_args = ''

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
bismark_methylation_extractor_args = params.bismark_methylation_extractor_args.replaceAll(/'/,"")
println ("[BISMARK METHYLATION EXTRACTOR MODULE, replaced] ARGS ARE: " + bismark_methylation_extractor_args)

process BISMARK_METHYLATION_EXTRACTOR {
	label 'bigMem'
	label 'multiCore'
		
    input:
	    path(reads)

	output:
	    path "C*",          emit: context_files
		path "*report.txt", emit: report
		path "*M-bias.txt",      emit: mbias

    script:
	cores = 2
	// readString = ""

	// Options we add are
	methXtract_options = bismark_methylation_extractor_args
	methXtract_options = methXtract_options + " --gzip "
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
	bismark_methylation_extractor --gzip -parallel ${cores} ${methXtract_options} ${reads}
	"""

}