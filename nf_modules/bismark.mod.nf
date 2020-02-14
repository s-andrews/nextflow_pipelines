nextflow.preview.dsl=2

params.bismark_args = ''
params.singlecell = ''
params.pbat = ''
params.verbose = false

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
bismark_args = params.bismark_args.replaceAll(/'/,"")
if (params.verbose){
	println ("[MODULE] BISMARK ARGS: " + bismark_args)
}

process BISMARK {
	
	label 'bigMem'
	label 'multiCore'
		
    input:
	    tuple val(name), path(reads)

	output:
	    path "*bam",        emit: bam
		path "*report.txt", emit: report

    script:
	cores = 2
	readString = ""

	// Options we add are
	bismark_options = bismark_args
	if (params.singlecell){
		bismark_options = bismark_options + " --non_directional "
	}
	
	if (params.pbat){
		bismark_options = bismark_options + " --pbat "
	}

	if (reads instanceof List) {
		readString = "-1 "+reads[0]+" -2 "+reads[1]
	}
	else {
		readString = reads
	}

	index = " --genome " + params.genome["bismark"]
	
	"""
	module load bismark
	bismark --parallel $cores $index $bismark_options $readString
	"""

}