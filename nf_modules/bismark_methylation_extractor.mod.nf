nextflow.preview.dsl=2

params.bismark_methylation_extractor_args = ''
params.singlecell = false
params.rrbs = false
params.verbose = false

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
bismark_methylation_extractor_args = params.bismark_methylation_extractor_args.replaceAll(/'/,"")
if (params.verbose){
	println ("[MODULE] BISMARK METHYLATION EXTRACTOR ARGS: " + bismark_methylation_extractor_args)
}

process BISMARK_METHYLATION_EXTRACTOR {
	label 'bigMem'
	label 'multiCore'
		
    input:
	    path(bam)

	output:
	    path "C*",          emit: context_files
		path "*report.txt", emit: report
		path "*M-bias.txt", emit: mbias
		path "*cov.gz",     emit: coverage
		
    script:
	cores = 2

	// Options we add are
	methXtract_options = bismark_methylation_extractor_args + " --gzip "
	
	if (params.singlecell){
		// println ("FLAG SINGLE CELL SPECIFIED: PROCESSING ACCORDINGLY")
	}

	isPE = isPairedEnd(bam)
	if (isPE){
		// not perform any ignoring behaviour for RRBS or single-cell libraries
		if (!params.rrbs && !params.singlecell){
			// default ignore parameters for paired-end libraries
			methXtract_options = methXtract_options + " --ignore_r2 2 "
		}
	}
	else{
		// println("File seems to be single-end")
	}

	// println ("Now running command: bismark_methylation_extractor -parallel ${cores} ${methXtract_options} ${bam}")
	"""
	module load bismark
	bismark_methylation_extractor --bedGraph --buffer 10G -parallel ${cores} ${methXtract_options} ${bam}
	"""

}


def isPairedEnd(bamfile) {

	// need to transform the nextflow.processor.TaskPath object to String
	bamfile = bamfile.toString()
	if (params.verbose){
		println ("Processing file: " + bamfile)
	}
	
	if (bamfile =~ /_pe/){
		if (params.verbose){
			println ("File is paired-end!")
		}
		return true
	}
	else{
	 	if (params.verbose){
			 println ("File is single-end")
		 }
		return false
	}
}