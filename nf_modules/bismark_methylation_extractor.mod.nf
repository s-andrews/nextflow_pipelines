nextflow.enable.dsl=2

process BISMARK_METHYLATION_EXTRACTOR {
	label 'bigMem'          // 20G
	label 'quadCore'        // 4 cores

	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	publishDir { outputdir },
		mode: "link", overwrite: true

    input:
	    tuple val (name), path(bam)
		val (outputdir)
		val (bismark_methylation_extractor_args)
		val (verbose)
		val (singlecell)
		val (rrbs)
		val (pbat)
		val (nonCG)
		val (emseq)

	output:
	    tuple val (name), path ("CpG*"),        emit: context_files_CG
		path "CH*",                             emit: context_files_nonCG
		path "*report.txt",                     emit: report
		path "*M-bias.txt",                     emit: mbias
		path "*cov.gz",                         emit: coverage

	script:

		if (verbose){
			println ("[MODULE] BISMARK METHYLATION EXTRACTOR ARGS: " + bismark_methylation_extractor_args)
		}

		def cores = 4

		// Options we add are
		def methXtract_options = bismark_methylation_extractor_args + " --gzip "

		if (singlecell){
			// println ("FLAG SINGLE CELL SPECIFIED: PROCESSING ACCORDINGLY")
		}

		if (nonCG){
			if (verbose){
				println ("FLAG nonCG specified: adding flag --CX ")
			}
			methXtract_options +=  " --CX "
		}

		def isPE = isPairedEnd(bam, verbose)
		if (isPE){
			// not perform any ignoring behaviour for RRBS or single-cell libraries
			if (!rrbs && !singlecell && !pbat && !emseq){
				// default ignore parameters for paired-end libraries
				methXtract_options +=  " --ignore_r2 2 "
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


def isPairedEnd(bamfile, verbose) {

	// need to transform the nextflow.processor.TaskPath object to String
	bamfile = bamfile.toString()
	if (verbose){
		println ("Processing file: " + bamfile)
	}

	if (bamfile =~ /_pe/){
		if (verbose){
			println ("File is paired-end!")
		}
		return true
	}
	else{
	 	if (verbose){
			 println ("File is single-end")
		 }
		return false
	}
}
