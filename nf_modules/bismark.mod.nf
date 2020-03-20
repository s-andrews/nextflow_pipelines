nextflow.preview.dsl=2

// parameters passed in by specialised pipelines
params.singlecell = ''
params.pbat = false


process BISMARK {
	
	label 'hugeMem'
	label 'multiCore'
	// label 'quadCore'
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bismark_args)
		val (verbose)

	output:
	    path "*bam",        emit: bam
		path "*report.txt", emit: report

	publishDir "$outputdir",
		mode: "link", overwrite: true


    script:
		cores = 1
		readString = ""

		if (verbose){
			println ("[MODULE] BISMARK ARGS: " + bismark_args)
		}

		// Options we add are
		bismark_options = bismark_args
		if (params.singlecell){
			bismark_options += " --non_directional "
		}
		
		if (params.pbat){
			bismark_options += " --pbat "
		}

		if (reads instanceof List) {
			readString = "-1 "+reads[0]+" -2 "+reads[1]
		}
		else {
			readString = reads
		}

		index = "--genome " + params.genome["bismark"]
		
		"""
		module load bismark
		bismark --parallel $cores $index $bismark_options $readString
		"""

}