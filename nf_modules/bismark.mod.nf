nextflow.enable.dsl=2

// parameters passed in by specialised pipelines
params.singlecell = ''
params.pbat = false
params.unmapped = false

process BISMARK {
	

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	// TODO: Fix memory requirements, probably with error handling...
	//label 'hugeMem'

	cpus { 5 }
  	memory { 20.GB * task.attempt }  
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
  	maxRetries 3
	
	// label 'mem40G'
	// label 'multiCore'
	// label 'quadCore'
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bismark_args)
		val (verbose)

	output:
	    path "*bam",        emit: bam
		path "*report.txt", emit: report
		val name,           emit: sample_name
		tuple val(unmapped_1_name), path ("*unmapped_reads_1.fq.gz"), optional: true, emit: unmapped1
		tuple val(unmapped_2_name), path ("*unmapped_reads_2.fq.gz"), optional: true, emit: unmapped2

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		cores = 1
		readString = ""
		// println ("WHO KEEPS CHANGING THE NAME 1: $name")

		if (verbose){
			println ("[MODULE] BISMARK ARGS: " + bismark_args)
		}

		// Options we add are
		bismark_options = bismark_args
		if (params.singlecell){
			bismark_options += " --non_directional "
		}
		else{
		
		}
		
		unmapped_1_name = ''
		unmapped_2_name = ''
		
		if (params.unmapped){
			bismark_options += " --unmapped "
			unmapped_1_name = name + "_unmapped_R1"
			unmapped_2_name = name + "_unmapped_R2"
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

		// add Genome build and aligner to output name
		if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
			bismark_name = name + "_" + params.genome["name"] + "_bismark_hisat2"
		}
		else{ // default is Bowtie 2
			bismark_name = name + "_" + params.genome["name"] + "_bismark_bt2"
		}
		// println ("Output basename: $bismark_name")
		// println ("WHO KEEPS CHANGING THE NAME 2: $name")

		"""
		module load bismark
		bismark --parallel $cores --basename $bismark_name $index $bismark_options $readString
		"""

}