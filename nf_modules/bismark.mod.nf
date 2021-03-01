nextflow.enable.dsl=2

// parameters passed in by specialised pipelines
params.singlecell = ''
params.pbat = false
params.unmapped = false
params.read_identity = ''

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
	    tuple val(name), path ("*bam"),        emit: bam
		path "*report.txt", emit: report
		// we always pass back the original name so we can use .join() later on, e.g. for bismark2bedGraph
		tuple val(name), path ("*unmapped_reads_1.fq.gz"), optional: true, emit: unmapped1
		tuple val(name), path ("*unmapped_reads_2.fq.gz"), optional: true, emit: unmapped2

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

		unmapped_name = ''	
			// add Genome build and aligner to output name
		if (params.read_identity == "1" || params.read_identity == "2"){
			// println ("FILENAME IS: $reads")
			if (params.read_identity == "1"){
				unmapped_name = name + "_unmapped_R1"
			}
			else{
				unmapped_name = name + "_unmapped_R2"
			}

			if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
				bismark_name = unmapped_name + "_" + params.genome["name"] + "_bismark_hisat2"
			}
			else{ // default is Bowtie 2
				bismark_name = unmapped_name + "_" + params.genome["name"] + "_bismark_bt2"
			}
		}
		else{
			if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
				bismark_name = name + "_" + params.genome["name"] + "_bismark_hisat2"
			}
			else{ // default is Bowtie 2
				bismark_name = name + "_" + params.genome["name"] + "_bismark_bt2"
			}
		}	
		// println ("Output basename: $bismark_name")
		
		"""
		module load bismark
		bismark --parallel $cores --basename $bismark_name $index $bismark_options $readString
		"""

}