nextflow.enable.dsl=2
params.no_output = false
//params.index = ""

// This is for the early versions of the TrAEL method did not incorporate the inline TrAEL barcodes. 
process INDEX_4BP_PREPROCESSING {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (index_preprocessing_args)
		val (supplied_indexes)
		val (verbose)

	output:
		path "*.log", emit: log 
        //tuple val(name), path ("*.fastq.gz"), emit: reads
		path ("*.fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	script:
		if (verbose){
			println ("[MODULE] WGS_4bpIndex PREPROCESSING ARGS: " + index_preprocessing_args)
		}

		//supplied_indexes = params.index

		// Run the eccDNA_preprocessing script	
		"""
		module load python
		/bi/apps/TrAELseq/latest/TrAEL-seq/WGS_index_preprocess.py --input_file ${reads} --index ${supplied_indexes} --count_all		
		"""
} 

process INDEX_4BP_PREPROCESSING_PE {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (index_preprocessing_args)
		val (verbose)

	output:
		path "*.log", emit: log 
        //tuple val(name), list ("*.fastq.gz"), emit: reads
		tuple val(name), path("*.fastq.gz"), emit: reads
		//path ("*.fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	script:
		if (verbose){
			println ("[MODULE] WGS_4bpIndex PREPROCESSING ARGS: " + index_preprocessing_args)
		}

		supplied_indexes = params.index
		r1 = reads[0]
		r2 = reads[1]

		println ("R1: " + r1)
		println ("R2: " + r2)

		// Run the preprocessing script	
		"""
		module load python
		/bi/apps/TrAELseq/latest/TrAEL-seq/WGS_index_preprocess_PE.py --input_fileR1 ${r1} --input_fileR2 ${r2} --index ${supplied_indexes} --count_all		
		"""
} 

process INDEX_4BP_PREPROCESSING_PE_TEST {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (index_preprocessing_args)
		val (verbose)

	output:
		path "*.log", emit: log 
        tuple val(name), path("*.fastq.gz"), emit: reads
		//path ("*.fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	script:
		if (verbose){
			println ("[MODULE] WGS_4bpIndex PREPROCESSING ARGS: " + index_preprocessing_args)
		}

		supplied_indexes = params.index
		r1 = reads[0]
		r2 = reads[1]

		println ("R1: " + r1)
		println ("R2: " + r2)

		// Run the preprocessing script	
		"""
		module load python
		/bi/apps/TrAELseq/latest/TrAEL-seq/WGS_index_preprocess_PE.py --input_fileR1 ${r1} --input_fileR2 ${r2} --index ${supplied_indexes} --count_all		
		"""
} 




def get_paired_names(fileList){

    println ("fileList = " + fileList)

    baseNames = [:]

    for (String s : fileList) {
        matcher = s =~ /^(.*)_(R?[1234]).(fastq|fq).gz$/
        println ("trying matches")
        //println (matcher[0])
        if (matcher.matches()) {
            if (! baseNames.containsKey(matcher[0][1])) {
                baseNames[matcher[0][1]] = []
            }
            baseNames[matcher[0][1]].add(matcher[0][2])
        }
    }

    patterns = []
    println("found anything in baseNames?")

    for (s in baseNames) {
        println("baseNames loop:")
        println (s)
        pattern = s.key+"_{"+s.value.join(",")+"}.{fastq,fq}.gz"
        patterns.add(pattern)
        println("#!!!!=============================!")
        println(pattern)
    }

    //file_ch = Channel.fromFilePairs(patterns, size:-1)
    //return(file_ch)
    return(patterns)
}