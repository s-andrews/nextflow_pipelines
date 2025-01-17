nextflow.enable.dsl=2
params.no_output = false
params.index = ""

// This is for the early versions of the TrAEL method did not incorporate the inline TrAEL barcodes. 
process INDEX_4BP_PREPROCESSING {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (index_preprocessing_args)
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

		supplied_indexes = params.index

		// Run the eccDNA_preprocessing script	
		"""
		module load python
		/bi/apps/TrAELseq/latest/TrAEL-seq/WGS_index_preprocess.py --input_file ${reads} --index ${supplied_indexes} --count_all		
		"""
} 
