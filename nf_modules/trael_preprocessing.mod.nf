nextflow.enable.dsl=2

process TRAEL_PREPROCESSING {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (trael_preprocessing_args)
		val (verbose)

	output:
		path "*.txt", optional: true, emit: stats 
		tuple val(name), path ("*UMIed*.fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] TrAEL-PREPROCESSING ARGS: " + trael_preprocessing_args)
		}

		// The TrAEL-seq preprocessing script works on all FastQ files in a folder		
		"""
		module load python
		/bi/apps/TrAELseq/TrAEL-seq/TrAELseq_preprocessing.py 		
		"""

}

process TRAEL_PREPROCESSING_INDEXING {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (trael_preprocessing_args)
		val (verbose)

	output:
		path "*.txt", optional: true, emit: stats 
		tuple val(name), path ("*UMIed*.fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] TrAEL-PREPROCESSING ARGS: " + trael_preprocessing_args)
		}

		// The TrAEL-seq preprocessing script works on all FastQ files in a folder		
		"""
		module load python
		/bi/apps/TrAELseq/TrAEL-seq/TrAELseq_preprocessing_UMIplusBarcode.py		
		"""

}

process TRAEL_SEQDEDUP {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (trael_preprocessing_args)
		val (verbose)

	output:
		path "*.txt", optional: true, emit: stats 
		tuple val(name), path ("*seqDedup*fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] TrAEL-PREPROCESSING ARGS: " + trael_preprocessing_args)
		}

		// The TrAEL-seq preprocessing script works on all FastQ files in a folder		
		"""
		module load python
		/bi/apps/TrAELseq/TrAEL-seq/TrAELseq_sequence_based_deduplication.py	
		"""

}

