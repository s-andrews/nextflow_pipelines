nextflow.enable.dsl=2


// This is for the early versions of the TrAEL method did not incorporate the inline TrAEL barcodes. 
process TRAEL_PREPROCESSING {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (trael_preprocessing_args)
		val (verbose)

	output:
		path "*.txt", optional: true, emit: stats 
        path ("*UMIed*.fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] TrAEL-PREPROCESSING ARGS: " + trael_preprocessing_args)
		}

		// Run the TrAELseq_preprocessing script	
		"""
		module load python
		/bi/apps/TrAELseq/latest/TrAEL-seq/TrAELseq_preprocessing.py ${reads}		
		"""
} 


// For later versions of the TrAEL method do incorporate the inline TrAEL barcodes. 
process TRAEL_PREPROCESSING_INDEXING {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (trael_preprocessing_args)
		val (verbose)

	output:
		path "*.txt", optional: true, emit: stats 
        path ("*UMIed*.fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] TrAEL-PREPROCESSING ARGS: " + trael_preprocessing_args)
		}

		/* Run the TrAELseq_preprocessing_UMIplusBarcode script
			This splits the input fastq file by the 9 TrAEL barcodes and T or no T at position 13, 
			so produces 20 output fastq files (9 + unassigned)*2 for T and noT. 
		*/
		"""
		module load python
		/bi/apps/TrAELseq/latest/TrAEL-seq/TrAELseq_preprocessing_UMIplusBarcode.py	${reads}
		"""

}


/* The output from the TrAEL preprocessing scripts are a set of paths.
We need to create new sample names for these as other nf modules e.g. bowtie2 require them.
*/
process SORT_TRAEL_NAMES {

    input:
    path(reads)

    output:
    tuple env(sample_id_index), path(reads), emit:reads

    script:
    """
    sample_id_index=${reads}
    sample_id_index=\${sample_id_index%_R1.fastq.gz}
    """
}

