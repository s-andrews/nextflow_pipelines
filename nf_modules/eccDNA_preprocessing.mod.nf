nextflow.enable.dsl=2
params.no_output = false


// This is for the early versions of the TrAEL method did not incorporate the inline TrAEL barcodes. 
process ECCDNA_PREPROCESSING {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (ecc_preprocessing_args)
		val (verbose)

	output:
		//path "*.txt", optional: true, emit: stats 
        tuple val(name), path ("*UMIed*.fastq.gz"), emit: reads

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	script:
		if (verbose){
			println ("[MODULE] eccDNA-PREPROCESSING ARGS: " + ecc_preprocessing_args)
		}

		// Run the eccDNA_preprocessing script	
		"""
		module load python
		/bi/apps/TrAELseq/latest/TrAEL-seq/eccDNAseq_preprocessing.py ${reads}		
		"""
} 
