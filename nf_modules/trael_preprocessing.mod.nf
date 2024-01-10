nextflow.enable.dsl=2

/* process TRAEL_PREPROCESSING {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (trael_preprocessing_args)
		val (verbose)

	output:
		path "*.txt", optional: true, emit: stats 
		tuple val(name), path ("*UMIed*.fastq.gz"), emit: reads
        //path ("*UMIed*.fastq.gz"), emit: reads

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
} */

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

		// The TrAEL-seq preprocessing script works on all FastQ files in a folder	- not anymore - might want to switch back to that though	
		"""
		module load python
		/bi/apps/TrAELseq/latest/TrAEL-seq/TrAELseq_preprocessing_UMIplusBarcode.py	${reads}
		"""

}

process SORT_TRAEL_NAMES {

    input:
    path(reads)

    output:
    tuple env(sample_id_index), path(reads), emit:reads
    //stdout

   // when:
   // reads.size() > 0

    script:
    """
    sample_id_index=${reads}
    sample_id_index=\${sample_id_index%_R1.fastq.gz}
    """
}

