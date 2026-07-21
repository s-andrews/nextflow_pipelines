nextflow.enable.dsl=2

process RESTRAININGORDER {
	
	tag "$bam" // Adds the file name to job submission

	label 'hugeMem'

	publishDir { outputdir },
		mode: "link", overwrite: true

    input:
	    tuple val(name), path(bam)
		val (outputdir)
		val (reStrainingOrder_args)
		val (verbose)

	output:
	    path "*.html",  emit: html
		path "*.txt", emit: stats

	script:
		if (verbose){
			println ("[MODULE] reStrainingOrder ARGS: " + reStrainingOrder_args)
		}

		// Options we add are
		def reStraining_options = reStrainingOrder_args
		reStraining_options +=  "--snp_file /bi/apps/reStrainingOrder/MGPv5_SNP_matrix_chr1.txt.gz "
	
		"""
		module load reStrainingOrder
		reStrainingOrder $reStraining_options $bam
		"""

}