nextflow.enable.dsl=2
params.no_output = false

process SAMTOOLS_SORT{	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	input:
		path(bam)
		val (outputdir)
		val (samtools_sort_args)
		val (verbose)

	output:
		// path "*report.txt", emit: report
		path "*bam",        emit: bam

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	
    script:
		samtools_sort_options = samtools_sort_args
		
		if (verbose){
			println ("[MODULE] SAMTOOLS SORT ARGS: " + samtools_sort_args)
		}
		
		// TODO: Find more elegant way to strip file ending of input BAM file

		"""
		module load samtools
		samtools sort $samtools_sort_options $bam -o ${bam}_sorted.bam 
		rename .bam_sorted _sorted *
    	"""
		
	
}

process SAMTOOLS_INDEX{	
    
	tag "$bam"     // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	input:
		path(bam)
		val (outputdir)
		val (samtools_index_args)
		val (verbose)

	output:
		path "*.bai",     emit: bai
    	
	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

    script:
		samtools_index_options = samtools_index_args
		
		if (verbose){
			println ("[MODULE] SAMTOOLS INDEX ARGS: " + samtools_index_args)
		}
		
		"""
		module load samtools
		samtools index $samtools_index_options $bam
		"""
		
	
}