nextflow.enable.dsl=2

process BAMCOVERAGE{	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	input:
		path(bam)
		// We don't specifically use these but we need to import the .bai
		// files for the bam files as deeptools will look for these.
		path(index)
		val (outputdir)
		val (bamcoverage_args)
		val (verbose)

	output:
		path "*bw", emit: bw

	publishDir "$outputdir",
		mode: "link", overwrite: true

	
    script:
		bamcoverage_options = bamcoverage_args
		
		if (verbose){
			println ("[MODULE] BAMCOVERAGE ARGS: " + bamcoverage_args)
		}
		
		// TODO: Find more elegant way to strip file ending of input BAM file

		"""
		module load deeptools
		bamCoverage $bamcoverage_options -b $bam -of bigwig -o ${bam}.bw 
		rename .bam.bw .bw *
    	"""
		
	
}

