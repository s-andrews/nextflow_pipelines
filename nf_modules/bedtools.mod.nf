nextflow.enable.dsl=2
params.no_output = false

process BEDTOOLS_INTERSECT{	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

    input:
		path(bam)
        path(blacklist)
		val (outputdir)
		val (verbose)

	output:
		path "*.bam",     emit: bam

        publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

    script:
		
        bam_out = bam.baseName + "_clean.bam" 
		
		"""
		module load bedtools
        bedtools intersect -v -abam ${bam} -b ${blacklist} > ${bam_out}
		"""	
}