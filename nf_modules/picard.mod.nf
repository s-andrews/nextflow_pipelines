nextflow.enable.dsl=2
params.no_output = false

process PICARD_ADD_REPLACE{	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'hugeMem' // 80GB

    input:
		path(bam)
		val (outputdir)
		val (verbose)

	output:
		path "*.bam",     emit: bam


    script:
		
        base = bam.baseName
        bam_out = bam.baseName + "_rg.bam" 
		
		"""
		module load picard
		picard AddOrReplaceReadGroups --INPUT ${bam} --OUTPUT ${bam_out} --RGID ${base} --RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM ${base}
		"""	
}


process PICARD_DEDUP{	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'hugeMem' // 80GB

    input:
		path(bam)
		val (outputdir)
		val (verbose)

	output:
		path "*.bam",     emit: bam


    script:
		
        base = bam.baseName
        bam_out = bam.baseName + "_dedup.bam" 
		
		"""
		module load picard
		picard MarkDuplicates --INPUT ${bam} --OUTPUT ${bam_out} --M ${base}_metrics.txt --REMOVE_DUPLICATES true
		"""	
}

