#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MERGE_BARCODES {

    tag "$name" // Adds name to job submission instead of (1), (2) etc.

    input: 
        tuple val(name), path(reads)
        val (outputdir)
        
    output:    
        tuple val(name), path ("*fastq"), emit: merged_fastq
        
    publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
        """
        cat ${reads}/*fastq > ${name}.fastq
        """
}
