#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process PYCHOPPER {

    tag "$name" // Adds name to job submission instead of (1), (2) etc.

    input: 
        tuple val(name), path(reads)
        val (outputdir)
        val(projectname)
             
    output:    
        path "*", emit: all
        
    publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
        """
        wc -l ${reads} > ${projectname}_${name}.log
        module load python
        module load hmmer
        module load ssub
        cdna_classifier.py -r report_${name}_${projectname}.pdf -u unclassified_${name}_${projectname}.fastq -w rescued_${name}.fastq -S stats_${name}_${projectname} ${reads} full_length_barcode_${name}.fastq
        """
}

