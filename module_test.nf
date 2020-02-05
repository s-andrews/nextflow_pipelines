#!/usr/bin/env nextflow

params.files = ""
params.outdir = "output"

// Enable modules
nextflow.preview.dsl=2
include './nf_modules/fastqc.mod.nf'
include './nf_modules/fastq_screen.mod.nf'

file_ch = Channel.fromPath(params.files) 


workflow {

    main:
        FASTQC(file_ch)
        FASTQ_SCREEN(file_ch)

        FASTQ_SCREEN.out.html.subscribe{
         println("Found fastq screen report $it")
        }

    publish:
        FASTQC.out to: params.outdir
        FASTQ_SCREEN.out.html to: params.outdir

}


file_ch.subscribe{
    println("Found $it")
}
