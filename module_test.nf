#!/usr/bin/env nextflow

params.files = ""
params.outdir = "output"

// Enable modules
nextflow.preview.dsl=2
include './nf_modules/fastqc.mod.nf'

file_ch = Channel.fromPath(params.files) 


workflow {

    main:
        FASTQC(file_ch)    

        FASTQC.out.subscribe{
         println("Found fastqc report $it")
        }

    publish:
        FASTQC.out to: params.outdir

}


file_ch.subscribe{
    println("Found $it")
}
