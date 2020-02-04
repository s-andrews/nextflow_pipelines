#!/usr/bin/env	nextflow

// Script parameters

params.files = ""
params.genome = ""
params.outdir = "."

params.index = "/bi/scratch/Genomes/Yeast/Saccharomyces_cerevisiae/R64-1-1/Saccharomyces_cerevisiae.R64-1-1"


file_ch = Channel.fromPath(params.files)

file_ch.into {fastqc_in_ch; trim_in_ch}

process runFastQC {
				publishDir params.outdir, mode: "copy", pattern: "*html"
				input:
				file query from fastqc_in_ch


				output:
				file "*html" into fastqc_out_ch

				"""
				module load fastqc
				fastqc ${query}
				"""

}

process runTrimGalore {
				publishDir params.outdir, mode: "copy", pattern: "*trimmed.fq.gz"
				publishDir params.outdir, mode: "copy", pattern: "*html"

				input:
				file query from trim_in_ch


				output:
				file "*trimmed*fq.gz" into trim_out_ch

				"""
				module load trim_galore
				module load fastqc
				trim_galore --gzip --fastqc ${query}
				"""

}


process runHisat2 {
				publishDir params.outdir, mode: "copy", pattern: "*bam"
				publishDir params.outdir, mode: "copy", pattern: "*_stats.txt"

				input:
				file query from trim_out_ch


				output:
				file "*bam" into hisat2_out_ch

				"""
				module load hisat2
				module load samtools
				hisat2 -x ${params.index} -p 2 -U ${query} | samtools view -bS -> ${query}_hisat2.bam
				"""

}




fastqc_out_ch.subscribe{ 
				println "Found report: $it"
}


hisat2_out_ch.subscribe{ 
				println "Found mapped file: $it"
}