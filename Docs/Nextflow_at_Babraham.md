# Nextflow Pipelines at the Babraham Institute

We are currently transitioning from our previous pipelining system (Clusterflow) to a new one based on [Nextflow](https://www.nextflow.io/docs/latest/index.html). We offer some preconfigured pipelines that generally discriminate between two different modes of operation: 

- data type specific, multi-step pipelines
- single program pipelines (formerly known as modules)

These pipelines are curated by the Babraham Bioinformatics Group, but you are of course welcome to write and use your own additional pipelines. If you need help getting started with Nextflow, please come and see any member of the Bioinformatics group who shall be happy to help.

## Pipelines:

Pipelines are supposed to work in a stream-lined and reproducible way every time they are run, and are designed so that users don't have to worry about specifying any of the plethora of options each tool provides. To this end, we try to run the individual programs of a pipeline with a pre-configured set of parameters that are (we find) sensible for the specific data type.

#### List of current pipelines:

##### nf_qc
    FastQC, FastQ Screen
##### nf_rnaseq
    FastQC, FastQ Screen, Trim Galore, trimmed FastQC, HISAT2
##### nf_chipseq
    FastQC, FastQ Screen, Trim Galore, trimmed FastQC, Bowtie2
##### nf_bisulfite_WGBS
    FastQC, FastQ Screen, Trim Galore, trimmed FastQC, Bismark, deduplicate, methylation extract, coverage file
##### nf_bisulfite_scBSseq
    FastQC, FastQ Screen, Trim Galore (5' clip), trimmed FastQC, Bismark, deduplicate, methylation extract, coverage file
##### nf_bisulfite_RRBS
    FastQC, FastQ Screen, Trim Galore, trimmed FastQC, trimmed FastQC, Bismark, methylation extract, coverage file


## Single Prorgam Pipelines:

#### List of current single program pipelines:
- nf_fastqc
- nf_fastq_screen
- nf_trim_galore
- nf_trim_galore_speciality
- nf_bowtie2
- nf_hisat2
- nf_bismark

In addition to the default parameters, each pipeline accepts a tool-specific add
All pre-configured pipelines take one additional argument, which has to be exactly in the following form to work:

```
--toolname_args "'--additional_option value --extra_flag etc.'"
```

So as an example, you could run specific trimming in Trim Galore like so:

```
--trim_galore_args "'--clip_r1 10 --clip_r2 10 --nextera'"
```

## In a few more details...

The worklows we are going to use here are based on the modules system introduced with [DSL2](https://www.nextflow.io/docs/latest/dsl2.html) (see also https://en.wikipedia.org/wiki/Domain-specific_language). In essence, we need a module for each program/tool, and then a separte workflow that defines the different steps that are carried out for given input files. Here is an example of the current RNA-seq workflow, which does the following consecutive steps for each FastQ file (single-end), or file  pair (paired-end):

- run FastQC or raw FastQ(s)
- run FastQ Screen species screen
- run Trim Galore to remove adapters and low quality base calls
- run FastQC again on the trimmed files
- take the trimmed FastQ files and align them to a genome using HISAT2

All output will be written to the working directory.

#### Example of an RNA-seq workflow
```nextflow
#!/usr/bin/env nextflow

// Enable modules
nextflow.preview.dsl=2

params.outdir = "."
params.genome = ""
params.verbose = false

params.fastqc_args = ''
params.fastq_screen_args = ''
params.trim_galore_args = ''
params.hisat2_args = ''

if (params.verbose){
    println ("[WORKFLOW] FASTQC ARGS: "           + params.fastqc_args)
    println ("[WORKFLOW] FASTQ SCREEN ARGS ARE: " + params.fastq_screen_args)
    println ("[WORKFLOW] TRIM GALORE ARGS: "      + params.trim_galore_args)
    println ("[WORKFLOW] HISAT2 ARGS ARE: "       + params.hisat2_args)
}

include './nf_modules/files.mod.nf'
include './nf_modules/genomes.mod.nf'
include './nf_modules/fastqc.mod.nf'                         params(fastqc_args: params.fastqc_args, verbose: params.verbose)
include  FASTQC as FASTQC2 from './nf_modules/fastqc.mod.nf' params(fastqc_args: params.fastqc_args, verbose: params.verbose)

include './nf_modules/fastq_screen.mod.nf'                   params(fastq_screen_args: params.fastq_screen_args, verbose: params.verbose)
include './nf_modules/trim_galore.mod.nf'                    params(trim_galore_args:  params.trim_galore_args, verbose: params.verbose)

include './nf_modules/hisat2.mod.nf' params(genome: getGenome(params.genome), hisat2_args:  params.hisat2_args, verbose: params.verbose)

file_ch = makeFilesChannel(args)

genome = getGenome(params.genome)


workflow {

    main:
        FASTQC(file_ch)
        FASTQ_SCREEN(file_ch)
        TRIM_GALORE(file_ch)
        FASTQC2(TRIM_GALORE.out.reads)
        HISAT2(TRIM_GALORE.out.reads)

    publish:
        FASTQC.out              to: params.outdir, mode: "link", overwrite: true
        FASTQ_SCREEN.out.html   to: params.outdir, mode: "link", overwrite: true
        FASTQ_SCREEN.out.png    to: params.outdir, mode: "link", overwrite: true
        FASTQ_SCREEN.out.report to: params.outdir, mode: "link", overwrite: true
        FASTQC2.out             to: params.outdir, mode: "link", overwrite: true
        TRIM_GALORE.out.reads   to: params.outdir, mode: "link", overwrite: true
        TRIM_GALORE.out.reports to: params.outdir, mode: 'link', overwrite: true
        HISAT2.out.bam          to: params.outdir, mode: "link", overwrite: true
        HISAT2.out.stats        to: params.outdir, mode: "link", overwrite: true

}


```


#### Example of a module (here the HISAT2 module)
```nextflow
nextflow.preview.dsl=2

params.hisat2_args = ''
params.verbose = false

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
hisat2_args = params.hisat2_args.replaceAll(/'/,"")
if (params.verbose){
	println ("[MODULE] HISAT2 ARGS: " + hisat2_args)
}

process HISAT2 {
	
	label 'bigMem'
	label 'multiCore'

    input:
	    tuple val(name), path(reads)

	output:
	    path "*bam",       emit: bam
		path "*stats.txt", emit: stats 

    script:
	cores = 8
	readString = ""
	hisat_options = hisat2_args

	// Options we add are
	hisat_options = hisat_options + " --no-unal --no-softclip "

	if (reads instanceof List) {
		readString = "-1 "+reads[0]+" -2 "+reads[1]
		hisat_options = hisat_options + " --no-mixed --no-discordant"
	}
	else {
		readString = "-U "+reads
	}
	index = params.genome["hisat2"]
	
	splices = " --known-splicesite-infile " + params.genome["hisat2_splices"]
	hisat_name = name + "_" + params.genome["name"]

	"""
	module load hisat2
	module load samtools
	hisat2 -p ${cores} ${hisat_options} -x ${index} ${splices} ${readString}  2>${hisat_name}_hisat2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${hisat_name}_hisat2.bam
	"""

}
```
