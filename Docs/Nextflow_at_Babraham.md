![./Images/nf_at_babraham.png](./Images/nf_at_babraham.png)


# Nextflow Pipelines at the Babraham Institute - a User Guide

#### Table of Contents
- [Multi-step pipelines](#pipelines)
- [Single program pipelines](#single-program-pipelines)
- [The Nextflow Config file](#the-nextflow-config-file)
- [Nextflow Dos and Don'ts](#nextflow-dos-and-donts)
  * [Logging styles]#logging-styles
  * [-bg - Executing jobs in the background]#executing-jobs-in-the-background
- [RNA-seq workflow in more detail](#RNA-seq-worklow-in-more-detail)
  * [Example Workflow](#example-workflow)
  * [Example Module](#example-module)


We have recently transitioned from our previous pipelining system [Clusterflow](https://clusterflow.io/) to a new one based on [Nextflow](https://www.nextflow.io/docs/latest/index.html). We offer some preconfigured pipelines that generally discriminate between two different modes of operation: 

- single program pipelines (formerly known as modules)
- data type specific, multi-step pipelines

These pipelines are curated by the Babraham Bioinformatics Group, but you are of course welcome to write and use your own additional pipelines. If you need help getting started with Nextflow, please come and see any member of the Bioinformatics group who shall be happy to help.

## Pipelines:

Pipelines are supposed to work in a stream-lined and reproducible way every time they are run, and are designed so that users don't have to worry about specifying any of the plethora of options each tool provides. To this end, we try to run the individual programs of a pipeline with a pre-configured set of parameters that are (we find) sensible for the specified data type.

### List of current pipelines:


- [RNA-seq workflow](#nf_rnaseq)
- [ChIP-seq/ATAC-seq workflow](#nf_chipseq)
- [QC workflow](#nf_qc)
- [Bisulfite-seq: WGBS workflow](#nf_bisulfite_WGBS)
- [Bisulfite-seq: PBAT workflow](#nf_bisulfite_PBAT)
- [Bisulfite-seq: RRBS workflow](#nf_bisulfite_RRBS)
- [Bisulfite-seq: single-cell BS-seq workflow](#nf_bisulfite_scBSseq)
- [Bisulfite-seq: single-cell NMT-seq workflow](#nf_bisulfite_scNMT)

#### nf_rnaseq

Here is an illustration of the RNA-seq workflow:

<img src="./Images/rnaseq_pipeline.png" width="600">

    FastQC
    FastQ Screen
    Trim Galore
    Trimmed FastQC
    HISAT2
    MultiQC
    
#### nf_chipseq
    FastQC
    FastQ Screen
    Trim Galore
    Trimmed FastQC
    Bowtie2
    MultiQC
    
#### nf_qc
    FastQC
    FastQ Screen
    MultiQC
    
#### nf_bisulfite_WGBS
    FastQC
    FastQ Screen [--bisulfite]
    Trim Galore
    Trimmed FastQC
    Bismark
    Deduplicate Bismark
    Methylation extract (coverage file) [--ignore_r2 2 for PE files]
    MultiQC
    
#### nf_bisulfite_scBSseq
    FastQC
    FastQ Screen [--bisulfite]
    Trim Galore [--clip_r1 6]
    trimmed FastQC
    Bismark [--non_directional]
    deduplicate Bismark
    Methylation extract (coverage file)
    bismark2summary
    MultiQC
    
#### nf_bisulfite_RRBS
    FastQC
    FastQ Screen [--bisulfite]
    Trim Galore [--rrbs]
    Trimmed FastQC
    Bismark
    Methylation extract (coverage file)
    bismark2summary
    MultiQC
    
#### nf_bisulfite_PBAT
    FastQC
    FastQ Screen [--bisulfite]
    Trim Galore [--clip_r1 9]
    Trimmed FastQC
    Bismark [--pbat]
    Deduplicate Bismark
    Methylation extract (coverage file)
    bismark2summary
    MultiQC

## Single Program Pipelines:

#### List of current single program pipelines:
- [FastQC module](#nf_fastqc)
- [FastQ Screen module](#nf_fastq_screen)
- [Trim Galore module](#nf_trim_galore)
- [Bowtie2 module](#nf_bowtie2)
- [HISAT2 module](#nf_hisat2)
- [Bismark module](#nf_bismark)

In addition to the pre-configured default parameters, each pipeline accepts a single tool-specific additional argument. For the purpose of constructing this extra agrument, all software tools are `lowercase only` (e.g. `fastqc`, not `FastQC`), followed by `_args`, followed by one or more additional options you would like to supply:

```
--toolname_args="--additional_option value --extra_flag etc"
```

So as an example, you could run specific trimming in Trim Galore like so:

```
--trim_galore_args="--clip_r1 10 --clip_r2 10 --nextera"
```

The `--toolname_args="..."` argument should enable experienced users to customise most tools to work in more specialised ways. It should however be stressed that it should be perfectly fine to run pre-configured pipelines such as `nf_chipseq` with no need to alter any parameters manually.




## Nextflow Dos and Don'ts


### Single-hyphen options are Nextflow options

I am sure there is a ton of interesting or useful Nextflow options which we are not aware of, but we will try to list a few of these interesting concepts here. 

#### Logging styles

By default, Nextflow launches an interactive pipeline that keeps overwriting status messages in place whenever a job status updates. This looks very tidy, and it is mesmerising to watch a job progress: 

<img src="./Images/interactive_log.png" width="500">

Sometimes, especially during development of a new pipeline, this neat logging mode may 'swallow' some useful print or debugging statements. Thus, for testing purposes you can use the option `-ansi-log false`. This will allow log messages to be shown, and use a single line for each process:

<img src="./Images/ansilog_false.png" width="500">

In both modes, the Nextflow process is running interactively, and presssing `ctrl + C`, or closing the terminal window (or close the laptop) will cause the entire Nextflow pipeline process to fail and abort.

#### Executing jobs in the background 

This option sends the entire workflow into the background, thus disconnecting it from the terminal session (similar to the `nohup` command). This option launches a daemon process (which will keep running on the headnode) that watches over your workflow, and submits new jobs to the SLURM queue as required. Use this option for big pipeline jobs, or whenever you do not want to watch the status progress yourself. Upon completion, the pipeline will send you an email with the job details. This option is HIGHLY RECOMMENDED!


#### arguments maye be swallowed!



- mention: Dynamic retries with back-off

#### The Nextflow log

Sometimes it is very informative to use `nextflow log` in work directory where you tried to execute one or more jobs. This brings up previously executed jobs in this folder, along with useful stats (e.g. whether the job suceeded or errored).

The nextflow log command lists the executions run in the current folder, here is an example:
```
TIMESTAMP          	DURATION  	RUN NAME         	STATUS	REVISION ID	SESSION ID                          	COMMAND  
2020-10-08 12:33:38	53m       	lonely_bhabha  	ERR   	8a59348cdc 	021addb3-61dc-47e2-b795-64a6a30945b3	nextflow nf_chipseq --genome GRCh38 *.fastq.gz        
2020-10-08 13:54:04	1h 11m 25s	maniac_laplace 	OK    	8a59348cdc 	021addb3-61dc-47e2-b795-64a6a30945b3	nextflow nf_chipseq --genome GRCh38 *.fastq.gz -resume
```

If you want to dig in deeper yourself, you can look at the hidden file `.nextflow.log` yourself (probably for debugging only).

- `-resume` (caching)

If a pipeline workflow has been interrupted or stopped (e.g. by accidentally closing a laptop), this option will attempt to resume the workflow at the point it got interrupted by using Nextflow's caching mechanism. This may save a lot of time.
				  
https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html
				  
<img src="./Images/caching_log.png" width="800">

- mention: `--list_genomes`
- mention: `fail strategy` (retry

#### A note on options on Nextflow:

Options in Nextflow have to be supplied **exactly** as they are expected: non-matching options are simply ignored! This means that there is no auto-completion, and typos/omissions/case errors will result in the option not getting used at all. So please take extra care when supplying additional options. As an example:

```
--fastQC_args "'--nogroup'"
--fastq_args "'--nogroup'"
--fastqc "'--nogroup'"
```

would all result in the same behavior: nothing.


## RNA-seq worklow in more detail

Our implementation of Nextflow pipelines implements the new (and experimental) module system introduced with [DSL2](https://www.nextflow.io/docs/latest/dsl2.html). In essence, processes for different software tools and/or required processing steps are defined as separate **modules**. These modules are then invoked in separate **workflows** that define the different steps that are carried out for a given set of input files. Here is an example of the current RNA-seq workflow, which does the following consecutive steps for each FastQ file (single-end), or file  pair (paired-end):

- run FastQC or raw FastQ(s)
- run FastQ Screen species screen
- run Trim Galore to remove adapters and low quality base calls
- run FastQC again, this time on the adapter-/quality trimmed files
- take the trimmed FastQ files and align them to a genome using HISAT2
- Once everything is complete - run MulitQC on all files of all samples

All output will be written to the working directory.

#### Example workflow

Here is an example of the current RNA-seq workflow:

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


#### Example module

Here is an example of the current HISAT2 module:

```nextflow
nextflow.enable.dsl=2

process HISAT2 {
	
    tag "$name" // Adds name to job submission instead of (1), (2) etc.

    label 'bigMem'
    label 'multiCore'

    input:
        tuple val(name), path(reads)
	val (outputdir)
	val (hisat2_args)
	val (verbose)

    output:
	path "*bam",       emit: bam
	path "*stats.txt", emit: stats 

	publishDir "$outputdir",
	mode: "link", overwrite: true

    script:
	
	if (verbose){
	    println ("[MODULE] HISAT2 ARGS: " + hisat2_args)
	}
	
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
	hisat2 -p ${cores} ${hisat_options} -x ${index} ${splices} ${readString}  2>${hisat_name}_hisat2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> 		${hisat_name}_hisat2.bam
	"""

}
```


## Credits
This documentation was written by Felix Krueger and Simon Andrews, part of the [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk) group.
<p align="center"> <img title="Babraham Bioinformatics" id="logo_img" src="./Images/bioinformatics_logo.png"></p>
