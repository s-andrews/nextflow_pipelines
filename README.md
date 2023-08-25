# Babraham Nextflow Pipelines

## Introduction

This repository contains the nextflow pipelines which are used for the bulk processing of our NGS data.  These pipelines are designed to work within the environment of our cluster, but could easily be made to work elsewhere.

These pipelines make it easy to run a standard processing of large numbers of sequences.  They generally come in one of two types:

1. A pipeline for the complete processing of a particular type of data (eg ```nf_rnaseq```) which will run several linked programs

2. A pipeline to run an individual program (eg ```nf_fastq_screen```) which can be useful when more limited processing is required

All of the pipelines are designed to be directly executable and many of them will just need to be given a set of fastq files to process.  Any that perform a mapping step will also need to supply details of the genome to process. More detailed technical information may be found at this [Nextflow @Babraham User Guide](./Docs/Nextflow_at_Babraham.md).


## Environment
The pipelines here are designed to work on a linux server which uses [environment modules](https://modules.readthedocs.io/en/latest/) to dynamically load software packages.  We do not include the analysis tools used by the pipelines and you must install these separately. Each of the pipelines will issue appropriate ```module load``` commands for the software it requires and will fail if these are not present.

You will need to have an installation of [nextflow](https://www.nextflow.io/) to run the pipelines, and this will need to be configured with an appropriate [executor](https://www.nextflow.io/docs/latest/executor.html).  Nextflow can talk to a variety of different queueing systems if running on a cluster, or can run in a local mode on a stand alone server.  Our default setup uses the [slurm](https://slurm.schedmd.com/) queueing system, but you can change this by editing the ```nextflow.config``` file we distribute.

If you are working on a system which has no infrastructure and you're not interested in setting anything up then you may prefer looking at pipelines such as those from [nf-core](https://github.com/nf-core) which offer a more complete infrastructure solution.


### Genomes
Any pipeline which maps data to a reference genome will require the genome to have been downloaded and indexed prior to running the pipeline.  Genome information is supplied to the pipelines in a series of config files found in the ```genomes.d``` folder.  Each of these is a series of key value pairs which tell the pipeline some basic information about the genome (its name and the species it comes from) and the location of the indices or other configuraton files for different programs.  For example our latest mouse file contains:

```
name	GRCm38
species	Mouse
fasta	/bi/scratch/Genomes/Mouse/GRCm38/
bismark	/bi/scratch/Genomes/Mouse/GRCm38/
cellranger  /bi/apps/cellranger/references/refdata-gex-mm10-2020-A/
bowtie	/bi/scratch/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38
bowtie2	/bi/scratch/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38
hisat2	/bi/scratch/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38
gtf	/bi/scratch/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38.90.gtf
hisat2_splices	/bi/scratch/Genomes/Mouse/GRCm38/Mus_musculus.GRCm38.90.hisat2_splices.txt
```

This means that within the pipelines you can just specify ```--genome GRCm38``` and all of the required information will then be extracted from this config file.

If you install these pipelines on your own system then you will need to prepare the genomes which you want to use and populate the corresponding files in ```genomes.d```

More detailed information please refer to the more detailed [Nextflow @Babraham User Guide](./Docs/Nextflow_at_Babraham.md).

# Problems and bugs
If you have any problems using these pipelines then feel free to file issues in our [Github repository](https://github.com/s-andrews/nextflow_pipelines/issues).  For other questions you can email [babraham.bioinformatics@babraham.ac.uk](mailto://babraham.bioinformatics@babraham.ac.uk)
