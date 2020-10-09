# Babraham Nextflow Pipelines

## Introduction

This repository contains the nextflow pipelines which are used for the bulk processing of our NGS data.  These pipelines are designed to work within the environment of our cluster, but could easily be made to work elsewhere.

These pipelines make it easy to run a standard processing of large numbers of sequences.  They generally come in one of two types:

1. A pipeline for the complete processing of a particular type of data (eg ```nf_rnaseq```) which will run several linked programs

2. A pipeline to run an individual program (eg ```nf_fastq_screen```) which can be useful when more limited processing is required

All of the pipelines are designed to be directly executable and many of them will just need to be given a set of fastq files to process.  Any that perform a mapping step will also need to supply details of the genome to process. More detailed technical information may be found at this [Nextflow @Babraham User Guide](./Docs/Nextflow_at_Babraham.md).


## Environment

The pipelines shown here are designed to operate within an environment where all of the required analysis tools and genomes are already available and where there is a queueing system for organising large numbers of jobs (Slurm is the default in our config). This means that the pipelines can run quickly and with minimal setup.  If you are working on a system which has no infrastructure and you're not interested in setting anything up then you would be better looking at pipelines such as those from [nf-core](https://github.com/nf-core) which offer a more complete infrastructure solution.

### Software
This repository provides our pipelines, but it doesn't include a copy of nextflow itself.  You will need to install this and add it to your PATH before trying to do anything with the pipelines. Please see here on how to [install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).

These pipelines assume that the software required to run the analysis components of these modules are already installed and are either already on the system ```PATH```, or are installed within [Environment modules](https://github.com/cea-hpc/modules) which the pipelines will attempt to load at runtime.

The pipelines will try to submit jobs to a local instance of [Slurm](https://slurm.schedmd.com/).  If you want to run them on the local machine, or if you have a different scheduler then you will need to change the ```executor``` line in the top level ```nextflow.config``` file.

### Genomes
Any pipeline which maps data to a reference genome will require the genome to have been downloaded and indexed prior to running the pipeline.  Genome information is supplied to the pipelines in a series of config files found in the ```genomes.d``` folder.  Each of these is a series of key value pairs which tell the pipeline some basic information about the genome (its name and the species it comes from) and the location of the indices or other configuraton files for different programs.  For example our latest mouse file contains:

```
name	GRCm38
species	Mouse
fasta	/bi/scratch/Genomes/Mouse/GRCm38/
bismark	/bi/scratch/Genomes/Mouse/GRCm38/
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
