# Nextflow Pipelines at the Babraham Institute

### General comments: 

We are currently transitioning from our previous pipelining system (Clusterflow) to a new one based on [Nextflow](https://www.nextflow.io/docs/latest/index.html). We offer some preconfigured pipelines that generally discriminate between two different modes of operation: 

- data type specific pipelines
- single program pipelines (formerly known as modules)

These pipelines are curated by the Babraham Bioinformatics Group, but you are of course welcome to write and use your own additional pipelines. If you need help getting started with Nextflow, please come and see any member of the Bioinformatics group who shall be happy to help.

## Pipelines:

Pipelines are supposed to work in a stream-lines and reproducible way every time they are run. We try to run them with default settings that we deem sensible.

List of current pipelines:

- nf_qc
- nf_rnaseq
- nf_chipseq
- nf_bisulfite_WGBS
- nf_bisulfite_scBSseq
- nf_bisulfite_RRBS



## Single Progam Pipelines:
- nf_fastqc
- nf_fastq_screen
- nf_trim_galore
- nf_trim_galore_speciality
- nf_bowtie2
- nf_hisat2
- nf_bismark


All preconfigured pipelines take one additional argument, which has to be exactly in the following form to work:

```
--toolname_args "'--additional_option value --extra_flag etc.'"
```

So as an example, you could run specific trimming in Trim Galore like so:

```
--trim_galore_args "'--clip_r1 10 --clip_r2 10 --nextera'"
```
