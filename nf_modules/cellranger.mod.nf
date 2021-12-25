nextflow.enable.dsl=2

// parameters passed in by specialised pipelines

process CELLRANGER_COUNT {
	
	tag "$sampleID" // Adds name to job submission instead of (1), (2) etc.
	label 'cellranger'

	cpus 16
  	memory { 35.GB * task.attempt }  
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
  	maxRetries 2
	
	// label 'mem40G'
	// label 'multiCore'
		
    input:
		val (sampleID)
		val (fastqs_dir)
		val (outputdir)
		val (cellranger_args)
		val (verbose)

	output:
	    path "${sampleID}/outs/*"

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		cores = 16
		readString = ""

		// These parameters are allowed for CellRanger count
		// Argument	Description
		
		// --id	A unique run ID string: e.g. sample345
		
		// --sample	Sample name as specified in the sample sheet supplied to cellranger mkfastq.
		// Can take multiple comma-separated values, which is helpful if the same library was sequenced on multiple flowcells and the sample name used (and therefore fastq file prefix) is not identical between them.
		// Doing this will treat all reads from the library, across flowcells, as one sample.
		// If you have multiple libraries for the sample, you will need to run cellranger count on them individually, and then combine them with cellranger aggr.
		// Allowable characters in sample names are letters, numbers, hyphens, and underscores.
		
		// --transcriptome	Path to the Cell Ranger compatible transcriptome reference e.g.
		// For a human-only sample, use /opt/refdata-gex-GRCh38-2020-A
		// For a human and mouse mixture sample, use /opt/refdata-gex-GRCh38-and-mm10-2020-A

		command = 'cellranger count '
		if (verbose){
			println ("[MODULE] Constructing the command:\n[MODULE]: >>$command<<")
		}
		command += " $cellranger_args "
		if (verbose){
			println ("[MODULE] CELLRANGER ARGS: " + cellranger_args)
		}

		command += " --localcores=$cores --localmem=${task.memory.toGiga()} --transcriptome=" + params.genome["cellranger"]
		if (verbose){
			println ("[MODULE]: >>$command<<")
		}
		
		command += " --sample=$sampleID"
		if (verbose){
			println ("[MODULE]: >>$command<<")
		}

		command += " --id=$sampleID"
		if (verbose){
			println ("[MODULE]: >>$command<<")
		}
		
		command += " --fastqs=$fastqs_dir"
		if (verbose){
			println ("[MODULE]: >>$command<<")
		}
		// Sample command:
		// ssub -o cellranger.log --email --j RANGER9 -c 16 -m 40G  cellranger count --transcriptome=/bi/apps/cellranger/references/refdata-gex-mm10-2020-A/ --fastqs=/bi/scratch/Felix/For_Carine/210315_multiome/MULTIseq/ --localcores=16 --localmem=32 --sample=SITTA12_gastr_d3_MULTIseq_3prime_RNA_A --id=SITTA12_gastr_d3_MULTIseq_3prime_RNA_A
		

		"""
		module load cellranger
		module load ssub

		$command
        
		"""
		
}