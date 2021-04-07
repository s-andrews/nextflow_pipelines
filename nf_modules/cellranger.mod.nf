nextflow.enable.dsl=2

// parameters passed in by specialised pipelines
params.singlecell = ''
params.pbat = false
params.read_identity = ''

process CELLRANGER {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	cpus { 16 }
  	memory { 35.GB * task.attempt }  
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
  	maxRetries 2
	
	// label 'mem40G'
	// label 'multiCore'
	// label 'quadCore'
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (cellranger_args)
		val (verbose)

	output:
	    // tuple val(name), path ("*bam"),        emit: bam
		// path "*report.txt", emit: report
		// we always pass back the original name so we can use .join() later on, e.g. for bismark2bedGraph
		// tuple val(name), path ("*unmapped_reads_1.fq.gz"), optional: true, emit: unmapped1
		// tuple val(name), path ("*unmapped_reads_2.fq.gz"), optional: true, emit: unmapped2

	publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		cores = 16
		readString = ""

		if (verbose){
			println ("[MODULE] CELLRANGER ARGS: " + cellranger_args)
		}

		// Options we add are
		// bismark_options = bismark_args
		println(reads)		
		
		// nf (params.unmapped){
		// 	bismark_options += " --unmapped "
		// 	unmapped_1_name = name + "_unmapped_R1"
		// 	unmapped_2_name = name + "_unmapped_R2"
		// }

		// if (params.pbat){
		// 	bismark_options += " --pbat "
		// }

		// if (reads instanceof List) {
		// 	readString = "-1 "+reads[0]+" -2 "+reads[1]
		// }
		// else {
		// 	readString = reads
		// }

		// index = "--genome " + params.genome["bismark"]

		// unmapped_name = ''	
		// 	// add Genome build and aligner to output name
		// if (params.read_identity == "1" || params.read_identity == "2"){
		// 	// println ("FILENAME IS: $reads")
		// 	if (params.read_identity == "1"){
		// 		unmapped_name = name + "_unmapped_R1"
		// 	}
		// 	else{
		// 		unmapped_name = name + "_unmapped_R2"
		// 	}

		// 	if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
		// 		bismark_name = unmapped_name + "_" + params.genome["name"] + "_bismark_hisat2"
		// 	}
		// 	else{ // default is Bowtie 2
		// 		bismark_name = unmapped_name + "_" + params.genome["name"] + "_bismark_bt2"
		// 	}
		// }
		// else{
		// 	if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
		// 		bismark_name = name + "_" + params.genome["name"] + "_bismark_hisat2"
		// 	}
		// 	else{ // default is Bowtie 2
		// 		bismark_name = name + "_" + params.genome["name"] + "_bismark_bt2"
		// 	}
		// }
		command = ''	
		println ("[MODULE] Constructing the command:\n[MODULE]: >>$command<<")
		command += " --localcores=$cores --localmem=32 --transcriptome=" + params.genome["cellranger"]
		println ("[MODULE]: >>$command<<")
		// Sample command:
		// ssub --mem=35G --cores 16 -o cellranger.log -j Ranger2 --email cellranger count --localmem=32 --localcores 16 --transcriptome=/bi/apps/cellranger/references/refdata-gex-mm10-2020-A/ --fastqs=/bi/group/bioinf/Stephen_Clark/210205_2more_10X_RNAseq/SIGAD11/ --id=SIGAD11_TET-TKO_Day8_EB --sample=SIGAD11_TET-TKO_Day8_EB
		//  --j RANGER9 -c 16 -m 40G  cellranger count --transcriptome=/bi/apps/cellranger/references/refdata-gex-mm10-2020-A/ --fastqs=/bi/scratch/Felix/For_Carine/210315_multiome/MULTIseq/ --localcores=16 --localmem=32 --sample=SITTA12_gastr_d3_MULTIseq_3prime_RNA_A --id=SITTA12_gastr_d3_MULTIseq_3prime_RNA_A
		"""
		module load cellranger
		module load ssub

		"""
		// cellranger count -- $cores --basename $bismark_name $index $bismark_options $readString
		
}