nextflow.enable.dsl=2
params.local = ''

process BOWTIE2 {
	
	//tag "$name" // Adds name to job submission instead of (1), (2) etc.
    tag "$reads" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
	label 'multiCore'
		
    input:
	    path(reads)
		val (outputdir)
		val (bowtie2_args)
		val (verbose)

	output:
	   // tuple val(name), path ("*bam"),        emit: bam
		//path "*stats.txt", emit: stats 

	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

        println("reads: " + reads)

		cores = 8
		readString = ""

		// Options we add are
		bowtie_options = bowtie2_args
		bowtie_options +=  " --no-unal " // We don't need unaligned reads in the BAM file
		
		if (params.local == '--local'){
			// println ("Adding option: " + params.local )
			bowtie_options += " ${params.local} " 
		}

		readString = "-U " + reads
		

		index = params.genome["bowtie2"]
		bowtie_name = reads + "_" + params.genome["name"]

		"""
		//module load bowtie2
		//module load samtools
		//bowtie2 -x ${index} -p ${cores} ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
		"""

}