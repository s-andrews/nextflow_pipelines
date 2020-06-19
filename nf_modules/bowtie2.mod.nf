nextflow.preview.dsl=2

process BOWTIE2 {
	
	label 'bigMem'
	label 'multiCore'
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bowtie2_args)
		val (verbose)

	output:
	    path "*bam",  emit: bam
		path "*stats.txt", emit: stats 

	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

		cores = 8
		readString = ""

		// Options we add are
		bowtie_options = bowtie2_args
		bowtie_options +=  " --no-unal " // We don't need unaligned reads in the BAM file
		
		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
			bowtie_options += " --no-discordant --no-mixed " // just output properly paired reads
		}
		else {
			readString = "-U " + reads
		}

		index = params.genome["bowtie2"]
		bowtie_name = name + "_" + params.genome["name"]

		"""
		module load bowtie2
		module load samtools
		bowtie2 -x ${index} -p ${cores} ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
		"""

}