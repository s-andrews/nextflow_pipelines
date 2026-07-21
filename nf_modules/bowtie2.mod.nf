nextflow.enable.dsl=2
params.local = ''
params.no_output = false

process BOWTIE2 {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
//	label 'hugeMem'
	label 'multiCore'

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bowtie2_args)
		val (verbose)
		val (genome)

	output:
	    tuple val(name), path ("*bam"),        emit: bam
		path "*stats.txt", emit: stats

	script:
		if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

		def cores = 8
		def readString = ""

		// Options we add are
		def bowtie_options = bowtie2_args
		bowtie_options +=  " --no-unal " // We don't need unaligned reads in the BAM file

		if (params.local == '--local'){
			// println ("Adding option: " + params.local )
			bowtie_options += " ${params.local} "
		}

		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
			bowtie_options += " --no-discordant --no-mixed " // just output properly paired reads
		}
		else {
			readString = "-U " + reads
		}


		def index = genome["bowtie2"]
		def bowtie_name = name + "_" + genome["name"]

		"""
		module load bowtie2
		module load samtools
		bowtie2 -x ${index} -p ${cores} ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
		"""

}
