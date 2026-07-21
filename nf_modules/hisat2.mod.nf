nextflow.enable.dsl=2

process HISAT2 {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
	label 'multiCore'

	publishDir { outputdir },
		mode: "link", overwrite: true

    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (hisat2_args)
		val (verbose)
		val (genome)

	output:
	    path "*bam",       emit: bam
		path "*stats.txt", emit: stats

    script:

		if (verbose){
			println ("[MODULE] HISAT2 ARGS: " + hisat2_args)
		}

		def cores = 8
		def readString = ""
		def hisat_options = hisat2_args

		// Options we add are
		hisat_options = hisat_options + " --no-unal --no-softclip --new-summary"

		if (reads instanceof List) {
			readString = "-1 "+reads[0]+" -2 "+reads[1]
			hisat_options = hisat_options + " --no-mixed --no-discordant"
		}
		else {
			readString = "-U "+reads
		}
		def index = genome["hisat2"]

		// TODO: need to add a check if the splice-site infile is present or not, and leave out this parameter otherwise
		def splices = " --known-splicesite-infile " + genome["hisat2_splices"]
		def hisat_name = name + "_" + genome["name"]

		"""
		module load hisat2
		module load samtools
		hisat2 -p ${cores} ${hisat_options} -x ${index} ${splices} ${readString}  2>${hisat_name}_hisat2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${hisat_name}_hisat2.bam
		"""

}
