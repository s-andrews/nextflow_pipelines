nextflow.preview.dsl=2

params.hisat2_args = ''

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
hisat2_args = params.hisat2_args.replaceAll(/'/,"")
// println ("[HISAT2 MODULE, replaced] ARGS ARE: " + hisat2_args)

process HISAT2 {
	
	label 'bigMem'
	label 'multiCore'

    input:
	    tuple val(name), path(reads)

	output:
	    path "*bam",       emit: bam
		path "*stats.txt", emit: stats 

    script:
	cores = 8
	readString = ""
	hisat_options = hisat2_args

	// Options we add are
	hisat_options = hisat_options + " --no-softclip "

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
	hisat2 -p ${cores} ${hisat_options} -x ${index} ${splices} ${readString}  2>${hisat_name}_hisat2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${hisat_name}_hisat2.bam
	"""

}