nextflow.preview.dsl=2


process HISAT2 {	
    input:
	    tuple val(name), path(reads)


	output:
	    path "*bam",  emit: bam
		path "*stats.txt", emit: stats 

    script:

	readString = ""

	// Options we add are
	// --sp 1000,1000 to stop soft clipping

	hisat_options = "--sp 1000,1000"

	if (reads instanceof List) {
		readString = "-1 "+reads[0]+" -2 "+reads[1]
		hisat_options = hisat_options+" --no-mixed --no-discordant"
	}
	else {
		readString = "-U "+reads
	}

	index = params.genome["hisat2"]
	splice = params.genome["hisat2_splices"]
	hisat_name = name+"_"+params.genome["name"]



	"""
	module load hisat2
	module load samtools
	hisat2 -x ${index} -p 2 ${hisat_options} ${readString}  2>${hisat_name}_hisat2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${hisat_name}_hisat2.bam
	"""

}