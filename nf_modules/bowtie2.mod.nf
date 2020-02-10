nextflow.preview.dsl=2


process BOWTIE2 {
	label 'bigMem'
	label 'multiCore'
		
    input:
	    tuple val(name), path(reads)


	output:
	    path "*bam",  emit: bam
		path "*stats.txt", emit: stats 

    script:

	readString = ""

	// Options we add are

	bowtie_options = ""

	if (reads instanceof List) {
		readString = "-1 "+reads[0]+" -2 "+reads[1]
	}
	else {
		readString = "-U "+reads
	}

	index = params.genome["bowtie2"]
	bowtie_name = name+"_"+params.genome["name"]



	"""
	module load bowtie2
	module load samtools
	bowtie2 -x ${index} -p 2 ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
	"""

}