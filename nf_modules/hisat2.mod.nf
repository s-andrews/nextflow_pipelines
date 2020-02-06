nextflow.preview.dsl=2


process HISAT2 {	
    input:
	    tuple val(name), path(reads)


	output:
	    path "*bam",  emit: bam
		path "*stats.txt", emit: stats 

    script:

	readString = ""

	if (reads instanceof List) {
		readString = "-1 reads[0] -2 reads[1]"
	}
	else {
		readString = "-U reads"
	}

	index = params.genome["hisat2"]
	splice = params.genome["hisat2_splices"]
	hisat_name = name+"_"+params.genome["name"]


	"""
	module load hisat2
	module load samtools
	hisat2 -x ${index} -p 2 ${readString}  2>${hisat_name}_stats.txt | samtools view -bS -> ${hisat_name}_hisat2.bam
	"""

}