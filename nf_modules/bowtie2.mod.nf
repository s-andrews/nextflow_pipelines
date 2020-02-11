nextflow.preview.dsl=2

params.bowtie2_args = ''

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
bowtie2_args = params.bowtie2_args.replaceAll(/'/,"")
println ("[BOWTIE2 MODULE, replaced] ARGS ARE: " + bowtie2_args)

process BOWTIE2 {
	label 'bigMem'
	label 'multiCore'
		
    input:
	    tuple val(name), path(reads)

	output:
	    path "*bam",  emit: bam
		path "*stats.txt", emit: stats 

    script:
	cores = 8
	readString = ""

	// Options we add are
	bowtie_options = bowtie2_args
	bowtie_options = bowtie_options + " --no-unal "
	
	if (reads instanceof List) {
		readString = "-1 "+reads[0]+" -2 "+reads[1]
	}
	else {
		readString = "-U "+reads
	}

	index = params.genome["bowtie2"]
	bowtie_name = name + "_" + params.genome["name"]

	"""
	module load bowtie2
	module load samtools
	bowtie2 -x ${index} -p ${cores} ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
	"""

}