nextflow.enable.dsl=2
params.no_output = false

process SAMTOOLS_SORT{

	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
		path(bam)
		val (outputdir)
		val (samtools_sort_args)
		val (verbose)

	output:
		// path "*report.txt", emit: report
		path "*bam",        emit: bam

    script:
		def samtools_sort_options = samtools_sort_args

		if (verbose){
			println ("[MODULE] SAMTOOLS SORT ARGS: " + samtools_sort_args)
		}

		// TODO: Find more elegant way to strip file ending of input BAM file

		"""
		module load samtools
		samtools sort $samtools_sort_options $bam -o ${bam}_sorted.bam
		rename .bam_sorted _sorted *
    	"""
}

// this is going in as a temporary measure so we can just easily pass in input that has a name and reads
process SAMTOOLS_SORT2{

	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
		tuple val(name), path(bam)
		val (outputdir)
		val (samtools_sort_args)
		val (verbose)

	output:
		// path "*report.txt", emit: report
		path "*bam",        emit: bam

    script:
		def samtools_sort_options = samtools_sort_args

		if (verbose){
			println ("[MODULE] SAMTOOLS SORT ARGS: " + samtools_sort_args)
		}

		// TODO: Find more elegant way to strip file ending of input BAM file

		"""
		module load samtools
		samtools sort $samtools_sort_options $bam -o ${bam}_sorted.bam
		rename .bam_sorted _sorted *
    	"""
}


process SAMTOOLS_INDEX{

	tag "$bam"     // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
		path(bam)
		val (outputdir)
		val (samtools_index_args)
		val (verbose)

	output:
		path "*.bai",     emit: bai

    script:
		def samtools_index_options = samtools_index_args

		if (verbose){
			println ("[MODULE] SAMTOOLS INDEX ARGS: " + samtools_index_args)
		}

		"""
		module load samtools
		samtools index $samtools_index_options $bam
		"""
}

// output is the index file and the sorted bam so that they stay together
process SAMTOOLS_INDEX2{

	tag "$bam"     // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
		path(bam)
		val (outputdir)
		val (samtools_index_args)
		val (verbose)

	output:
		path "*.bai",     emit: bai
		path (bam), 	  emit: bam

    script:
		def samtools_index_options = samtools_index_args

		if (verbose){
			println ("[MODULE] SAMTOOLS INDEX ARGS: " + samtools_index_args)
		}

		"""
		module load samtools
		samtools index $samtools_index_options $bam
		"""
}


process SAM2BAM{

	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
		tuple val(name), path(sam)
		val (outputdir)

	output:
		path "*bam",        emit: bam

    script:

		// we don't want to use the 'name' field as this probably doesn't have the genome info
		// so we use baseName to get the filename without the suffix
		def bam = sam.baseName + ".bam"

		"""
		module load samtools
		samtools view -b $sam -o ${bam}
    	"""
}

// filter on mapq > 30
process SAMTOOLS_FILT{

	tag "$bamfile" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	//publishDir { outputdir },
	//	mode: "link", overwrite: true, enabled: !params.no_output  // change this if we wnat the file to be written out rather than just in the work folder

	input:
		path(bamfile)
		val (outputdir)

	output:
		path "*bam",        emit: bam

    script:

		def bam_out = bamfile.baseName + "_mapq30.bam"

		"""
		module load samtools
		samtools view -b $bamfile -q 30 -o ${bam_out}
		"""
}
