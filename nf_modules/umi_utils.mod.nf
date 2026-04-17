nextflow.enable.dsl=2
params.no_output = false

process WHITESPACE_TO_UNDERSCORE {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

    input:
		tuple val(name), path(reads)

	output:
        tuple val(name), path ("*fastq.gz"), emit: reads


    script:

        out_nameR1 = name + "_str_R1.fastq.gz"
		in_R1 = reads[0]

		out_nameR2 = name + "_str_R2.fastq.gz"
		in_R2 = reads[1] 
			
		"""
		zcat ${in_R1} | sed 's/\\s\\+/_/g' - | gzip > ${out_nameR1}
		zcat ${in_R2} | sed 's/\\s\\+/_/g' - | gzip > ${out_nameR2}
		"""	
}


process EMBED_UMI_PE {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

    input:
		tuple val(name), path(reads)

	output:
        tuple val(name), path ("*fastq.gz"), emit: reads


    script:

        out_nameR1 = name + "_str_R1.fastq.gz"
		in_R1 = reads[0]

		out_nameR2 = name + "_str_R2.fastq.gz"
		in_R2 = reads[1] 

		"""
		zcat ${in_R1} | sed -E 's/^([^ ]+) ([^ ]+) ([^:]+):([^ ]+)/\\1_\\4 \\2 \\3/' | gzip > ${out_nameR1}
		zcat ${in_R2} | sed -E 's/^([^ ]+) ([^ ]+) ([^:]+):([^ ]+)/\\1_\\4 \\2 \\3/' | gzip > ${out_nameR2}
    	"""					
}

// The sed version was working but it is very slow. Will try again with awk. 
// doesn't make much difference
// 		zcat ${in_R1} | sed -E 's/^([^ ]+) ([^ ]+) ([^:]+):([^ ]+)/\\1_\\4 \\2 \\3/' | gzip > ${out_nameR1}
//		zcat ${in_R2} | sed -E 's/^([^ ]+) ([^ ]+) ([^:]+):([^ ]+)/\\1_\\4 \\2 \\3/' | gzip > ${out_nameR2}

// awk 'NR%4==1 && \$NF ~ /:/ {split(\$NF,a,":"); \$1=\$1"_"a[2]; sub(":"a[2]"$","",\$NF)} 1'

process EMBED_UMI_SE {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

    input:
		tuple val(name), path(reads)

	output:
        tuple val(name), path ("*fastq.gz"), emit: reads


    script:

        out_nameR1 = name + "_str_R1.fastq.gz"
		in_R1 = reads[0]
			
		"""
		zcat ${in_R1} | sed -E 's/^([^ ]+) ([^ ]+) ([^:]+):([^ ]+)/\\1_\\4 \\2 \\3/' | gzip > ${out_nameR1}
		"""	
}

//zcat ${in_R1} |  awk 'NR%4==1 {split($NF,a,":"); $1=$1"_"a[2]; sub(":"a[2]"$","",$NF)} 1' | gzip > ${out_nameR1}


process UMI_TOOLS_DEDUP {

	label 'bigMem' // 20GB

	input:
		path(bam)
		path(index) // so it's in the work folder for umi_tools to access
		val (outputdir)

	output:
        path ("*bam"), emit: bam
		path ("*log"), emit: log

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output


	script:

	bam_out = bam.baseName + "_dedup.bam" 

	"""
	module load umitools
	umi_tools dedup --log ${bam}_dedup.log --umi-separator="_" -I ${bam} > ${bam_out}
	"""

}

// I think we need to capture the output

//for i in *_sorted.bam; do INFILE=$i; OUTFILE=$(sed 's/_sorted.bam/_sorted_dedup.bam/' <<< $INFILE); ssub --mem 20G -o ${OUTFILE} -e ${OUTFILE}.err umi_tools dedup --log ${INFILE}_dedup.log --umi-separator=":" -I $INFILE; done