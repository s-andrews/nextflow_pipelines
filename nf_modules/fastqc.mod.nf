nextflow.preview.dsl=2

process FASTQC {	
	input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastqc_args)
		val (verbose)

	output:
	    tuple val(name), path ("*fastqc*"), emit: all
		path "*.zip",  emit: report
	
	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:

		if (verbose){
			println ("[MODULE] FASTQC ARGS: "+ fastqc_args)
		}

		"""
		module load fastqc
		fastqc $fastqc_args -q -t 2 ${reads}
		"""

}