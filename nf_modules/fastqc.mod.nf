nextflow.enable.dsl=2
params.nogroup = false
params.no_output = false

process FASTQC {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	// paired-end files get 8 cores, single-end files get 4
	cpus { reads instanceof List ? 8 : 4 }

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (fastqc_args)
		val (verbose)

	output:
	    tuple val(name), path ("*fastqc*"), emit: all
		path "*.zip",  emit: report

	script:

		if (params.nogroup){
			// println ("ADDING --nogroup: " + fastqc_args)
			fastqc_args += " --nogroup "
		}

		if (verbose){
			println ("[MODULE] FASTQC ARGS: "+ fastqc_args)
		}

		"""
		module load fastqc
		fastqc $fastqc_args -q -t ${task.cpus} ${reads}
		"""
}
