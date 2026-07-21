nextflow.enable.dsl=2

process UMIBAM {

	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	// dynamic directive to increase memory as required
	cpus 1
	memory { 20.GB * task.attempt }
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
  	maxRetries 5

	publishDir { outputdir },
		mode: "link", overwrite: true

	input:
	    tuple val(name), path(bam)
		val (outputdir)
		val (umibam_args)
		val (verbose)
		val (dual)

	output:
		path "*report.txt", emit: report
		tuple val(name), path ("*bam"),        emit: bam


    script:
		if (verbose){
			println ("[MODULE] UMIBAM ARGS: " + umibam_args)
		}

		def umibam_opts = umibam_args

		// Specialised Epigenetic Clock Processing
		if (dual){
			umibam_opts += " --double_umi "
		}
		else{
			umibam_opts += " --umi "
		}

		"""
		module load UmiBam

		UmiBam $umibam_opts  $bam

		rename UMI_d d *
		"""

		// The output files should be renamed so that they bismark2report picks up everything

		// renaming files using Bash
		// for f in *UMI_dedup* ; do mv "\$f" "\${f/UMI_/}" ; done

		// renaming using rename (works on our cluster)
		// rename UMI_d d *

		// A third option should be the saveAs directive (https://www.nextflow.io/docs/latest/process.html#publishdir)
		// unclear to me at the moment though how this would work exactly
}
