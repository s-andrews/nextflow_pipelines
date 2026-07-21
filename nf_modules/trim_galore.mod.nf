nextflow.enable.dsl=2

params.single_end = false
params.no_output = false


process TRIM_GALORE {

	tag "$name"                         // Adds name to job submission instead of (1), (2) etc.

	label 'multiCore'                    // sets cpus = 8

	// dynamic directive
	memory { 10.GB * task.attempt }
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
	maxRetries 2

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
	    tuple val (name), path (reads)
		val (outputdir)
		val (trim_galore_args)
		val (verbose)
		val (singlecell)
		val (rrbs)
		val (pbat)
		val (clock)
		val (three_prime_clip_R1)
		val (three_prime_clip_R2)

	output:
	    tuple val(name), path ("*fq.gz"), emit: reads
		path "*trimming_report.txt", optional: true, emit: report

    script:
		if (verbose){
			println ("[MODULE] TRIM GALORE ARGS: " + trim_galore_args)
		}

		def pairedString = ""
		if (params.single_end){
			// paired-end mode may be overridden, see e.g. TrAEL-seq Indexing
		}
		else{
			if (reads instanceof List) {
				pairedString = "--paired"
			}
		}

		def trim_galore_opts = trim_galore_args

		// Set multi-core
		trim_galore_opts += " -j 8 "

		// Specialised Epigenetic Clock Processing
		if (clock){
			trim_galore_opts += " --breitling "
		}
		else{
			if (singlecell){
				trim_galore_opts += " --clip_r1 6 "
				if (pairedString == "--paired"){
					trim_galore_opts += " --clip_r2 6 "
				}
			}

			if (rrbs){
				trim_galore_opts = trim_galore_opts + " --rrbs "
			}

			if  (pbat){
				trim_galore_opts = trim_galore_opts + " --clip_r1 $pbat "
				if (pairedString == "--paired"){
					trim_galore_opts = trim_galore_opts + " --clip_r2 $pbat "
				}
			}

			// Second step of Clock processing:
			if (three_prime_clip_R1 && three_prime_clip_R2){
				trim_galore_opts +=	" --three_prime_clip_R1 ${three_prime_clip_R1} --three_prime_clip_R2 ${three_prime_clip_R2} "
			}
		}

		"""
		module load trim_galore
		module load fastqc
		trim_galore $trim_galore_opts ${pairedString} ${reads}
		"""

}
