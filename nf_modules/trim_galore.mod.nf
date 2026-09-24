nextflow.enable.dsl=2

params.single_end = false
params.no_output = false

// Opt-in only (used by nf_houseley_run_processing): when the caller already knows a sample is
// paired-end, this fails loudly instead of silently aligning as single-end if 'reads' somehow
// only resolves to one file at this task's own staging time -- a known Nextflow/SLURM issue
// where a real pair intermittently loses one read during staging, not a genuine SE sample.
params.expect_paired = false


process TRIM_GALORE {

	tag "$name"                         // Adds name to job submission instead of (1), (2) etc.

	label 'multiCore'                    // sets cpus = 8

	// dynamic directive
	memory { 10.GB * task.attempt }
	// Retries with backoff, then 'ignore's rather than aborting the whole run -- see
	// params.expect_paired above for why.
	errorStrategy {
		if (task.attempt <= 2) {
			sleep(Math.pow(2, task.attempt) * 30 as long)
			return 'retry'
		}
		return 'ignore'
	}
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
		def pairingCheck = ""
		if (params.single_end){
			// paired-end mode may be overridden, see e.g. TrAEL-seq Indexing
		}
		else{
			if (reads instanceof List) {
				pairedString = "--paired"
			}
			else if (params.expect_paired) {
				// Fails the actual submitted task (not a Groovy-level throw here), so this goes
				// through the normal retry-then-ignore errorStrategy above rather than aborting
				// the whole run outright.
				pairingCheck = """
				echo "[PAIRING ERROR] TRIM_GALORE (${name}) expected paired-end input but only received one file (${reads}). This looks like the known Nextflow/SLURM issue where one read of a pair goes missing during task preparation -- not a genuine single-end sample." >&2
				exit 1
				""".stripIndent()
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
		${pairingCheck}
		module load trim_galore
		module load fastqc
		trim_galore $trim_galore_opts ${pairedString} ${reads}
		"""

}
