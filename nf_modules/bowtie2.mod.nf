nextflow.enable.dsl=2
params.local = ''
params.no_output = false
params.single_end = false

// See nf_modules/trim_galore.mod.nf for why this exists (opt-in, used by nf_houseley_run_processing).
params.expect_paired = false

process BOWTIE2 {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'
//	label 'hugeMem'
	label 'multiCore'

	// Retries with backoff, then 'ignore's rather than aborting the whole run, so one bad sample
	// (e.g. hitting the SLURM multi-file-staging issue documented in nf_modules/umi_utils.mod.nf's
	// EMBED_UMI_PE) doesn't take the rest of the group down with it.
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
	    tuple val(name), path(reads)
		val (outputdir)
		val (bowtie2_args)
		val (verbose)
		val (genome)

	output:
	    tuple val(name), path ("*bam"),        emit: bam
		path "*stats.txt", emit: stats

	script:
		if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

		def cores = 8
		def readString = ""
		def pairingCheck = ""

		// Options we add are
		def bowtie_options = bowtie2_args
		bowtie_options +=  " --no-unal " // We don't need unaligned reads in the BAM file

		if (params.local == '--local'){
			// println ("Adding option: " + params.local )
			bowtie_options += " ${params.local} "
		}

		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
			bowtie_options += " --no-discordant --no-mixed " // just output properly paired reads
		}
		else {
			readString = "-U " + reads
			// Only applies when the caller hasn't already said single-end is expected (e.g.
			// nf_traelseq defaults params.single_end = true for its R1-only TrAEL data) --
			// otherwise this would wrongly flag genuinely single-end samples as corrupted pairs.
			if (params.expect_paired && ! params.single_end) {
				// See nf_modules/trim_galore.mod.nf for why this is a real exit in the submitted
				// script rather than a Groovy-level throw -- lets it retry-then-ignore normally.
				pairingCheck = """
				echo "[PAIRING ERROR] BOWTIE2 (${name}) expected paired-end input but only received one file (${reads}). This looks like the known Nextflow/SLURM issue where one read of a pair goes missing during task preparation -- not a genuine single-end sample." >&2
				exit 1
				""".stripIndent()
			}
		}


		def index = genome["bowtie2"]
		def bowtie_name = name + "_" + genome["name"]

		"""
		${pairingCheck}
		module load bowtie2
		module load samtools
		bowtie2 -x ${index} -p ${cores} ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
		"""

}
