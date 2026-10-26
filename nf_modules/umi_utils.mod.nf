nextflow.enable.dsl=2
params.no_output = false

process WHITESPACE_TO_UNDERSCORE {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
	maxRetries 2

    input:
		tuple val(name), path(reads)

	output:
        tuple val(name), path ("*fastq.gz"), emit: reads


    script:

        def out_nameR1 = name + "_str_R1.fastq.gz"
		def in_R1 = reads[0]

		def out_nameR2 = name + "_str_R2.fastq.gz"
		def in_R2 = reads[1]

		"""
		set -o pipefail
		zcat ${in_R1} | sed 's/\\s\\+/_/g' - | gzip > ${out_nameR1}
		zcat ${in_R2} | sed 's/\\s\\+/_/g' - | gzip > ${out_nameR2}
		"""
}


process EMBED_UMI_PE {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem' // 20GB

	// zcat occasionally hits a transient read failure on a large file that's just landed on a
	// network filesystem (e.g. NFS attribute-cache staleness across SLURM nodes). Without
	// pipefail, only gzip's exit code was checked, so a failed zcat still produced a "successful"
	// but empty/truncated output file -- silently corrupting data instead of failing loudly. Now
	// it fails immediately and retries with backoff, the same pattern already used by TRIM_GALORE.
	//
	// On large real groups, one specific sample has repeatedly (and reproducibly, across many
	// separate runs) had only ONE of its two input files actually staged into the task's work
	// directory by Nextflow itself -- confirmed by inspecting the work dir directly, not just the
	// generated script. This has been ruled out as a pairing bug in our own code (exhaustively
	// verified correct in isolation, at full 30-sample scale, against the real data) and ruled out
	// as a submission-burst issue (still happens with maxForks capping concurrency). It only ever
	// happens under real SLURM dispatch, pointing at something in Nextflow's own grid-executor
	// staging, which is out of our control to fix directly. So retries still get 2 attempts (in
	// case it is transient for some samples), but once exhausted this now 'ignore's rather than
	// aborting the whole run -- one bad sample no longer takes the other ~30 down with it. The
	// failure is still loud and visible in the log either way; the affected sample can be re-run
	// on its own afterwards.
	errorStrategy {
		if (task.attempt <= 2) {
			sleep(Math.pow(2, task.attempt) * 30 as long)
			return 'retry'
		}
		return 'ignore'
	}
	maxRetries 2
	maxForks 5

    input:
		tuple val(name), path(reads)

	output:
        tuple val(name), path ("*fastq.gz"), emit: reads


    script:

        def out_nameR1 = name + "_str_R1.fastq.gz"
		def in_R1 = reads[0]

		def out_nameR2 = name + "_str_R2.fastq.gz"
		def in_R2 = reads[1]

		// EMBED_UMI_PE is inherently paired-only (see EMBED_UMI_SE for single-end), so unlike the
		// opt-in params.expect_paired check elsewhere, a missing second file is always wrong here.
		// Fails with a clear message instead of the cryptic 'zcat null'/'gzip: null.gz' pipefail
		// would otherwise produce.
		def pairingCheck = (in_R2 == null) ? """
			echo "[PAIRING ERROR] EMBED_UMI_PE (${name}) expected 2 input files (R1 and R2) but only received 1 (${in_R1}). This looks like the known Nextflow/SLURM issue where one read of a pair goes missing during task preparation -- not a genuine problem with this sample." >&2
			exit 1
			""".stripIndent() : ""

		"""
		${pairingCheck}
		set -o pipefail
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

	// See EMBED_UMI_PE above for why this is 'ignore' rather than the retry-only default, and why
	// maxForks is capped.
	errorStrategy {
		if (task.attempt <= 2) {
			sleep(Math.pow(2, task.attempt) * 30 as long)
			return 'retry'
		}
		return 'ignore'
	}
	maxRetries 2
	maxForks 5

    input:
		tuple val(name), path(reads)

	output:
        tuple val(name), path ("*fastq.gz"), emit: reads


    script:

        def out_nameR1 = name + "_str_R1.fastq.gz"
		def in_R1 = reads[0]

		"""
		set -o pipefail
		zcat ${in_R1} | sed -E 's/^([^ ]+) ([^ ]+) ([^:]+):([^ ]+)/\\1_\\4 \\2 \\3/' | gzip > ${out_nameR1}
		"""
}

//zcat ${in_R1} |  awk 'NR%4==1 {split($NF,a,":"); $1=$1"_"a[2]; sub(":"a[2]"$","",$NF)} 1' | gzip > ${out_nameR1}


process UMI_TOOLS_DEDUP {

	label 'bigMem' // 20GB

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
		path(bam)
		path(index) // so it's in the work folder for umi_tools to access
		val (outputdir)
		val (paired_end)
		val (umi_separator)

	output:
        path ("*bam"), emit: bam
		path ("*log"), emit: log


	script:

	def options = ""
	if (paired_end) {
		options += "--paired"
	}

	def bam_out = bam.baseName + "_dedup.bam"

	"""
	module load umitools
	umi_tools dedup --log ${bam}_dedup.log ${options} --umi-separator="${umi_separator}" -I ${bam} > ${bam_out}
	"""

}

// I think we need to capture the output

//for i in *_sorted.bam; do INFILE=$i; OUTFILE=$(sed 's/_sorted.bam/_sorted_dedup.bam/' <<< $INFILE); ssub --mem 20G -o ${OUTFILE} -e ${OUTFILE}.err umi_tools dedup --log ${INFILE}_dedup.log --umi-separator=":" -I $INFILE; done
