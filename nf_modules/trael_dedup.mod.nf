nextflow.enable.dsl=2

params.no_output = false

// Used only by nf_traelseq_paired. Trypanosome paired-end TrAEL mapping should only deduplicate
// reads that map uniquely once the paired end is taken into account -- multi-mapping reads are
// dropped here rather than kept in a lowqual fallback bucket (unlike the single-end UMIBAM2,
// which keeps low-mapq reads and dedups them by UMI+sequence alone). The filtered, unique
// properly-paired reads produced here go on to UMI_TOOLS_DEDUP (nf_modules/umi_utils.mod.nf) --
// the custom umibam2_PE_unique.py fallback-bucket logic isn't needed for this path, since
// multi-mappers have already been dropped by the time dedup runs.
params.unique_mapq = 20

process FILTER_UNIQUE_PAIRS {

	tag "$bam" // Adds name to job submission instead of (1), (2) etc.
	label 'bigMem'

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
		path(bam)
		val (outputdir)
		val (mapq)
		val (verbose)

	output:
		path ("*unique.bam"), emit: bam

	script:
		if (verbose){
			println ("[MODULE] FILTER_UNIQUE_PAIRS MAPQ threshold: " + mapq)
		}

		def bam_out = bam.baseName + "_unique.bam"

		// A pair is kept only as a whole fragment, never as one mate on its own: 'samtools view
		// -q' filters individual BAM records, so on a pair where Bowtie2 scored the two mates
		// differently (common -- TrAEL's R1 is short and aggressively trimmed, and can score
		// lower/higher than R2 independently of whether the fragment as a whole is uniquely
		// placed), it silently keeps the passing mate and drops the other. That orphan then goes
		// into UMI_TOOLS_DEDUP's --paired mode as a lone, unpaired record -- with its SAM flags
		// still claiming a mapped mate that no longer exists in the file, so nothing errors, and
		// the missing mate (e.g. R2) is easy to miss until someone checks read counts.
		//
		// Keep the fragment if EITHER mate reaches the MAPQ threshold, per the whole point of
		// aligning as pairs: even a mate that looks ambiguous alone (short/repetitive) can be
		// uniquely resolved once Bowtie2's own insert-size/orientation constraint (--no-discordant
		// --no-mixed, already enforced upstream in BOWTIE2) ties it to a confident mate.
		"""
		module load samtools
		samtools view -f 2 ${bam} | awk -v m=${mapq} '\$5+0 >= m {print \$1}' | sort -u > keep_names.txt
		samtools view -b -f 2 -N keep_names.txt ${bam} -o ${bam_out}
		"""
}
