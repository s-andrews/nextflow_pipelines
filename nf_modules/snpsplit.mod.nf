nextflow.enable.dsl=2

process SNP_SPLIT_GENOME_PREP {

	tag "$strain" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem'

	input:
		val (outputdir)
		path (vcf)
		val (genome)
		val (strain)
		val (strain2)

	output:
		path ("*"),  emit: all
	
	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:

		if(vcf.exists()) {
			vcf_file = file(vcf)
			vcf_file_path = vcf_file.resolve()
    		println("Using vcf file ${vcf_file_path}")
		} else {
			println("\nCouldn't find vcf from module, exiting...\n")
			exit 1
		}

		// check this again
		genome_path = "${genome}"

		strain_args = "--strain ${strain}"

		if (strain2 == ""){
		 	println ("\nUsing one strain for genome preparation: " + strain)
			//println("strain args = " + strain_args)	 
		} else {
			strain_args += " --strain2 ${strain2}"
			println("\nUsing two strains for genome preparation: " + strain + " and " + strain2)
			//println("strain args = " + strain_args)
		}
 	

		"""
		module load snpsplit
		SNPsplit_genome_preparation --vcf_file ${vcf_file_path} --reference_genome ${genome_path} ${strain_args}		
		"""
}

process SNP_SPLIT {

	label 'bigMem'

	input:
		val (outputdir)
		path (snp_file)
		path (mapped_bam)
		val (bisulfite)
		val (hic)

	output:
		path ("*"),  emit: all
	
	publishDir "$outputdir",
		mode: "link", overwrite: true

	script:

		if(snp_file.exists()) {
			//snps = file(snp_file)
			//snp_file_path = snps.resolve()
			//println("Using snp file ${snp_file_path}")
			snp_file_path = snp_file
			println("Using snp file ${snp_file_path}")
		} else {
			println("\nCouldn't find snp_file from module, exiting...\n")
			exit 1
		}

		if(mapped_bam.exists()) {
			//bam = file(mapped_bam)
			//bam_file_path = bam.resolve()
			//println("Using bam file ${bam_file_path}")
			bam_file_path = mapped_bam
		} else {
			println("\nCouldn't find bam_file from module, exiting...\n")
			exit 1
		}

		snpsplit_args = ""

		if (bisulfite) {
			snpsplit_args += "--bisulfite"
		}
		if (hic) {
			snpsplit_args += "--hic"
		}

		"""
		module load snpsplit
		SNPsplit --snp_file ${snp_file_path} ${bam_file_path} ${snpsplit_args}		
		"""

}
