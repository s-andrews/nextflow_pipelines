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
 

		//echo "nothing here" > guava.log 
		//SNPsplit_genome_preparation --vcf_file ${vcf_file_path} --reference_genome /bi/scratch/Genomes/Mouse/GRCm39/chromosomes/ ${strain_args}			

		"""
		module load snpsplit
		SNPsplit_genome_preparation --vcf_file ${vcf_file_path} --reference_genome ${genome_path} ${strain_args}		
		"""
}

