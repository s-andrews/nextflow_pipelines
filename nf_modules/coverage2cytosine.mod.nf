nextflow.preview.dsl=2

params.singlecell = false
params.rrbs       = false
params.verbose    = false
params.pbat       = false
params.nome       = false

process COVERAGE2CYTOSINE {
	label 'bigMem'
	label 'multiCore'
		
    input:
	    path(bam)
		val (outputdir)
		val (coverage2cytosine_args)
		val (verbose)

	output:
	    path "C*",          emit: context_files
		path "*report.txt", emit: report
		path "*M-bias.txt", emit: mbias
		path "*cov.gz",     emit: coverage
	
	publishDir "$outputdir",
		mode: "link", overwrite: true
    
	script:
		
		if (verbose){
			println ("[MODULE] BISMARK COVERAGE2CYTOSINE ARGS: " + coverage2cytosine_args)
		}

		// Options we add are
		cov2cyt_options = coverage2cytosine_args + " --gzip "
		
		if (params.nome){
			println ("FLAG NOMe-Seq SPECIFIED: PROCESSING ACCORDINGLY")
		}

		
		println ("Now running command: coverage2cytosine ${cov2cyt_options} ${bam}")
		"""
		module load bismark
		bismark_methylation_extractor --bedGraph --buffer 10G -parallel ${cores} ${methXtract_options} ${bam}
		"""
		// for i in *cov.gz; do ssub -o log.$i --mem=20G --email coverage2cytosine --genome /bi/scratch/Genomes/Mouse/GRCm38/ --gzip --nome --output $i.NOME $i; done
		
}