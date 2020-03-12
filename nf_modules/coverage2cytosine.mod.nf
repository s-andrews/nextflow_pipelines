nextflow.preview.dsl=2

params.singlecell = false
params.rrbs       = false
params.verbose    = false
params.pbat       = false
params.nome       = false

genome = params.genome["bismark"]

process COVERAGE2CYTOSINE {
	label 'hugeMem'
		
    input:
	    path(coverage_file)
		val (outputdir)
		val (coverage2cytosine_args)
		val (verbose)

	output:
	    path "*{report.txt.gz,report.txt}", emit: report
		path "*{.cov.gz,.cov}",             emit: coverage
	
	publishDir "$outputdir",
		mode: "link", overwrite: true
    
	script:
		
		// removing the file extension from the input file name 
		// (https://www.nextflow.io/docs/latest/script.html#removing-part-of-a-string)
		outfile_basename = coverage_file.toString()  // Important to convert nextflow.processor.TaskPath object to String first
		outfile_basename = (outfile_basename - ~/.bismark.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov$/)

		if (verbose){
			println ("[MODULE] BISMARK COVERAGE2CYTOSINE ARGS: " + coverage2cytosine_args)
			println ("Bismark Genome is: " + genome)
		}

		// Options we add are
		cov2cyt_options = coverage2cytosine_args + " --gzip "
		
		if (params.nome){
			if (verbose){
				println ("NOMe-seq outfile basename: $outfile_basename")
			}
			cov2cyt_options += " --nome"
		}

		
		if (verbose){
			println ("Now running command: coverage2cytosine --genome $genome $cov2cyt_options --output ${outfile_basename} $coverage_file ")
		}

		"""
		module load bismark
		coverage2cytosine --genome $genome $cov2cyt_options --output ${outfile_basename} $coverage_file
		"""
		
		
}