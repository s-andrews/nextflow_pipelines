nextflow.enable.dsl=2

process COVERAGE2CYTOSINE {
	tag "$coverage_file" // Adds name to job submission instead of (1), (2) etc.

	// dynamic directive
	memory { 20.GB * task.attempt }
	errorStrategy { sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry' }
	maxRetries 5

	publishDir { outputdir },
		mode: "link", overwrite: true

    input:
	    path(coverage_file)
		val (outputdir)
		val (coverage2cytosine_args)
		val (verbose)
		val (genome)
		val (nome)

	output:
	    path "*{report.txt.gz,report.txt}", emit: report
		path "*{.cov.gz,.cov}",             emit: coverage
		path "*cytosine_context_summary.txt", optional: true, emit: summary

	script:

		def bismark_genome = genome["bismark"]

		// removing the file extension from the input file name
		// (https://www.nextflow.io/docs/latest/script.html#removing-part-of-a-string)
		def outfile_basename = coverage_file.toString()  // Important to convert nextflow.processor.TaskPath object to String first
		outfile_basename = (outfile_basename - ~/.bismark.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov.gz$/)
		outfile_basename = (outfile_basename - ~/.cov$/)

		if (verbose){
			println ("[MODULE] BISMARK COVERAGE2CYTOSINE ARGS: " + coverage2cytosine_args)
			println ("Bismark Genome is: " + bismark_genome)
		}

		// Options we add are
		def cov2cyt_options = coverage2cytosine_args + " --gzip "

		if (nome){
			if (verbose){
				println ("NOMe-seq outfile basename: $outfile_basename")
			}
			cov2cyt_options += " --nome"
		}


		if (verbose){
			println ("Now running command: coverage2cytosine --genome $bismark_genome $cov2cyt_options --output ${outfile_basename} $coverage_file ")
		}

		"""
		module load bismark
		coverage2cytosine --genome $bismark_genome $cov2cyt_options --output ${outfile_basename} $coverage_file
		"""


}
