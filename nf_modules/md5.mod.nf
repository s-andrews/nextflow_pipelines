nextflow.enable.dsl=2
params.no_output = false

process CALC_MD5{	
    
	tag "$file" // Adds name to job submission instead of (1), (2) etc.

	publishDir { outputdir },
		mode: "link", overwrite: true, enabled: !params.no_output

	input:
		path(file)
		val (outputdir)
		val (verbose)

	output:
		// path "*report.txt", emit: report
		path "*md5",        emit: md5


    script:
				
		// TODO: Find more elegant way to strip file ending of input BAM file

		"""
		md5sum $file > ${file}.md5 
    	"""
		
	
}

