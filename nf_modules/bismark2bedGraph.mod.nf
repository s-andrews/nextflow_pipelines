nextflow.enable.dsl=2

process BISMARK2BEDGRAPH {

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'bigMem' // 20G

	publishDir { outputdir },
		mode: "link", overwrite: true

    input:
	    tuple val (name), path(reads)
		val (outputdir)
		val (bismark2bedGraph_args)
		val (verbose)
		val (dirty_harry)

	output:
	    path "*cov.gz",        emit: coverage
		path "*bedGraph.gz",   emit: bedGraph

    script:

		if (verbose){
			println ("[MODULE] BISMARK2BEDGRAPH ARGS: " + bismark2bedGraph_args)
		}

		// Options we add are
		def bismark2bedGraph_options = bismark2bedGraph_args

		def output_name
		if (dirty_harry){
			output_name = name + "_DH.bedGraph.gz"  // Dirty Harry
		}
		else{
			output_name = name + ".bedGraph.gz"
		}
		// println ("Output name: $output_name")
		// println ("Input names: $reads")

		def all_reads = reads

		"""
		module load bismark
		bismark2bedGraph --buffer 15G -o $output_name $bismark2bedGraph_options $all_reads
		"""

}
