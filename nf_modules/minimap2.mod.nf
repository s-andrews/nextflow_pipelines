nextflow.enable.dsl=2
params.local = ''
params.no_output = true
mm_opts = ""

// TODO: add in splice/no splice options
// check memory usage for full size files

process MINIMAP2 {
	
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	label 'hugeMem'
	label 'multiCore'
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		path (genome_index)
		val (splice)
		//val (verbose)

	output:
	    tuple val(name), path ("*sam"),        emit: sam

	publishDir "$outputdir",
		mode: "link", overwrite: true, enabled: !params.no_output

	script:

		if (splice) {
			mm_opts = "-x splice"
		} else {
			mm_opts = "-x map-ont" // This is the default option for -x so we don't strictly need this.
		}

		index = genome_index.baseName
		outname = name + "_" + index + ".sam"

		"""
		module load minimap
		minimap2 -t 32 --sam-hit-only -o ${outname} -a ${mm_opts} -uf ${genome_index} ${reads} 		
		"""

}

//minimap2 -t 32 --sam-hit-only -o ${name}.sam -a ${mm_opts} -uf ${genome_index} ${reads} 		
// trying it without the y
// minimap2 -y -t 32 --sam-hit-only -o ${name}.sam -ax splice -uf ${genome_index} ${reads}
