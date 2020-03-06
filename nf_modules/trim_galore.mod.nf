nextflow.preview.dsl=2

// params.trim_galore_args = ''
params.singlecell = ''
params.rrbs = ''
params.pbat = ''


process TRIM_GALORE {	
    
	input:
	    tuple val (name), path (reads)
		val (outputdir)
		val (trim_galore_args)
		val (verbose)

	output:
	    tuple val(name), path ("*fq.gz"), emit: reads
		path "*trimming_report.txt", optional: true, emit: report
		
	publishDir "$outputdir",
		mode: "link", overwrite: true


    script:
		if (verbose){
			println ("[MODULE] TRIM GALORE ARGS: " + trim_galore_args)
		}
		
		// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
		// This is only a temporary workaround until Paolo has fixed the Nextflow bug.
		// https://github.com/nextflow-io/nextflow/issues/1519
		// trim_galore_args = params.trim_galore_args.replaceAll(/'/,"")
		trim_galore_args = trim_galore_args.replaceAll(/'/,"")


		pairedString = ""
		if (reads instanceof List) {
			pairedString = "--paired"
		}

		if (params.singlecell){
			trim_galore_args = trim_galore_args + " --clip_r1 6 "
			if (pairedString == "--paired"){
				trim_galore_args = trim_galore_args + " --clip_r2 6 "
			}
		}

		if (params.rrbs){
			trim_galore_args = trim_galore_args + " --rrbs "
		}
		
		if  (params.pbat){
			trim_galore_args = trim_galore_args + " --clip_r1 $params.pbat "
			if (pairedString == "--paired"){
				trim_galore_args = trim_galore_args + " --clip_r2 $params.pbat "
			}
		}

		"""
		module load trim_galore
		trim_galore $trim_galore_args ${pairedString} ${reads}
		"""

}