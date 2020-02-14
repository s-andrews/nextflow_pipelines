nextflow.preview.dsl=2

params.trim_galore_args = ''
params.singlecell = ''
params.rrbs = ''
params.pbat = ''
params.verbose = false

// We need to replace single quotes in the arguments so that they are not getting passed in as a single string
trim_galore_args = params.trim_galore_args.replaceAll(/'/,"")
if (params.verbose){
	println ("[MODULE] TRIM GALORE ARGS: " + trim_galore_args)
}

process TRIM_GALORE {	
    
	input:
	    tuple val(name), path(reads)

	output:
	    tuple val(name), path ("*fq.gz"), emit: reads
		tuple val(name), path ("*trimming_report.txt"), emit: reports
		
    script:

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