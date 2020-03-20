nextflow.preview.dsl=2

params.singlecell = ''
params.rrbs       = ''
params.pbat       = ''
params.clock      = false
// For Epigenetic Clock Processing
params.three_prime_clip_R1 = ''
params.three_prime_clip_R2 = ''


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
		
		pairedString = ""
		if (reads instanceof List) {
			pairedString = "--paired"
		}
		
		// Specialised Epigenetic Clock Processing		
		if (params.clock){
			trim_galore_args += " --breitling "	
		}
		else{
			if (params.singlecell){
				trim_galore_args += " --clip_r1 6 "
				if (pairedString == "--paired"){
					trim_galore_args += " --clip_r2 6 "
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

			// Second step of Clock processing:
			if (params.three_prime_clip_R1 && params.three_prime_clip_R2){
				trim_galore_args +=	" --three_prime_clip_R1 ${params.three_prime_clip_R1} --three_prime_clip_R2 ${params.three_prime_clip_R2} "
			}
		}

		"""
		module load trim_galore
		trim_galore $trim_galore_args ${pairedString} ${reads}
		"""

}