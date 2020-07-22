nextflow.enable.dsl=2

params.dual = false

process UMIBAM {	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	input:
	    path(bam)
		val (outputdir)
		val (umibam_args)
		val (verbose)

	output:
		path "*report.txt", emit: report
		path "*bam",        emit: bam

	publishDir "$outputdir",
		mode: "link", overwrite: true


    script:
		if (verbose){
			println ("[MODULE] UMIBAM ARGS: " + umibam_args)
		}
		
		// Specialised Epigenetic Clock Processing		
		if (params.dual){
			umibam_args += " --double_umi "	
		}
		
		"""
		module load UmiBam
		
		UmiBam $umibam_args  $bam
		
		rename UMI_d d *
		"""
		
		// The output files should be renamed so that they bismark2report picks up everything
		
		// renaming files using Bash
		// for f in *UMI_dedup* ; do mv "\$f" "\${f/UMI_/}" ; done
		
		// renaming using rename (works on our cluster)
		// rename UMI_d d *

		// A third option should be the saveAs directive (https://www.nextflow.io/docs/latest/process.html#publishdir)
		// unclear to me at the moment though how this would work exactly
}