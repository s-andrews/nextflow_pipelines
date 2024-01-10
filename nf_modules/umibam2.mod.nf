nextflow.enable.dsl=2

params.dual = false

// This only works for single ended data

process UMIBAM2 {	
    
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	// dynamic directive to increase memory as required
	cpus = 1
	memory { 20.GB * task.attempt }  
	errorStrategy { sleep(Math.pow(2, task.attempt) * 30 as long); return 'retry' }
  	maxRetries 5
	  
	input:
	   // tuple val(name), path(bam)
	    path(bam)
		path(bai) // an index file is required even if it's not used directly
		val (outputdir)
		val (umibam_args)
		val (verbose)

	output:
		path "*report.txt", emit: report
		//tuple val(name), path ("*bam"),        emit: bam
		path ("*bam"),        emit: bam

	publishDir "$outputdir",
		mode: "link", overwrite: true


    script:
		if (verbose){
			println ("[MODULE] UMIBAM ARGS: " + umibam_args)
		}
			
		"""
		module load python3
		python /bi/apps/TrAELseq/latest/TrAEL-seq/umibam2.py $bam
		"""
		
		// The output files should be renamed so that they bismark2report picks up everything
		
		// renaming files using Bash
		// for f in *UMI_dedup* ; do mv "\$f" "\${f/UMI_/}" ; done
		
		// renaming using rename (works on our cluster)
		// rename UMI_d d *

		// A third option should be the saveAs directive (https://www.nextflow.io/docs/latest/process.html#publishdir)
		// unclear to me at the moment though how this would work exactly
}