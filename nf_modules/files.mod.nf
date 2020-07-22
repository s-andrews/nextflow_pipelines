#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.single_end = false

def makeFilesChannel(fileList) {    
    
    file_ch = Channel.fromFilePairs(
                getFileBaseNames(fileList),
                ,size:-1
                )

    return(file_ch)
}

def getFileBaseNames(fileList) {

    baseNames = [:]

    bareFiles = []

    for (String s : fileList) {
       
        if (params.single_end){
            matcher = s =~ /^(.*).(fastq|fq).gz$/

            if (matcher.matches()) {
                bareFiles.add(matcher[0][1])
            }
        }
        else{

            // let's make a distinction for paired-end files already trimmed with Trim Galore
            // Paired-end files trimmed with Trim Galore follow the folling pattern:
            // lane1_TGGTTGTT_small_test_L001_R1_val_1.fq.gz
            // lane1_TGGTTGTT_small_test_L001_R3_val_2.fq.gz

            if (s =~ /_val_/){
                // println ("Input file '$s' looks like a Trim Galore paired-end file")
                matcher = s =~ /^(.*)_(R?[1234])_val_[12].(fastq|fq).gz$/
                // in the above example, this would identify the following basename:
                // lane1_TGGTTGTT_small_test_L001
                // println (matcher[0])
                if (matcher.matches()) {
                    if (! baseNames.containsKey(matcher[0][1])) {
                        baseNames[matcher[0][1]] = []
                    }
                    baseNames[matcher[0][1]].add(matcher[0][2])
                }
                else {
                    matcher = s =~ /^(.*).(fastq|fq).gz$/

                    if (matcher.matches()) {
                        bareFiles.add(matcher[0][1])
                    }
                }
            
            }
            else{ // not Trim Galore processed
                matcher = s =~ /^(.*)_(R?[1234]).(fastq|fq).gz$/

                // println (matcher[0])
                if (matcher.matches()) {
                    if (! baseNames.containsKey(matcher[0][1])) {
                        baseNames[matcher[0][1]] = []
                    }
                    baseNames[matcher[0][1]].add(matcher[0][2])
                }
                else {
                    matcher = s =~ /^(.*).(fastq|fq).gz$/

                    if (matcher.matches()) {
                        bareFiles.add(matcher[0][1])
                    }
                }
            }
        }

    }

    patterns = []
    for (s in baseNames) {
        pattern = s.key+"_{"+s.value.join(",")+"}.{fastq,fq}.gz"
        patterns.add(pattern)
    }
    for (s in bareFiles) {
        pattern = s+".{fastq,fq}.gz"
        patterns.add(pattern)
    }


    return(patterns)
}
