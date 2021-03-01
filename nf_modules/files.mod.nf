#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def makeFilesChannel(fileList) {    
    
    // def meta = [:]
    file_ch = Channel.fromFilePairs( getFileBaseNames(fileList), size:-1)
        //.map { it -> [ [meta.id = it[0]], it[1]] }
            
        // .map { meta.id = it[0]}
        //.map { meta.files = it[1]}
        //.view()
        // .subscribe onNext: { println it }, onComplete: { println 'Done' }

    // meta.each { key, val -> 
    //   println ("Key: $key = Files: $val")
    // }
    return(file_ch)

    // TODO: changing the input meta-data to a data structure that will be available throughout the workflow.
    // Inspired by NF-Core and Harshil
// // Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
// def create_fastq_channels(LinkedHashMap row) {
//     def meta = [:]
//     meta.id           = row.sample
//     meta.single_end   = row.single_end.toBoolean()

//     def array = []
//     if (!file(row.fastq_1).exists()) {
//         exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
//     }
//     if (meta.single_end) {
//         array = [ meta, [ file(row.fastq_1) ] ]
//     } else {
//         if (!file(row.fastq_2).exists()) {
//             exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
//         }
//         array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
//     }
//     return array    
// }

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
        // println("$pattern")
    }
    for (s in bareFiles) {
        pattern = s+".{fastq,fq}.gz"
        patterns.add(pattern)
    }

    return(patterns)
}
