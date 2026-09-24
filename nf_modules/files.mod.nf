#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def makeFilesChannel(fileList) {

    def patterns = getFileBaseNames(fileList)

    // Channel.fromFilePairs is called with size:-1 because a single call has to handle both
    // single-end patterns (1 file) and paired-end patterns (2+ files) together. That disables
    // Nextflow's own "expected N files" check, so a pattern that only matches some of its files
    // (e.g. a transiently-unreadable symlink on a busy network filesystem) is silently accepted
    // as a short list rather than raising an error -- which then surfaces much later downstream
    // as a cryptic failure (e.g. a 'null' file being fed into whatever runs next). We know from
    // the pattern itself how many files it should have matched, so check that explicitly here
    // and fail loudly and immediately if it doesn't line up.
    def expectedCounts = [:]
    patterns.each { pattern ->
        def m = pattern =~ /^(.*)_\{([^}]*)\}\./
        if (m.find()) {
            expectedCounts[m.group(1)] = m.group(2).split(',').length
        }
    }

    // def meta = [:]
    def file_ch = Channel.fromFilePairs( patterns, size:-1)
        .map { key, files ->
            def expected = expectedCounts[key]
            if (expected != null && files.size() != expected) {
                throw new IllegalStateException(
                    "[FILE PAIRING ERROR] Expected ${expected} file(s) for '${key}' but only found " +
                    "${files.size()}: ${files}. This usually means one of a set of paired fastq files " +
                    "was missing or unreadable when the pipeline started (e.g. a broken symlink, or a " +
                    "slow/stale network filesystem) -- please check the input files for this sample and re-run."
                )
            }
            [key, files]
        }
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

    def baseNames = [:]

    def bareFiles = []

    fileList.each { s ->

        def matcher

        if (params.single_end){
            matcher = s =~ /^(.*).(fastq|fq).gz$/

            if (matcher.matches()) {
                bareFiles.add(matcher[0][1])
            }
        }
        else{

            // let's make a distinction for paired-end files already trimmed with Trim Galore
            // Paired-end files trimmed with Trim Galore follow the following pattern:
            // lane1_TGGTTGTT_small_test_L001_R1_val_1.fq.gz
            // lane1_TGGTTGTT_small_test_L001_R3_val_2.fq.gz

            if (s =~ /_val_/){
                // println ("Input file '$s' looks like a Trim Galore paired-end file")
                matcher = s =~ /^(.*)_(R?[1234]_val_[12]).(fastq|fq).gz$/
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

    def patterns = []
    baseNames.each { s ->
        // println (s)
        def pattern = s.key+"_{"+s.value.join(",")+"}.{fastq,fq}.gz"
        patterns.add(pattern)
        // println("$pattern")
    }
    bareFiles.each { s ->
        def pattern = s+".{fastq,fq}.gz"
        patterns.add(pattern)
    }

    return(patterns)
}
