#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.files = args

def getFileBaseNames(fileList) {

    baseNames = [:]

    for (String s : fileList) {
        matcher = s =~ /^(.*)_([1234]).(fastq|fq).gz$/

        if (matcher.matches()) {
            if (! baseNames.containsKey(matcher[0][1])) {
                baseNames[matcher[0][1]] = []
            }
            baseNames[matcher[0][1]].add(matcher[0][2])
        }

    }

    patterns = []
    for (s in baseNames) {
        pattern = s.key+"_{"+s.value.join(",")+"}.{fastq,fq}.gz"
        patterns.add(pattern)
    }

    return(patterns)
}

file_ch = Channel.fromFilePairs(
    getFileBaseNames(params.files),
    ,size:-1)

file_ch.subscribe {println("Found $it")}