#!/usr/bin/env nextflow
nextflow.preview.dsl=2


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
        matcher = s =~ /^(.*)_(R?[1234]).(fastq|fq).gz$/

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
