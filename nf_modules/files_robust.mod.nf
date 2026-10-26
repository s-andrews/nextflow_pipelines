nextflow.enable.dsl=2

// Houseley-only alternative to makeFilesChannel() in files.mod.nf.
//
// makeFilesChannel() hands off to Channel.fromFilePairs(), which re-derives R1/R2 pairs by
// globbing the filesystem itself rather than trusting the order/completeness of the file list
// it was given. In real Houseley runs (large groups, many samples dispatched together) this has
// repeatedly and reproducibly dropped one read of a pair for a specific sample, silently handing
// a 1-file "pair" downstream instead of erroring. This version never re-globs anything -- it
// builds pairs directly from the filenames it's given, so pairing can't depend on filesystem
// glob timing/ordering at all. Only used by pipelines dispatched from nf_houseley_run_processing
// (opted in via --robust_file_pairing), so behaviour for every other pipeline/user is unchanged.

def makeFilesChannelRobust(fileList) {

    def groups = groupFilesRobust(fileList)

    def pairs = groups.collect { key, readMap ->
        // Sort by read tag (R1 < R2 < R3 < R4, or R1_val_1 < R2_val_2) so files always end up
        // in the right order regardless of what order they arrived in on the command line.
        def orderedFiles = readMap.keySet().sort().collect { readTag -> file(readMap[readTag]) }
        // A key with only one file (a genuinely single-end sample -- either a bare filename with
        // no R1/R2 tag, or an R1 with no matching R2, e.g. a Houseley run where the sequencer only
        // produced R1 for every sample) must come out as a bare path, not a 1-element list.
        // Downstream modules decide single-end vs paired-end with 'reads instanceof List', which a
        // 1-element list still satisfies -- that wrongly sent single-end samples down the paired
        // code path (e.g. Bowtie2's '-1 file -2 <nothing>'), rather than falling back to a plain
        // single-end command.
        [key, orderedFiles.size() == 1 ? orderedFiles[0] : orderedFiles]
    }

    return Channel.fromList(pairs)
}

// Mirrors the classification regexes in files.mod.nf's getFileBaseNames(), but records the
// actual filename against each read tag (not just the tag), since pairing is finished here
// rather than being handed off to a filesystem glob afterwards.
def groupFilesRobust(fileList) {

    def groups = [:]     // key -> [readTag: filename]
    def bareFiles = []

    fileList.each { s ->

        // Match against the basename only, so a directory prefix in the input (e.g. dispatching
        // on 'trim_galore/sample_R1.fq.gz' after the Houseley Vazyme pre-trim step) doesn't leak
        // into the pairing key. Downstream processes use that key as the sample "name" to build
        // their own output filenames -- a key with a '/' in it makes them try to write into a
        // subdirectory that doesn't exist in their own isolated work directory. The full path
        // (with any prefix) is still what's stored as the value, since that's what's needed to
        // actually resolve/stage the file.
        def baseName = s.contains('/') ? s.substring(s.lastIndexOf('/') + 1) : s

        def matcher

        if (params.single_end) {
            matcher = baseName =~ /^(.*)\.(fastq|fq)\.gz$/
            if (matcher.matches()) {
                bareFiles.add(s)
            }
            return
        }

        if (baseName =~ /_val_/) {
            matcher = baseName =~ /^(.*)_(R?[1234]_val_[12])\.(fastq|fq)\.gz$/
        } else {
            matcher = baseName =~ /^(.*)_(R?[1234])\.(fastq|fq)\.gz$/
        }

        if (matcher.matches()) {
            def key = matcher[0][1]
            def readTag = matcher[0][2]
            if (! groups.containsKey(key)) {
                groups[key] = [:]
            }
            groups[key][readTag] = s
        } else {
            matcher = baseName =~ /^(.*)\.(fastq|fq)\.gz$/
            if (matcher.matches()) {
                bareFiles.add(s)
            }
        }
    }

    bareFiles.each { s ->
        def baseName = s.contains('/') ? s.substring(s.lastIndexOf('/') + 1) : s
        def key = baseName.replaceFirst(/\.(fastq|fq)\.gz$/, "")
        groups[key] = ["": s]
    }

    return groups
}
