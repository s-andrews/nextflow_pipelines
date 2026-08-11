#!/usr/bin/env nextflow
nextflow.enable.dsl=2


def getGenome(name) {

    // Find a file with the same name as the genome in our genomes.d directory

    def scriptDir = workflow.projectDir

    // die gracefully if the user specified an incorrect genome
    def fileName = scriptDir.toString() + "/genomes.d/" + name + ".genome"
    def testFile = new File(fileName)
    if (!testFile.exists()) {

        // We'll try the users home genomes directory

        scriptDir = new File(System.getProperty("user.home") + "/genomes.d/")

        // die gracefully if the user specified an incorrect genome
        fileName = scriptDir.toString() + "/" + name + ".genome"
        testFile = new File(fileName)
        if (!testFile.exists()) {
            println("\nFile >>$fileName<< does not exist. Listing available genomes...\n")
            listGenomes()
        }
    }

    def genomeFH = new File (fileName).newInputStream()

    def genomeValues = [:]  // initialising map. name is also part of each .genome file

    genomeFH.eachLine {
        def sections =  it.split("\\s+",2)
        genomeValues[sections[0]] = sections[1]
    }

    return genomeValues

}

// Names of every available genome (both the shared genomes.d and the user's own ~/genomes.d),
// without the .genome suffix. Shared by listGenomes() and resolveGenome().
def availableGenomeNames() {

    def names = []

    def scriptDir = new File(workflow.projectDir.toString() + "/genomes.d/")
    if (scriptDir.exists()) {
        scriptDir.list().each { file ->
            if (file =~ /\.genome$/) {
                names.add(file.replaceFirst(/\.genome$/, ""))
            }
        }
    }

    def homeScriptDir = new File(System.getProperty("user.home") + "/genomes.d/")
    if (homeScriptDir.exists()) {
        homeScriptDir.list().each { file ->
            if (file =~ /\.genome$/) {
                names.add(file.replaceFirst(/\.genome$/, ""))
            }
        }
    }

    return names.unique().sort()
}

def listGenomes(){

    println ("These genomes are currently available to choose from:")
    println ("=====================================================")
    def scriptDir = workflow.projectDir + "/genomes.d/"
    // println (scriptDir) // last slash is consumed
    def allFiles = scriptDir.list()

    allFiles.sort().each { file ->

        if( file =~ /.genome$/){

            def genomeFH = new File(scriptDir.toString() + "/$file").newInputStream()
            def name = file.replaceFirst(/.genome/, "")

            println (name)
            genomeFH.eachLine {
                if (params.verbose){
                    println ("\t$it")
                }
            }
        }
    }

    // We'll repeat this for the genomes.d directory in the users home directory
    def homeScriptDir = new File(System.getProperty("user.home") + "/genomes.d/")
    // println (homeScriptDir) // last slash is consumed

    if (homeScriptDir.exists()) {
        def homeFiles = homeScriptDir.list()

        homeFiles.sort().each { file ->

            if( file =~ /.genome$/){

                def genomeFH = new File(homeScriptDir.toString() + "/$file").newInputStream()
                def name = file.replaceFirst(/.genome/, "")

                println (name)
                genomeFH.eachLine {
                    if (params.verbose){
                        println ("\t$it")
                    }
                }
            }
        }
    }


    println ("\nTo see this list of available genomes with more detailed information about paths and indexes,\nplease re-run the command including '--list_genomes --verbose'\n\n")

    System.exit(1)
}

// Standard Levenshtein (edit) distance between two strings, case-insensitive.
def levenshteinDistance(String a, String b) {

    a = a.toLowerCase()
    b = b.toLowerCase()

    def m = a.length()
    def n = b.length()
    def dist = (0..m).collect { (0..n).collect { 0 } }

    (m + 1).times { i -> dist[i][0] = i }
    (n + 1).times { j -> dist[0][j] = j }

    m.times { di ->
        def i = di + 1
        n.times { dj ->
            def j = dj + 1
            def cost = (a.charAt(i - 1) == b.charAt(j - 1)) ? 0 : 1
            dist[i][j] = Math.min(
                Math.min(dist[i - 1][j] + 1, dist[i][j - 1] + 1),
                dist[i - 1][j - 1] + cost
            )
        }
    }

    return dist[m][n]
}

/* Resolve a (possibly mistyped) genome name against genomes.d.
 *
 *  1. Exact match -> used directly.
 *  2. An rDNA<->rRNA swap that produces an exact match -> auto-accepted (this specific typo is
 *     common enough in Houseley group sample sheets that we don't bother the user about it).
 *  3. Otherwise, the closest available name(s) by edit distance are offered: if an interactive
 *     terminal is attached the user is asked to confirm; otherwise (e.g. running with -bg, or
 *     from a script) we can't safely guess, so it's treated as unresolved.
 *  4. No close-enough match at all -> unresolved.
 *
 * Returns [resolved: true/false, name: <the genome name that was actually used, or null>,
 *          genome: <the loaded genome map, or null>, message: <human-readable explanation>]
 */
def resolveGenome(String requestedName) {

    def available = availableGenomeNames()

    if (available.contains(requestedName)) {
        return [resolved: true, name: requestedName, genome: getGenome(requestedName),
                message: "'${requestedName}' found."]
    }

    // rDNA <-> rRNA is a very common substitution in this group's sample sheets
    def swapped = requestedName.contains("rDNA") ? requestedName.replace("rDNA", "rRNA")
                : requestedName.contains("rRNA") ? requestedName.replace("rRNA", "rDNA")
                : null

    if (swapped && available.contains(swapped)) {
        println("[GENOME] '${requestedName}' not found -- using '${swapped}' instead (rDNA/rRNA auto-correction).")
        return [resolved: true, name: swapped, genome: getGenome(swapped),
                message: "'${requestedName}' auto-corrected to '${swapped}'."]
    }

    // Otherwise look for the closest name(s) by edit distance. Genome names are conventionally
    // "<build>_<extras...>" (e.g. GRCm38_rRNA_lambda) -- so two extra signals beat raw edit
    // distance on its own: (1) a candidate sharing the requested name's leading "<build>" token,
    // and (2) a candidate that retains ALL of the requested name's other tokens too (e.g. for
    // "GRCm38_lambda", prefer "GRCm38_rRNA_lambda" -- which keeps both "GRCm38" and "lambda" --
    // over "GRCm38_rRNA", which drops "lambda", even though the latter has a smaller raw edit
    // distance; and don't let "GRCh38_T4_lambda" outscore either just for being shorter).
    def requestedTokens = requestedName.toLowerCase().split('_')
    def requestedBuild = requestedTokens[0]
    def threshold = Math.max(6, (requestedName.length() * 0.4).intValue())
    def scored = available.collect { candidate ->
        def candidateTokens = candidate.toLowerCase().split('_') as List
        def retainedCount = requestedTokens.findAll { candidateTokens.contains(it) }.size()
        [name: candidate, distance: levenshteinDistance(requestedName, candidate),
         sameBuild: candidateTokens.contains(requestedBuild), retainedCount: retainedCount]
    }
    .findAll { it.distance <= threshold }
    .sort { [it.sameBuild ? 0 : 1, -it.retainedCount, it.distance] }

    if (!scored) {
        return [resolved: false, name: null, genome: null,
                message: "'${requestedName}' not found in genomes.d, and no close match was found."]
    }

    def best = scored[0].name
    def console = System.console()

    if (console == null) {
        return [resolved: false, name: null, genome: null,
                message: "'${requestedName}' not found in genomes.d. Closest match: '${best}' -- " +
                         "re-run interactively (not with -bg) to confirm, or fix the sample sheet."]
    }

    def answer = console.readLine("[GENOME] '${requestedName}' not found. Did you mean '${best}'? [y/N]: ")
    if (answer?.trim()?.toLowerCase() in ["y", "yes"]) {
        return [resolved: true, name: best, genome: getGenome(best),
                message: "'${requestedName}' confirmed by user as '${best}'."]
    }

    return [resolved: false, name: null, genome: null,
            message: "'${requestedName}' not found in genomes.d; suggested '${best}' was declined."]
}
