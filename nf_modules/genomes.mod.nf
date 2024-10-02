#!/usr/bin/env nextflow
nextflow.enable.dsl=2


def getGenome(name) {

    // Find a file with the same name as the genome in our genomes.d directory

    scriptDir = workflow.projectDir
    
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

    genomeFH = new File (fileName).newInputStream()

    genomeValues = [:]  // initialising map. name is also part of each .genome file

    genomeFH.eachLine {
        sections =  it.split("\\s+",2)
        genomeValues[sections[0]] = sections[1]
    }

    return genomeValues

}

def listGenomes(){
    
    println ("These genomes are currently available to choose from:")
    println ("=====================================================")
    scriptDir = workflow.projectDir + "/genomes.d/"
    // println (scriptDir) // last slash is consumed
    allFiles = scriptDir.list()
    
    for( def file : allFiles.sort() ) {
        
        if( file =~ /.genome$/){

            genomeFH = new File(scriptDir.toString() + "/$file").newInputStream()
            name = file.replaceFirst(/.genome/, "")
        
            println (name)
            genomeFH.eachLine {
                if (params.verbose){
                    println ("\t$it")
                }
            }
        }
    }

    // We'll repeat this for the genomes.d directory in the users home directory 
    scriptDir = new File(System.getProperty("user.home") + "/genomes.d/")
    // println (scriptDir) // last slash is consumed

    if (scriptDir.exists()) {
        allFiles = scriptDir.list()
        
        for( def file : allFiles.sort() ) {
            
            if( file =~ /.genome$/){

                genomeFH = new File(scriptDir.toString() + "/$file").newInputStream()
                name = file.replaceFirst(/.genome/, "")
            
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

