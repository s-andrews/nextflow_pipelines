#!/usr/bin/env nextflow
nextflow.preview.dsl=2


def getGenome(name) {

    // Find a file with the same name as the genome
    // in our genomes.d directory

    scriptDir = workflow.projectDir

    genomeFH = new File(scriptDir.toString()+"/genomes.d/"+name+".genome").newInputStream()

    genomeValues = ["name" : name]

    genomeFH.eachLine {
        sections =  it.split("\\s+",2)
        genomeValues[sections[0]] = sections[1]
    }

    return genomeValues

}

