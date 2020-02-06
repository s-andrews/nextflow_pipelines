#!/usr/bin/env nextflow
nextflow.preview.dsl=2


def getGenome(name) {

    // Find a file with the same name as the genome
    // in our genomes.d directory

    // TODO: Figure out where our script is so we can 
    // look relative to that. Kludge it for now

    scriptDir = "/bi/home/andrewss/nextflow_pipelines"

    genomeFH = new File(scriptDir+"/genomes.d/"+name+".genome").newInputStream()

    genomeValues = ["name" : name]

    genomeFH.eachLine {
        sections =  it.split("\\s+",2)
        genomeValues[sections[0]] = sections[1]
    }

    return genomeValues

}

