#!/usr/bin/env nextflow

// DSL 2 syntax
nextflow.preview.dsl=2

// parameters
params.help = false
params.read_path  = "${workflow.projectDir}/data"

// help message
def helpMessage() {
    log.info"""
    =================================================================
    Usage: ${workflow.projectDir}/main.nf --read_path PATH/OF/READS
    =================================================================
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}

// Create channel for reads
ch_reads = Channel
    .fromFilePairs(params.read_path + '/**{1,2}.f*q*', flat: true)

ch_reads.view()