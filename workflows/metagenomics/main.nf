#!/usr/bin/env nextflow

// DSL 2 syntax
nextflow.preview.dsl=2

// parameters
params.help = false
params.read_path  = "${workflow.projectDir}/data"

// parameters decont
params.decont_refpath = '/data/nucleotide/'
params.decont_index   = 'hg19.fa'
params.decont_outdir  = './pipeline_output/decont_out'
ch_bwa_idx = file(params.decont_refpath)
// parameters kraken2                                              // ***
params.kraken2_refpath = '/data/minikraken2_v2_8GB_201904_UPDATE/' // ***
params.kraken2_outdir = './pipeline_output/kraken2_out'            // ***
ch_kraken_idx = file(params.kraken2_refpath)                       // ***

include './decont' params(index: "$params.decont_index", outdir: "$params.decont_outdir")
include './kraken2' params(outdir: "$params.kraken2_outdir")       // ***


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

workflow{
    DECONT(ch_bwa_idx, ch_reads)
    KRAKEN2(ch_kraken_idx, DECONT.out[0])                          // ***
    BRACKEN(ch_kraken_idx, KRAKEN2.out[0], Channel.from('s', 'g')) // ***
}