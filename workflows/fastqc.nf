#!/usr/bin/env nextflow

process fastqc{
    input:
    path input

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -q $input
    """

    workflow {
        Channel.fromPath("*.fastq.gz") | fastqc
    }
}