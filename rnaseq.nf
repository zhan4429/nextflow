#!/usr/bin/env nextflow
params.reads = "$baseDir/input/fastq/ggal_gut_{1,2}.fq"
params.readlength = 100
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"
params.STAR_genome_index_dir=/scratch/brown/zhan4429/reference/STAR

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
 
    INDEX(params.transcriptome)
    FASTQC(read_pairs_ch)
    QUANT(INDEX.out, read_pairs_ch)
}

process INDEX {
  input:
    path transcriptome from params.transcriptome_fa
  """
  salmon index -t ${transcriptome} -i transcripts_index
  """
  output:
    path 'transcripts_index'
}

process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir params.outdir
 
    input:
    tuple val(sample_id), path(reads)
 
    output:
    path "fastqc_${sample_id}_logs"
 
    script:
    """
    fastqc.sh "$sample_id" "$reads"
    fastqc -o $fastqc_out -t $threads $FASTQ1 $FASTQ2
    """
}