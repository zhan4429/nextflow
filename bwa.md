## align
```
process align_sample{
    input:
    file 'reference.fa' from genome_ch
    file 'sample.fq' from reads_ch

    output:
    file 'sample.bam' into bam_ch

    script:
    """
    bwa mem reference.fa sample.fq \
        | samtools sort -o sample.bam
    """
}

```


## index
```
process index_sample {

    input:
    file 'sample.bam' from bam_ch
    
    output:
    file 'sample.bai' into bai_ch

    script:
    """
    samtools index sample.bam
    """
}
```