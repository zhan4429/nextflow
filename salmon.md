## quant
```
process QUANT {
    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path pair_id

    script:
    """
    salmon quant -i $index \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o $pari_id
    """
}
```