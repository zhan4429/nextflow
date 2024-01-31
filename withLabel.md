## withLabel
The `withLabel` selectors allow the configuration of all processes annotated with a label directive as shown below:

```
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
        queue = 'long'
    }
}
```

The above configuration example assigns 16 cpus, 64 Gb of memory and the long queue to all processes annotated with the `big_mem` label.

## label
The `label` directive allows the annotation of processes with mnemonic identifier of your choice. For example:

process bigTask {
  label 'big_mem'

  '''
  <task script>
  '''
}


## Example in institution config file
```
process {
  executor = 'slurm'
  clusterOptions =  '-N 1 -n 16 -t 02:00:00'
  withLabel: blast { module = 'blast+' }
  withLabel: software_check { module = 'blast+:parallel' }
}
```

