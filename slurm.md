### How to use slurm ?
Nextflow is designed to work on many executors such as SGE, SLURM, ... or even on clouds such as Kubernates, Amazon, ...

On our cluster, we have the SLURM batch scheduler. To active it, create a file named `nextflow.config` in current directory and write the following lines:

```
process.executor = 'slurm'
```

### Run the workflow
```
nextflow run nextflow-io/hello
```

In another terminal you can check the execution with :
```
squeue -u <username>
```

## Load multiple modules
```
module = 'bioinfo/bwa-0.7.15:bioinfo/samtools-1.8'

# OR

module = ['bioinfo/bwa-0.7.15','bioinfo/samtools-1.8']
```

## ceres.config from IOWA state
```
process {
  executor = 'slurm'
  clusterOptions =  '-N 1 -n 16 -t 02:00:00'
  withLabel: blast { module = 'blast+' }
  withLabel: software_check { module = 'blast+:parallel' }
}
```