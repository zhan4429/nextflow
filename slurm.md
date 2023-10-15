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