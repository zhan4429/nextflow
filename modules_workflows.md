## Modules
A process that can be used within different pipelines and is as atomic as possible i.e. cannot be split into another module.
e.g. a module file containing the process definition for a single tool such as FASTQC.
## Sub-workflow
A chain of multiple modules that offer a higher-level of functonality within the context of a pipeline.
e.g. a sub-workflow to sort, index and run some basic stats on a BAM file.

## Workflow
An end-to-end pipeline - this can either be implemented using a large monolithic script as wthin Nextflow DSL1, or by using a combinatoin of Nextflow DSL2 individual modules and sub-workflows. 
e.g. from one or more inputs to a series of final outputs.
