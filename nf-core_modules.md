# Module
## Terminology
Module:
- Atomic process
- Can be used within pipelines
- e.g. a module file containing a single tool such as FastQC.

```
nf-core modules
```

## Install tools
```
nf-core modules install --help
```

## Find a module
```
nf-core modules list
nf-core modules list remote
nf-core modules list local
```

## Install a module
```
nf-core modules install [OPTIONS] <pipeline directory> <tool name>
```

## Remove a module
```
nf-core modules remove [OPTIONS] <pipeline directory> <tool name>
nf-core modules remove . fastqc
```

## Check the modules
Check for upstream updates
```
nf-core modules lint [OPTIONS] <pipeline directory> <tool name>
nf-core modules lint . fastqc
```

## Update the modules
Check for upstream updates
```
nf-core update [OPTIONS] <pipeline directory> <tool name>
nf-core modules update . fastqc
```

## Use a module
1. Install module or create local one.
2. Adapt conf/modules.config:
    a. Customize tool parameters
        i. args (Tool1)
        ii args2 (Tool 2)
```
params {
   modules {
       'picard_markduplicates' {
} }
```
3. Add includeConfig 'conf/modules.config' to nextflow.config 
4. `INCLUDE` module into (sub-)workflows.

## Combine subworkflows to workflows
