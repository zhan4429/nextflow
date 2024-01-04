Here's how to use both `-profile` and `-c` when running nf-core pipelines:

## Understand their roles:

`-profile`: Specifies a configuration profile that defines settings for memory, CPU, and other resources. This optimizes pipeline execution for different environments (e.g., local, cluster, cloud).
`-c`: Overrides specific configuration parameters within the chosen profile. This allows fine-tuning without modifying the profile itself.

## Order of usage
Always use `-profile` before `-c`. This ensures the profile is loaded first, providing a base for any overrides.

Example usage:
```
nextflow run nf-core/<pipeline> -profile <profile_name> -c config.configName=value
```

