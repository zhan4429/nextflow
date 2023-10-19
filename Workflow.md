## Workflow definition
We can connect processes to create our pipeline inside a workflow scope. The workflow scope starts with the keyword `workflow`, followed by an optional name and finally the workflow body delimited by curly brackets `{}`.

### Implicit workflow
A workflow definition which does not declare any name is assumed to be the main workflow, and it is implicitly executed. Therefore it’s the entry point of the workflow application.

### Invoking processes with a workflow

To combined multiple processes invoke them in the order they would appear in a workflow. When invoking a process with multiple inputs, provide them in the same order in which they are declared in the input block of the process.

For example:
```
//workflow_01.nf
nextflow.enable.dsl=2

process INDEX {
    input:
      path transcriptome
    output:
      path 'index'
    script:
      """
      salmon index -t $transcriptome -i index
      """
}

 process QUANT {
    input:
      each  path(index)
      tuple(val(pair_id), path(reads))
    output:
      path pair_id
    script:
      """
      salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
      """
}

workflow {
    transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz',checkIfExists: true)
    read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz',checkIfExists: true)

    //index process takes 1 input channel as a argument
    index_ch = INDEX(transcriptome_ch)

    //quant channel takes 2 input channels as arguments
    QUANT( index_ch, read_pairs_ch ).view()
}
```
In this example, the INDEX process is invoked first and the QUANT process second. The INDEX output channel, assigned to the variable index_ch, is passed as the first argument to the QUANT process. The read_pairs_ch channel is passed as the second argument.

### Process composition
Processes having matching input-output declaration can be composed so that the output of the first process is passed as input to the following process.

We can therefore rewrite the previous workflow example as follows:
```
[..truncated..]

workflow {
  transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
  read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')

  // pass INDEX process as a parameter to the QUANT process
  QUANT(INDEX(transcriptome_ch),read_pairs_ch ).view()
}
```

### Process outputs
In the previous examples we have connected the INDEX process output to the QUANT process by;

Assigning it to a variable `index_ch = INDEX( transcriptome_ch )` and passing it to the QUANT process as an argument.
Calling the process as an argument within the QUANT process, QUANT( INDEX( transcriptome_ch ), read_pairs_ch )
A process’s output channel can also be accessed calling the process and then using the out attribute for the respective process object.

For example:
```
[..truncated..]

workflow {
    transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
    read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
    
    //call INDEX process
    INDEX(transcriptome_ch)

    // INDEX process output accessed using the `out` attribute
    QUANT(INDEX.out,read_pairs_ch)
    QUANT.out.view()
}
```
When a process defines two or more output channels, each of them can be accessed using the list element operator e.g. `out[0]`, `out[1]`, or using named outputs.

### Process named output
It can be useful to name the output of a process, especially if there are multiple outputs.

The process output definition allows the use of the `emit`: option to define a named identifier that can be used to reference the channel in the external scope.

The scope is the part of the Nextflow script where a named variable is accessible.

For example, in the script below we name the output from the INDEX process as salmon_index using the emit: option. We can then reference the output as INDEX.out.salmon_index in the workflow scope.

```
//workflow_02.nf
nextflow.enable.dsl=2

process INDEX {

  input:
  path transcriptome

  output:
  path 'index', emit: salmon_index

  script:
  """
  salmon index -t $transcriptome -i index
  """
}

process QUANT {
   input:
     each  path(index)
     tuple(val(pair_id), path(reads))
   output:
     path pair_id
   script:
     """
     salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
     """
}

workflow {
  transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
  read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
  
  //call INDEX process
  INDEX(transcriptome_ch)
  
  //access INDEX object named output
  QUANT(INDEX.out.salmon_index,read_pairs_ch).view()
}

```

### Accessing script parameters
A workflow component can access any variable and parameter defined in the outer scope:

For example:
```
//workflow_03.nf
[..truncated..]

params.transcriptome = 'data/yeast/transcriptome/*.fa.gz'
params.reads = 'data/yeast/reads/ref1*_{1,2}.fq.gz'

workflow {
  transcriptome_ch = channel.fromPath(params.transcriptome)
  read_pairs_ch = channel.fromFilePairs(params.reads)
  
  INDEX(transcriptome_ch)
  QUANT(INDEX.out.salmon_index,read_pairs_ch).view()
}
```
In this example `params.transcriptome` and `params.reads` can be accessed inside the workflow scope.

## Key points

- A Nextflow workflow is defined by invoking processes inside the workflow scope.

- A process is invoked like a function inside the workflow scope passing any required input parameters as arguments. e.g. INDEX(transcriptome_ch).

- Process outputs can be accessed using the out attribute for the respective process. Multiple outputs from a single process can be accessed using the `[]` or , if specified , the output name.

