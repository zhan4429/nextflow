n Nextflow, a process is the basic processing primitive to execute a user script.

The process definition starts with the keyword `process`, followed by process name and finally the process body delimited by curly brackets. The process body must contain a string which represents the command or, more generally, a script that is executed by it. A basic process looks like the following example:

```
process doMoreThings {
  shell '/bin/bash', '-euo', 'pipefail'

  """
  blastp -db $db -query query.fa -outfmt 6 > blast_result
  cat blast_result | head -n 10 | cut -f 2 > top_hits
  blastdbcmd -db $db -entry_batch top_hits > sequences
  """
}
```

Since Nextflow uses the same Bash syntax for variable substitutions in strings, you must manage them carefully depending on whether you want to evaluate a Nextflow variable or a Bash variable.

When you need to access a system environment variable in your script, you have two options.

If you don’t need to access any Nextflow variables, you can define your script block with single-quotes:
```
process printPath {
  '''
  echo The path is: $PATH
  '''
}
```
Otherwise, you can define your script with double-quotes and escape the system environment variables by prefixing them with a back-slash `\` character, as shown in the following example:

```
process doOtherThings {
  """
  blastp -db \$DB -query query.fa -outfmt 6 > blast_result
  cat blast_result | head -n $MAX | cut -f 2 > top_hits
  blastdbcmd -db \$DB -entry_batch top_hits > sequences
  """
}
```
In this example, `$MAX` is a Nextflow variable that must be defined elsewhere in the pipeline script. Nextflow replaces it with the actual value before executing the script. Meanwhile, `$DB` is a Bash variable that must exist in the execution environment, and Bash will replace it with the actual value during execution.


## Process
A process is the way Nextflow executes commands you would run on the command line or custom scripts.

A process can be thought of as a particular step in a workflow, e.g. an alignment step in RNA-seq analysis. Processes are independent of each other (don’t require any another process to execute) and can not communicate/write to each other. Data is passed between processes via input and output Channels.

For example, below is the command line you would run to create a index for the yeast transcriptome to be used with the salmon aligner:

### Process definition

The process definition starts with keyword `process`, followed by process name, in this case INDEX, and finally the process body delimited by curly brackets `{}`. The process body must contain a string which represents the command or, more generally, a script that is executed by it.

```
process INDEX {
  script:
  "salmon index -t ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i ${projectDir}/data/yeast/salmon_index --kmerLen 31"
}
```

To add the process to a workflow add a workflow block, and call the process like a function. 
Note: As we are using DSL2 we need to include `nextflow.enable.dsl=2` in the script.

```
//process_index.nf
nextflow.enable.dsl=2

process INDEX {
  script:
  "salmon index -t ${projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i data/yeast/salmon_index --kmerLen 31"
}

workflow {
  //process is called like a function in the workflow block
  INDEX()
}
```
### Definition blocks
The previous example was a simple process with no defined inputs and outputs that ran only once. To control inputs, outputs and how a command is executed a process may contain five definition blocks:

1. directives - 0, 1, or more: allow the definition of optional settings that affect the execution of the current process e.g. the number of cpus a task uses and the amount of memory allocated.
2. inputs - 0, 1, or more: Define the input dependencies, usually channels, which determines the number of times a process is executed.
3. outputs - 0, 1, or more: Defines the output channels used by the process to send results/data produced by the process.
4. when clause - optional: Allows you to define a condition that must be verified in order to execute the process.
5. script block - required: A statement within quotes that defines the commands that are executed by the process to carry out its task.
The syntax is defined as follows:

```
process < NAME > {
  [ directives ]        
  input:                
  < process inputs >
  output:               
  < process outputs >
  when:                 
  < condition >
  [script|shell|exec]:  
  < user script to be executed >
}
```

### Script

At minimum a process block must contain a script block.

The script block is a String “statement” that defines the command that is executed by the process to carry out its task. These are normally the commands you would run on a terminal.

**A process contains only one script block, and it must be the last statement when the process contains input and output declarations.**

The script block can be a simple one line string in quotes e.g.
```
nextflow.enable.dsl=2

process PROCESSBAM {
    script:
    "samtools sort -o ref1.sorted.bam ${projectDir}/data/yeast/bams/ref1.bam"
}

workflow {
  PROCESSBAM()
}
```

Or, for commands that span multiple lines you can encase the command in triple quotes `"""`.

For example:

```
//process_multi_line.nf
nextflow.enable.dsl=2

process PROCESSBAM {
    script:
    """
    samtools sort -o ref1.sorted.bam ${projectDir}/data/yeast/bams/ref1.bam
    samtools index ref1.sorted.bam
    samtools flagstat ref1.sorted.bam
    """
}

workflow {
  PROCESSBAM()
}
```
By default the process command is interpreted as a Bash script. However any other scripting language can be used just simply starting the script with the corresponding Shebang declaration. For example:
```
//process_python.nf
nextflow.enable.dsl=2

process PYSTUFF {
  script:
  """
  #!/usr/bin/env python
  import gzip

  reads = 0
  bases = 0

  with gzip.open('${projectDir}/data/yeast/reads/ref1_1.fq.gz', 'rb') as read:
      for id in read:
          seq = next(read)
          reads += 1
          bases += len(seq.strip())
          next(read)
          next(read)

  print("reads", reads)
  print("bases", bases)
  """
}

workflow {
  PYSTUFF()
}

```
//process_rscript.nf
nextflow.enable.dsl=2

process RSTUFF {
  script:
  """
  #!/usr/bin/env Rscript
  library("ShortRead")
  countFastq(dirPath="data/yeast/reads/ref1_1.fq.gz")
  """
}

workflow {
  RSTUFF()
}
```
This allows the use of a different programming languages which may better fit a particular job. However, for large chunks of code it is suggested to save them into separate files and invoke them from the process script.

```
nextflow.enable.dsl=2

process PYSTUFF {

  script:
  """
  myscript.py
  """
}

workflow {
  PYSTUFF()
}
```

Scripts such as the one in the example above, myscript.py, can be stored in a bin folder at the same directory level as the Nextflow workflow script that invokes them, and given execute permission. Nextflow will automatically add this folder to the PATH environment variable. To invoke the script in a Nextflow process, simply use its filename on its own rather than invoking the interpreter e.g. `myscript.py` instead of `python myscript.py`. 



