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



### Script parameters

The command in the script block can be defined dynamically using Nextflow variables e.g. `${projectDir}`. To reference a variable in the script block you can use the `$` in front of the Nextflow variable name, and additionally you can add `{}` around the variable name e.g. `${projectDir}`.

#### Variable substitutions

Similar to bash scripting Nextflow uses the “`$`” character to introduce variable substitutions. The variable name to be expanded may be enclosed in braces `{variable_name}`, which are optional but serve to protect the variable to be expanded from characters immediately following it which could be interpreted as part of the name. It is a good rule of thumb to always use the `{}` syntax.
In the example below the variable kmer is set to the value 31 at the top of the Nextflow script. The variable is referenced using the $kmer syntax within the multi-line string statement in the script block. A Nextflow variable can be used multiple times in the script block.
```
//process_script.nf
nextflow.enable.dsl=2

kmer = 31

process INDEX {

  script:
  """
  salmon index \
  -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  -i index \
  --kmer $kmer
  echo "kmer size is $kmer"
  """
}

workflow {
  INDEX()
}
```

In most cases we do not want to hard code parameter values. We saw in the parameter episode the use of a special Nextflow variable params that can be used to assign values from the command line. You would do this by adding a key name to the params variable and specifying a value, like `params.keyname = value`.

In the example below we define the variable params.kmer with a default value of 31 in the Nextflow script.
```
//process_script_params.nf
nextflow.enable.dsl=2

params.kmer = 31

process INDEX {

  script:
  """
  salmon index \
  -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  -i index \
  --kmer $params.kmer
  echo "kmer size is $params.kmer"
  """
}

workflow {
  INDEX()
}
```

Remember, we can change the default value of kmer to 11 by running the Nextflow script using the command below. Note: parameters to the workflow have two hyphens `--`.

```
nextflow run process_script_params.nf --kmer 11
```

#### Bash variables

Nextflow uses the same Bash syntax for variable substitutions, `$variable`, in strings. However, Bash variables need to be escaped using `\` character in front of `\$variable` name.

In the example below we will set the bash variable KMERSIZE to the value of $params.kmer, and then use KMERSIZE in our script block.

```
//process_escape_bash.nf
nextflow.enable.dsl=2

process INDEX {

  script:
  """
  #set bash variable KMERSIZE
  KMERSIZE=$params.kmer
  salmon index -t $projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmer \$KMERSIZE
  echo "kmer size is $params.kmer"
  """
}

params.kmer = 31

workflow {
  INDEX()
}
```


#### Shell

Another alternative is to use a `shell` block definition instead of script. When using the shell statement Bash variables are referenced in the normal way `$my_bash_variable`; However, the shell statement uses a different syntax for Nextflow variable substitutions: `!{nextflow_variable}`, which is needed to use both Nextflow and Bash variables in the same script.

For example in the script below that uses the shell statement we reference the Nextflow variables as `!{projectDir}` and `!{params.kmer}`, and the Bash variable as `${KMERSIZE}`.

```
//process_shell.nf
nextflow.enable.dsl=2

params.kmer = 31

process INDEX {

  shell:
  '''
  #set bash variable KMERSIZE
  KMERSIZE=!{params.kmer}
  salmon index -t !{projectDir}/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -i index --kmer ${KMERSIZE}
  echo "kmer size is !{params.kmer}"
  '''
}

workflow {
  INDEX()
}
```


#### Conditional script execution

Sometimes you want to change how a process is run depending on some condition. In Nextflow scripts we can use conditional statements such as the `if` statement or any other expression evaluating to boolean value true or false.

#### If statement

The if statement uses the same syntax common to other programming languages such Java, C, JavaScript, etc.

```
if( < boolean expression > ) {
    // true branch
}
else if ( < boolean expression > ) {
    // true branch
}
else {
    // false branch
}
```

For example, the Nextflow script below will use the if statement to change which index is created depending on the Nextflow variable `params.aligner`.

```
//process_conditional.nf
nextflow.enable.dsl=2

params.aligner = 'kallisto'
params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
params.kmer = 31

process INDEX {
  script:
  if( params.aligner == 'kallisto' ) {
    """
    echo indexed using kallisto
    kallisto index -i index -k $params.kmer $params.transcriptome
    """
  }  
  else if( params.aligner == 'salmon' ) {
    """
    echo indexed using salmon
    salmon index -t $params.transcriptome -i index --kmer $params.kmer
    """
  }  
  else {
    """
    echo Unknown aligner $params.aligner
    """
  }  
}

workflow {
  INDEX()
}
```

```
nextflow run process_conditional.nf -process.echo --aligner kallisto
N E X T F L O W  ~  version 21.04.0
Launching `juggle_processes.nf` [cheeky_shirley] - revision: 588f20ae5a
executor >  local (1)
[84/c44f25] process > INDEX [100%] 1 of 1 ✔
indexed using kallisto
```

## Inputs

Processes are isolated from each other but can communicate by sending values and files via Nextflow channels from input and into output blocks.

The input block defines which channels the process is expecting to receive input from. The number of elements in input channels determines the process dependencies and the number of times a process executes.

![Channels-process](images/channel-process.png)

You can only define one input block at a time and it must contain one or more input declarations.

The input block follows the syntax shown below:

```
input:
  <input qualifier> <input name>
```

The input qualifier declares the type of data to be received.

### Input qualifiers

- val: Lets you access the received input value by its name as a variable in the process script.
- env: Lets you use the input value to set an environment variable named as the specified input name.
- path: Lets you handle the received value as a file, staging the file properly in the execution context.
- stdin: Lets you forward the received value to the process stdin special file.
- tuple: Lets you handle a group of input values having one of the above qualifiers.
- each: Lets you execute the process for each entry in the input collection. 