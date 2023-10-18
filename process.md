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


#### Input values
The `val` qualifier allows you to receive value data as input. It can be accessed in the process script by using the specified input name, as shown in the following example:

```
//process_input_value.nf
nextflow.enable.dsl=2

process PRINTCHR {

  input:
  val chr

  script:
  """
  echo processing chromosome $chr
  """
}

chr_ch = Channel.of( 1..22, 'X', 'Y' )

workflow {

  PRINTCHR(chr_ch)
}

```

```
$ nextflow run process_input_value.nf -process.echo
N E X T F L O W  ~  version 21.04.0
Launching `juggle_processes.nf` [wise_kalman] - revision: 7f90e1bfc5
executor >  local (24)
[b1/88df3f] process > PRINTCHR (24) [100%] 24 of 24 ✔
processing chromosome 3
processing chromosome 1
processing chromosome 2
..truncated...
```

In the above example the process is executed 24 times; each time a value is received from the queue channel `chr_ch` it is used to run the process.

##### Channel order
The channel guarantees that items are delivered in the same order as they have been sent, but since the process is executed in a parallel manner, there is no guarantee that they are processed in the same order as they are received.

#### Input files
When you need to handle files as input you need the `path` qualifier. Using the `path` qualifier means that Nextflow will stage it in the process execution directory, and it can be accessed in the script by using the name specified in the input declaration.

The input file name can be defined dynamically by defining the input name as a Nextflow variable and referenced in the script using the  `$variable_name` syntax.

For example in the script below we assign the variable name read to the input files using the `path` qualifier. The file is referenced using the variable substitution syntax `${read}` in the script block:

```
//process_input_file.nf
nextflow.enable.dsl=2

process NUMLINES {
    input:
    path read

    script:
    """
    printf '${read} '
    gunzip -c ${read} | wc -l
    """
}

reads_ch = Channel.fromPath( 'data/yeast/reads/ref*.fq.gz' )

workflow {
  NUMLINES(reads_ch)
}
```

```
$ nextflow run process_input_file.nf -process.echo
[cd/77af6d] process > NUMLINES (1) [100%] 6 of 6 ✔
ref1_1.fq.gz 58708

ref3_2.fq.gz 52592

ref2_2.fq.gz 81720

ref2_1.fq.gz 81720

ref3_1.fq.gz 52592

ref1_2.fq.gz 58708
```

The input name can also be defined as user specified filename inside `quotes`. For example in the script below the name of the file is specified as `sample.fq.gz` in the input definition and can be referenced by that name in the script block.

```
//process_input_file_02.nf
nextflow.enable.dsl=2

process NUMLINES {
    input:
    path 'sample.fq.gz'

    script:
    """
    printf 'sample.fq.gz '
    gunzip -c sample.fq.gz | wc -l
    """
}

reads_ch = Channel.fromPath( 'data/yeast/reads/ref*.fq.gz' )

workflow {
  NUMLINES(reads_ch)
}
```

```
$ nextflow run process_input_file_02.nf -process.echo
[d2/eb0e9d] process > NUMLINES (1) [100%] 6 of 6 ✔
sample.fq.gz 58708

sample.fq.gz 52592

sample.fq.gz 81720

sample.fq.gz 81720

sample.fq.gz 52592

sample.fq.gz 58708
```

##### File Objects as inputs
When a process declares an input file the corresponding channel elements must be file objects i.e. created with the `path` helper function from the file specific channel factories e.g. `Channel.fromPath` or `Channel.fromFilePairs`.

#### Add input channel 
Add an input channel to the script below that takes the reads channel as input. FastQC is a quality control tool for high throughput sequence data.
```
//process_exercise_input.nf
nextflow.enable.dsl=2

process FASTQC {
   //add input channel
   
   script:
   """
   mkdir fastqc_out
   fastqc -o fastqc_out ${reads}
   ls -1 fastqc_out
   """
}
reads_ch = Channel.fromPath( 'data/yeast/reads/ref1*_{1,2}.fq.gz' )

workflow {
  FASTQC(reads_ch)
}
```
Then run your script using
```
nextflow run process_exercise_input.nf -process.echo
```
Solution
```
//process_exercise_input_answer.nf
nextflow.enable.dsl=2
process FASTQC {
   input:
   path reads

   script:
   """
   mkdir fastqc_out
   fastqc -o fastqc_out ${reads}
   ls -1 fastqc_out
   """
}
reads_ch = Channel.fromPath( 'data/yeast/reads/ref1*_{1,2}.fq.gz' )

workflow {
  FASTQC(reads_ch)
}
```

```
N E X T F L O W  ~  version 21.04.0
Launching `process_exercise_input_answer.nf` [jovial_wescoff] - revision: e3db00a4dc
executor >  local (2)
[d9/559a27] process > FASTQC (2) [100%] 2 of 2 ✔
Analysis complete for ref1_1.fq.gz
ref1_1_fastqc.html
ref1_1_fastqc.zip

Analysis complete for ref1_2.fq.gz
ref1_2_fastqc.html
ref1_2_fastqc.zip

```

#### Combining input channels
A key feature of processes is the ability to handle inputs from multiple channels. However it’s important to understand how the number of items within the multiple channels affect the execution of a process.

Consider the following example:
```
//process_combine.nf
nextflow.enable.dsl=2

process COMBINE {
  input:
  val x
  val y

  script:
  """
  echo $x and $y
  """
}

num_ch = Channel.of(1, 2, 3)
letters_ch = Channel.of('a', 'b', 'c')

workflow {
  COMBINE(num_ch, letters_ch)
}
```
```
$ nextflow run process_combine.nf -process.echo
```
Both channels contain three elements, therefore the process is executed three times, each time with a different pair:

2 and b

1 and a

3 and c

What is happening is that the process waits until it receives an input value from all the queue channels declared as input.

When this condition is verified, it uses up the input values coming from the respective queue channels, runs the task. This logic repeats until one or more queue channels have no more content. The process then stops.

**What happens when not all channels have the same number of elements?**

For example:
```
//process_combine_02.nf
nextflow.enable.dsl=2

process COMBINE {
  input:
  val x
  val y

  script:
  """
  echo $x and $y
  """
}

ch_num = Channel.of(1, 2)
ch_letters = Channel.of('a', 'b', 'c', 'd')

workflow {
  COMBINE(ch_num, ch_letters)
}
```

```
$ nextflow run process_combine_02.nf -process.echo
```
In the above example the process is executed only two times, because when a queue channel has no more data to be processed it stops the process execution.
```
2 and b
1 and a
```

#### Value channels and process termination
Note however that value channels, `Channel.value`, do not affect the process termination.

To better understand this behaviour compare the previous example with the following one:
```
//process_combine_03.nf
nextflow.enable.dsl=2

process COMBINE {
  input:
  val x
  val y

  script:
  """
  echo $x and $y
  """
}
ch_num = Channel.value(1)
ch_letters = Channel.of('a', 'b', 'c')

workflow {
  COMBINE(ch_num, ch_letters)
}
```

```
$ nextflow run process_combine_03.nf -process.echo
```
In this example the process is run three times.

```
1 and b
1 and a
1 and c
```
#### Combining input channels
Write a nextflow script process_exercise_combine.nf that combines two input channels
```
transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz')
kmer_ch = channel.of(21)
```
And include the command below in the script directive

```
  script:
  """
  salmon index -t $transcriptome -i index -k $kmer .
  """
```
Solution
```
// process_exercise_combine_answer.nf
nextflow.enable.dsl=2
process COMBINE {
 input:
 path transcriptome
 val kmer

 script:
 """
 salmon index -t $transcriptome -i index -k $kmer
 """
}

transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz', checkIfExists: true)
kmer_ch = channel.of(21)

workflow {
  COMBINE(transcriptome_ch, kmer_ch)
}
```

