In Nextflow, a process is the basic processing primitive to execute a user script.

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

#### Input repeaters
We saw previously that by default the number of times a process runs is defined by the queue channel with the fewest items. However, the each qualifier allows you to repeat the execution of a process for each item in a list or a queue channel, every time new data is received.

For example if we can fix the previous example by using the input qualifer each for the letters queue channel:

```
//process_repeat.nf
nextflow.enable.dsl=2

process COMBINE {
  input:
  val x
  each y

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
$ nextflow run process_repeat.nf -process.echo
The process will run eight times.

2 and d
1 and a
1 and c
2 and b
2 and c
1 and d
1 and b
2 and a
```

Extend the script `process_exercise_repeat.nf` by adding more values to the kmer queue channel e.g. (21, 27, 31) and running the process for each value.

```
//process_exercise_repeat.nf
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
How many times does this process run?
```
//process_exercise_repeat_answer.nf
nextflow.enable.dsl=2

process COMBINE {
  input:
  path transcriptome
  each kmer
 
  script:
  """
  salmon index -t $transcriptome -i index -k $kmer
  """
}

transcriptome_ch = channel.fromPath('data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz', checkIfExists: true)
kmer_ch = channel.of(21, 27, 31)

workflow {
  COMBINE(transcriptome_ch, kmer_ch)
}
```
```
nextflow run process_exercise_repeat.nf -process.echo
```
This process runs three times.

## Key points
- A Nextflow process is an independent step in a workflow

- Processes contain up to five definition blocks including: `directives`, `inputs`, `outputs`, `when clause` and finally a `script` block.

- The script block contains the commands you would like to run.

- A process should have a script but the other four blocks are optional

- Inputs are defined in the input block with a type qualifier and a name.


## Output
We have seen how to input data into a process; now we will see how to output files and values from a process.

The output declaration block allows us to define the channels used by the process to send out the files and values produced.

An output block is not required, but if it is present it can contain one or more output declarations.

The output block follows the syntax shown below:

```
output:
  <output qualifier> <output name>
  <output qualifier> <output name>
  ...
```

### Output values
Like the input, the type of output data is defined using type qualifiers.

The `val` qualifier allows us to output a value defined in the script.

Because Nextflow processes can only communicate through channels, **if we want to share a value input into one process as input to another process we would need to define that value in the output declaration block** as shown in the following example:

```
//process_output_value.nf
nextflow.enable.dsl=2

process METHOD {
  input:
  val x

  output:
  val x

  script:
  """
  echo $x > method.txt
  """
}
// Both 'Channel' and 'channel' keywords work to generate channels.
// However, it is a good practice to be consistent through the whole pipeline development
methods_ch = channel.of('salmon', 'kallisto')

workflow {
  METHOD(methods_ch)
  // use the view operator to display contents of the channel
  METHOD.out.view({ "Received: $it" })
}

```
#### Output
```
[d8/fd3f9d] process > METHOD (2) [100%] 2 of 2 ✔

Received: salmon
Received: kallisto
```

#### Output files
If we want to capture a file instead of a value as output we can use the `path` qualifier that can capture one or more files produced by the process, over the specified channel.

```
//process_output_file.nf
nextflow.enable.dsl=2

methods_ch = channel.of('salmon', 'kallisto')

process METHOD {
  input:
  val x

  output:
  path 'method.txt'

  """
  echo $x > method.txt
  """
}

workflow {
  METHOD(methods_ch)
  // use the view operator to display contents of the channel
  METHOD.out.view({ "Received: $it" })
}
```
##### output
```
executor >  local (2)
[13/85ecb6] process > METHOD (1) [100%] 2 of 2 ✔

Received: /Users/ggrimes2/Downloads/nf-training/work/ea/f8998abfdcbfc6038757a60806c00a/method.txt
Received: /Users/ggrimes2/Downloads/nf-training/work/13/85ecb69a9bc85370e749254646984d/method.txt
```

In the above example the process METHOD creates a file named method.txt in the work directory containing the method name.

Since a file parameter using the same name, method.txt, is declared in the output block, when the task is completed that file is sent over the output channel.

A downstream operator, such as `.view` or a process declaring the same channel as input will be able to receive it.

#### Multiple output files
When an output file name contains a `*` or `?` metacharacter it is interpreted as a pattern match. This allows us to capture multiple files into a list and output them as a one item channel.

For example, here we will capture the files fastqc.html and directory fastqc.zip produced as results from FastQC in the output channel.
```
//process_output_multiple.nf
nextflow.enable.dsl=2

process FASTQC {
  input:
  path read

  output:
  path "fqc_res/*"

  script:
  """
  mkdir fqc_res
  fastqc $read -o fqc_res
  """
}

read_ch = channel.fromPath("data/yeast/reads/ref1*.fq.gz")

workflow {
  FASTQC(read_ch)
  FASTQC.out.view()
}
```

```
$ nextflow run process_output_multiple.nf
[3e/86de98] process > FASTQC (2) [100%] 2 of 2 ✔
[work/64/9de164199568800a4609c3af78cf71/fqc_res/ref1_2_fastqc.html, work/64/9de164199568800a4609c3af78cf71/fqc_res/ref1_2_fastqc.zip]
Analysis complete for ref1_2.fq.gz

[work/3e/86de9868ecf321702f2df0b8ccbfd3/fqc_res/ref1_1_fastqc.html, work/3e/86de9868ecf321702f2df0b8ccbfd3/fqc_res/ref1_1_fastqc.zip]
Analysis complete for ref1_1.fq.gz
```



Note: There are some caveats on glob pattern behaviour:
- Input files are not included in the list of possible matches.
- Glob pattern matches against both files and directories path.
- When a two stars pattern ** is used to recurse through subdirectories, only file paths are matched i.e. directories are not included in the result list.

#### Grouped inputs and outputs
So far we have seen how to declare multiple input and output channels, but each channel was handling only one value at time. However Nextflow can handle groups of values using the tuple qualifiers.

In tuples the first item is the grouping key and the second item is the list.
```
[group_key,[file1,file2,...]]
```

When using channel containing a tuple, such a one created with .filesFromPairs factory method, the corresponding input declaration must be declared with a tuple qualifier, followed by definition of each item in the tuple.
```
//process_tuple_input.nf
nextflow.enable.dsl=2

process TUPLEINPUT{
  input:
  tuple val(sample_id), path(reads)
  
  script:
  """
  echo $sample_id
  echo $reads
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  TUPLEINPUT(reads_ch)
}

```

##### outputs
```
ref1
ref1_1.fq.gz ref1_2.fq.gz
```
In the same manner an output channel containing tuple of values can be declared using the tuple qualifier following by the definition of each tuple element in the tuple.

In the code snippet below the output channel would contain a tuple with the grouping key value as the Nextflow variable sample_id and a list containing the files matching the following pattern "*FP*.fq.gz".
```
output:
  tuple val(sample_id), path("*FP*.fq.gz")
```
An example can be seen in this script below.
```
//process_tuple_io_fastp.nf
nextflow.enable.dsl=2

process FASTP {
  input:
  tuple val(sample_id), path(reads)
  
  output:
  tuple val(sample_id), path("*FP*.fq.gz")
  
  script:
  """
  fastp \
   -i ${reads[0]} \
   -I ${reads[1]} \
   -o ${sample_id}_FP_R1.fq.gz \
   -O ${sample_id}_FP_R2.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  FASTP(reads_ch)
  FASTP.out.view()
}
```

```
nextflow run process_tuple_io_fastp.nf
```
The output is now a tuple containing the sample id and the two processed fastq files.
```
[ref1, [work/9b/fcca75db83718a15f7a95caabfbc15/ref1_FP_R1.fq.gz, work/9b/fcca75db83718a15f7a95caabfbc15/ref1_FP_R2.fq.gz]]
```
#### Conditional execution of a process
The when declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value; `true` or `false`.

It is useful to enable/disable the process execution depending on the state of various inputs and parameters.

In the example below the process CONDITIONAL will only execute when the value of the chr variable is less than or equal to 5:

```
//process_when.nf
nextflow.enable.dsl=2

process CONDITIONAL {
  input:
  val chr

  when:
  chr <= 5

  script:
  """
  echo $chr
  """
}

chr_ch = channel.of(1..22)

workflow {
  CONDITIONAL(chr_ch)
}
```
```
4

5

2

3

1
```

## Directives
Directive declarations allow the definition of optional settings, like the number of cpus and amount of memory, that affect the execution of the current process without affecting the task itself.

**They must be entered at the top of the process body, before any other declaration blocks (i.e. input, output, etc)**.

**Note: You do not use = when assigning a value to a directive.**

Directives are commonly used to define the amount of computing resources to be used or extra information for configuration or logging purpose.

For example:

```
//process_directive.nf
nextflow.enable.dsl=2

process PRINTCHR {
  tag "tagging with chr$chr"
  cpus 1
  echo true

  input:
  val chr

  script:
  """
  echo processing chromosome: $chr
  echo number of cpus $task.cpus
  """
}

chr_ch = channel.of(1..22, 'X', 'Y')

workflow {
  PRINTCHR(chr_ch)
}

```

##### Outout
```
processing chromosome: 1
number of cpus 1

processing chromosome: 2
number of cpus 1

processing chromosome: 6
number of cpus 1
[..truncated..]
```
The above process uses the three directives, `tag`, `cpus` and `echo`.

The `tag` directive to allow you to give a custom tag to each process execution. This tag makes it easier to identify a particular task (executed instance of a process) in a log file or in the execution report.

The second directive `cpus` allows you to define the number of CPUs required for each task.

The third directive `echo true` prints the stdout to the terminal.

We use the Nextflow task.cpus variable to capture the number of cpus assigned to a task. This is frequently used to specify the number of threads in a multi-threaded command in the script block.

Another commonly used directive is memory specification: memory.

A complete list of directives is available at this link.

Another example:
```
//process_directives_answer.nf
nextflow.enable.dsl=2

process FASTQC {
  tag "$sample_id"
  cpus 2
  
  input:
  tuple val(sample_id), path(reads)
  
  output:
  tuple val(sample_id), path("fastqc_out")
 
  script:
  """
  mkdir fastqc_out
  fastqc $reads -o fastqc_out -t 1
  """
}

read_pairs_ch = Channel.fromFilePairs('data/yeast/reads/ref*_{1,2}.fq.gz')

workflow {
  FASTQC(read_pairs_ch)
  FASTQC.out.view()
}
N E X T F L O W  ~  version 21.04.0
Launching `process_exercise_directives.nf` [sad_rosalind] - revision: 2ccbfa4937
executor >  local (3)
[90/de1125] process > FASTQC (ref1) [100%] 3 of 3 ✔
[ref2, work/ea/9e6a341b88caf8879e8d18b77049c8/fastqc_out]
[ref3, work/94/d059b816a9ec3d868f2924c26813e7/fastqc_out]
[ref1, work/90/de11251d362f494d6650789d9f8c1d/fastqc_out]
```
### Organising outputs
#### PublishDir directive
Nextflow manages intermediate results from the pipeline’s expected outputs independently.

Files created by a process are stored in a task specific working directory which is considered as temporary. Normally this is under the `work` directory, which can be deleted upon completion.

The files you want the workflow to return as results need to be defined in the process output block and then the output directory specified using the directive `publishDir`. More information here.

**Note: A common mistake is to specify an output directory in the publishDir directive while forgetting to specify the files you want to include in the output block**.

```
publishDir <directory>, parameter: value, parameter2: value ...
```

For example if we want to capture the results of the QUANT process in a `results/quant` output directory we need to define the files in the output and specify the location of the results directory in the `publishDir` directive:
```
//process_publishDir.nf
nextflow.enable.dsl=2

process QUANT {
  publishDir "results/quant"
  
  input:
  tuple val(sample_id), path(reads)
  each index
  
  output:
  tuple val(sample_id), path("${sample_id}_salmon_output")
  
  script:
  """
  salmon quant -i $index \
   -l A \
   -1 ${reads[0]} \
   -2 ${reads[1]} \
   -o ${sample_id}_salmon_output
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

workflow {
  QUANT(reads_ch, index_ch)
  QUANT.out.view()
}
```
```
$ nextflow run process_publishDir.nf
N E X T F L O W  ~  version 21.04.0
Launching `process_publishDir.nf` [friendly_pauling] - revision: 9b5c315893
executor >  local (1)

[48/f97234] process > QUANT (1) [100%] 1 of 1 ✔
[ref1, work/48/f97234d7185cbfbd86e2f11c1afab5/ref1_salmon_output]
We can use the UNIX command tree to examine the contents of the results directory.
```
```
tree results
results/
└── quant
    └── ref1_salmon_output -> work/48/f97234d7185cbfbd86e2f11c1afab5/ref1_salmon_output
```
In the above example, the publishDir "results/quant", creates **a symbolic link** -> to the output files specified by the process salmon_quant to the directory path results/quant.

##### publishDir
The publishDir output is relative to the path the pipeline run has been launched. Hence, it is a good practice to use implicit variables (https://www.nextflow.io/docs/latest/script.html?highlight=projectdir#script-implicit-variables) like projectDir to specify publishDir value.

#### publishDir parameters
The publishDir directive can take optional parameters, for example the mode parameter can take the value "copy" to specify that you wish to copy the file to output directory rather than just a symbolic link to the files in the working directory. Since the working directory is generally deleted on completion of a pipeline, it is safest to use mode: "`copy`" for results files. The `default mode (symlink)` is helpful for checking intermediate files which are not needed in the long term.

```
publishDir "results/quant", mode: "copy"
```
Full list here (https://www.nextflow.io/docs/latest/process.html#publishdir).

### Manage semantic sub-directories
You can use more than one publishDir to keep different outputs in separate directories. To specify which files to put in which output directory use the parameter pattern with the a glob pattern that selects which files to publish from the overall set of output files.

In the example below we will create an output folder structure in the directory results, which contains a separate sub-directory for bam files, pattern:"*.bam" , and a salmon output directory, ${sample_id}_salmon_output". Remember, we need to specify the files we want to copy as outputs.

```
//process_publishDir_semantic.nf
nextflow.enable.dsl=2

process QUANT {
  publishDir "results/bams", pattern: "*.bam", mode: "copy"
  publishDir "results/quant", pattern: "${sample_id}_salmon_output", mode: "copy"

  input:
  tuple val(sample_id), path(reads)
  path index
  
  output:
  tuple val(sample_id), path("${sample_id}.bam")
  path "${sample_id}_salmon_output"
  
  script:
  """
  salmon quant -i $index \
   -l A \
   -1 ${reads[0]} \
   -2 ${reads[1]} \
   -o ${sample_id}_salmon_output \
   --writeMappings | samtools sort | samtools view -bS -o ${sample_id}.bam
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')
index_ch = Channel.fromPath('data/yeast/salmon_index')

workflow {
  QUANT(reads_ch, index_ch)
}
```

```
$ nextflow run process_publishDir_semantic.nf
N E X T F L O W  ~  version 21.04.0
Launching `process_publishDir_semantic.nf` [golden_poisson] - revision: 421a604840

executor >  local (1)
[be/950786] process > QUANT (1) [100%] 1 of 1 ✔
```

We can now use the tree command to examine the results directory.

```
$ tree results
results/
├── bams
│   └── ref1.bam
└── quant
    └── ref1_salmon_output
        ├── aux_info
        │   ├── ambig_info.tsv
        │   ├── expected_bias.gz
        │   ├── fld.gz
        │   ├── meta_info.json
        │   ├── observed_bias.gz
        │   └── observed_bias_3p.gz
        ├── cmd_info.json
        ├── libParams
        │   └── flenDist.txt
        ├── lib_format_counts.json
        ├── logs
        │   └── salmon_quant.log
        └── quant.sf

6 directories, 12 files
```

## Key Points
- Outputs to a process are defined using the output blocks.

- You can group input and output data from a process using the tuple qualifer.

- The execution of a process can be controlled using the when declaration and conditional statements.

- Files produced within a process and defined as output can be saved to a directory using the `publishDir` directive.



