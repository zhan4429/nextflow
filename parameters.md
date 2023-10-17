## Pipeline parameters
The Nextflow `wc.nf` script defines a pipeline parameter `params.input`. Pipeline parameters enable you to change the input to the workflow at runtime, via the command line or a configuration file, so they are not hard-coded into the script.

Pipeline parameters are declared in the workflow by prepending the prefix params, separated by the dot character, to a variable name e.g., params.input.

Their value can be specified on the command line by prefixing the parameter name with a double dash character, e.g., `--input`.

In the script wc.nf the pipeline parameter `params.input` was specified with a value of "data/yeast/reads/ref1_1.fq.gz".

To process a different file, e.g. data/yeast/reads/ref2_2.fq.gz, in the wc.nf script we would run:

```
nextflow run wc.nf --input 'data/yeast/reads/ref2_2.fq.gz'

N E X T F L O W  ~  version 21.04.0
Launching `wc.nf` [gigantic_woese] - revision: 8acb5cb9b0
executor >  local (1)
[26/3cf986] process > NUM_LINES (1) [100%] 1 of 1 ✔
ref2_2.fq.gz 81720
```

We can also use wild cards to specify multiple input files (This will be covered in the channels episode). In the example below we use the `*` to match any sequence of characters between ref2_ and .fq.gz. Note: If you use wild card characters on the command line you must enclose the value in quotes.

```
$ nextflow run wc.nf --input 'data/yeast/reads/ref2_*.fq.gz'
```

This runs the process NUM_LINES twice, once for each file it matches.

```
N E X T F L O W  ~  version 21.04.0
Launching `wc.nf` [tender_lumiere] - revision: 8acb5cb9b0
executor >  local (2)
[cc/b6f793] process > NUM_LINES (1) [100%] 2 of 2 ✔
ref2_2.fq.gz 81720

ref2_1.fq.gz 81720
```

### Adding a parameter to a script
To add a pipeline parameter to a script prepend the prefix params, separated by a dot character ., to a variable name e.g., `params.input`.

Let’s make a copy of the wc.nf script as wc-params.nf and add a new input parameter.
```
$ cp wc.nf wc-params.nf
```

To add a parameter sleep with the default value 2 to wc-params.nf we add the line:

```
params.sleep = 2
```

Note: You should always add a sensible default value to the pipeline parameter. We can use this parameter to add another step to our NUM_LINES process.

```
script:
 """
 sleep ${params.sleep}
 printf '${read} '
 gunzip -c ${read} | wc -l
 """
```

This step, sleep ${params.sleep}, will add a delay for the amount of time specified in the params.sleep variable, by default 2 seconds. To access the value inside the script block we use {variable_name} syntax e.g. ${params.sleep}.

We can now change the sleep parameter from the command line, For Example:

```
nextflow run wc-params.nf --sleep 10
```

### Parameter File
If we have many parameters to pass to a script it is best to create a parameters file. Parameters are stored in JSON or YAML format. JSON and YAML are data serialization languages, that are a way of storing data objects and structures, such as the params object in a file.

The `-params-file` option is used to pass the parameters file to the script.

For example the file wc-params.json contains the parameters sleep and input in JSON format.

```
{
  "sleep": 5,
  "input": "data/yeast/reads/etoh60_1*.fq.gz"
}
```

To run the wc-params.nf script using these parameters we add the option `-params-file` and pass the file wc-params.json:

```
$ nextflow run wc-params.nf -params-file wc-params.json
N E X T F L O W  ~  version 21.04.0
Launching `wc-params.nf` [nostalgic_northcutt] - revision: 2f86c9ac7e
executor >  local (2)
[b4/747eaa] process > NUM_LINES (1) [100%] 2 of 2 ✔
etoh60_1_2.fq.gz 87348

etoh60_1_1.fq.gz 87348
```

### Key points
- Pipeline parameters are specified by prepending the prefix `params` to a variable name, separated by dot character.

- To specify a pipeline parameter on the command line for a Nextflow run use `--variable_name` syntax.

- You can add parameters to a JSON or YAML formatted file and pass them to the script using option `-params-file`.