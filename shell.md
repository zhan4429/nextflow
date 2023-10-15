### Shell
The shell block is a string expression that defines the script that is executed by the process. It is an alternative to the Script definition with one important difference: it uses the exclamation mark `!` character, instead of the usual dollar `$` character, to denote Nextflow variables.

This way, it is possible to use both Nextflow and Bash variables in the same script without having to escape the latter, which makes process scripts easier to read and maintain. For example:
```
process myTask {
    input:
    val str

    shell:
    '''
    echo "User $USER says !{str}"
    '''
}

workflow {
    Channel.of('Hello', 'Hola', 'Bonjour') | myTask
}
```
In the above example, $USER is treated as a Bash variable, while !{str} is treated as a Nextflow variable.

Shell script definitions require the use of single-quote `'` delimited strings. When using double-quote `"` delimited strings, dollar variables are interpreted as Nextflow variables as usual. 

Variables prefixed with `!` must always be enclosed in curly brackets, i.e. `!{str}` is a valid variable whereas `!str` is ignored.

Shell scripts support the use of the Template mechanism. The same rules are applied to the variables defined in the script template.