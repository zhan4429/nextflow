## module
Environment Modules is a package manager that allows you to dynamically configure your execution environment and easily switch between multiple versions of the same software tool.

If it is available in your system you can use it with Nextflow in order to configure the processes execution environment in your pipeline.

In a process definition you can use the module directive to load a specific module version to be used in the process execution environment. For example:

```
process basicExample {
  module 'ncbi-blast/2.2.27'

  """
  blastp -query <etc..>
  """
}
```

You can repeat the module directive for each module you need to load. Alternatively multiple modules can be specified in a single module directive by separating all the module names by using a `:` (colon) character as shown below:

```
 process manyModules {

   module 'ncbi-blast/2.2.27:t_coffee/10.0:clustalw/2.1'

   """
   blastp -query <etc..>
   """
}
```