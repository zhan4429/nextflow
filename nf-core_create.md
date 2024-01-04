https://nf-co.re/docs/contributing/tutorials/creating_with_nf_core

# Creating a new pipeline

If you can't find a proper pipeline in community, you could create a pipeline by your
self. You can re-use modules in which calculations steps are defined
by the community. In such way, you can avoid to write a full pipeline from yourself.

The minimal set of files required to have a pipeline is to have locally
**main.nf**, **nextflow.config** and **modules.json** inside your project folder.
You should have also a **modules** directory inside your project::

```
  $ mkdir -p my-new-pipeline/modules
  $ cd my-new-pipeline
  $ touch main.nf nextflow.config modules.json README.md
```

Next you have to edit modules.json in order to have minimal information:
```
{
  "name": "<your pipeline name>",
  "homePage": "<your pipeline repository URL>",
  "repos": { }
}
```
Without this requisites you will not be able to add community modules to your pipelines using `nf-core/tools`.

## Browsing modules list

You can get a list of modules by using nf-core/tools (see here how you can install it):
```
$ nf-core modules list remote
```

### Adding a module to a pipeline

You can download and add a module to your pipeline using nf-core/tools:
```
$ nf-core modules install --dir . fastqc
```

The `--dir .` option is optional, the default installation path is the `CWD` (that need to be your pipeline source directory).

### List all modules in a pipeline
You can have a full list of installed modules using:
```
$ nf-core modules list local
```

### Update a pipeline module

You can update a module simple by calling:

```
$ nf-core modules update fastqc
```
