## Reports and visualisations
Nextflow has a number of reports and visualisations it can provide you automatically for any given pipeline you have. Let’s have Nextflow create a few of them for us by executing the following command:
```
nextflow run main.nf -with-report -with-timeline -with-dag dag.png
```
After successful executing, you will find three more files in your current directory: `report.html`, `timeline.html` and `dag.png`. The first file contains a workflow report, which includes various information regarding execution such as runtime, resource usage and details about the different processes. The second file contains a timeline for how long each individual process took to execute, while the last contains a visualisation of the workflow itself.

