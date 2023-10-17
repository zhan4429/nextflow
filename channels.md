### Channels
Earlier we learnt that channels are the way in which Nextflow sends data around a workflow. Channels connect processes via their inputs and outputs. Channels can store multiple items, such as files (e.g., fastq files) or values. The number of items a channel stores determines how many times a process will run using that channel as input.
**Note: When the process runs using one item from the input channel, we will call that run a task.** 

### Why use Channels?
Channels are how Nextflow handles file management, allowing complex tasks to be split up, run in parallel, and reduces ‘admin’ required to get the right inputs to the right parts of the pipeline.

![Channels](images/channel-files.png)
