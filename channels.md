## Channels
Earlier we learnt that channels are the way in which Nextflow sends data around a workflow. Channels connect processes via their inputs and outputs. Channels can store multiple items, such as files (e.g., fastq files) or values. The number of items a channel stores determines how many times a process will run using that channel as input.
**Note: When the process runs using one item from the input channel, we will call that run a task.** 

## Why use Channels?
Channels are how Nextflow handles file management, allowing complex tasks to be split up, run in parallel, and reduces ‘admin’ required to get the right inputs to the right parts of the pipeline.

![Channels](images/channel-files.png)

Channels are asynchronous, which means that outputs from a set of processes will not necessarily be produced in the same order as the corresponding inputs went in. However, the first element into a channel queue is the first out of the queue (First in - First out). This allows processes to run as soon as they receive input from a channel. Channels only send data in one direction, from a producer (a process/operator), to a consumer (another process/operator).

## Channel types
Nextflow distinguishes between two different kinds of channels: `queue channels` and `value channels`.

### Queue channel
Queue channels are a type of channel in which data is consumed (used up) to make input for a process/operator. Queue channels can be created in two ways:

1. As the outputs of a process.
2. Explicitly using channel factory methods such as `Channel.of` or `Channel.fromPath`.

#### DSL1
**In Nextflow DSL1 queue channels can only be used once in a workflow, either connecting workflow input to process input, or process output to input for another process. In DSL2 we can use a queue channel multiple times.**

### Value channels
The second type of Nextflow channel is a value channel. A value channel is bound to a single value. A value channel can be used an unlimited number times since its content is not consumed. This is also useful for processes that need to reuse input from a channel, for example, a reference genome sequence file that is required by multiple steps within a process, or by more than one process.

### Creating Channels using Channel factories
Channel factories are used to explicitly create channels. In programming, factory methods (functions) are a programming design pattern used to create different types of objects (in this case, different types of channels). They are implemented for things that represent more generalised concepts, such as a Channel.

Channel factories are called using the `Channel.<method>` syntax, and return a specific instance of a `Channel`.

### The value Channel factory
The value factory method is used to create a value channel. Values are put inside parentheses () to assign them to a channel.

For example:

ch1 = Channel.value( 'GRCh38' )
ch2 = Channel.value( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] )
ch3 = Channel.value( ['chr1' : 248956422, 'chr2' : 242193529, 'chr3' : 198295559] )
Creates a value channel and binds a string to it.
Creates a value channel and binds a list object to it that will be emitted as a single item.
Creates a value channel and binds a map object to it that will be emitted as a single item.
The value method can only take 1 argument, however, this can be a single list or map containing several elements.

Reminder:

A List object can be defined by placing the values in square brackets [] separated by a comma.
A Map object is similar, but with key:value pairs separated by commas.
myList = [1776, -1, 33, 99, 0, 928734928763]
myMap = [ p1 : "start", q2 : "end" ]