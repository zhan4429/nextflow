## Operators
In the Channels episode we learnt how to create Nextflow channels to enable us to pass data and values around our workflow. If we want to modify the contents or behaviour of a channel, Nextflow provides methods called operators. We have previously used the view operator to view the contents of a channel. There are many more operator methods that can be applied to Nextflow channels that can be usefully separated into several groups:

- Filtering operators: reduce the number of elements in a channel.
- Transforming operators: transform the value/data in a channel.
- Splitting operators: split items in a channel into smaller chunks.
- Combining operators: join channel together.
- Forking operators: split a single channel into multiple channels.
- Maths operators: apply simple math function on channels.
- Other: such as the view operator.
In this episode you will see examples, and get to use different types of operators.

## Using Operators
To use an operator, the syntax is the channel name, followed by a dot `.` , followed by the operator name and brackets `()`.

```
channel_obj.<operator>()
```

### view

The `view` operator prints the items emitted by a channel to the console appending a new line character to each item in the channel.
```
ch = channel.of('1', '2', '3')
ch.view()
```
We can also chain together the channel factory method `.of` and the operator `.view()` using the dot notation.

```
ch = channel.of('1', '2', '3').view()
```

To make code more readable we can split the operators over several lines. **The blank space between the operators is ignored and is solely for readability**.

```
ch = channel
      .of('1', '2', '3')
      .view()
```
prints:
```
1
2
3
```

#### Closures

An optional closure `{}` parameter can be specified to customise how items are printed.

Briefly, a closure is a block of code that can be passed as an argument to a function. In this way you can define a chunk of code and then pass it around as if it were a string or an integer. By default the parameters for a closure are specified with the groovy keyword `$it (‘it’ is for ‘item’)`.

For example here we use the the view operator and apply a closure to it, to add a chr prefix to each element of the channel using string interpolation.
```
ch = channel
  .of('1', '2', '3')
  .view({ "chr$it" })
```
It prints:
```
chr1
chr2
chr3
```
**Note: the `view()` operator doesn’t change the contents of the channel object**.
```
ch = channel
  .of('1', '2', '3')
  .view({ "chr$it" })

ch.view()  
```
```
chr1
chr2
chr3
1
2
3
```
### Filtering operators

We can reduce the number of items in a channel by using filtering operators.

The filter operator allows you to get only the items emitted by a channel that satisfy a condition and discard all the others. The filtering condition can be specified by using either:

- a regular expression
- a literal value
- a data type qualifier, e.g. Number (any integer,float …), String, Boolean
- or any boolean statement.
- Data type qualifier

Here we use the filter operator on the chr_ch channel specifying the data type qualifier Number so that only numeric items are returned. The Number data type includes both integers and floating point numbers. We will then use the view operator to print the contents.
```
chr_ch = channel.of( 1..22, 'X', 'Y' )
autosomes_ch =chr_ch.filter( Number )
autosomes_ch.view()
```
To simplify the code we can chain multiple operators together, such as filter and view using a `.`.

The previous example could be rewritten like: The blank space between the operators is ignored and is used for readability.
```
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter( Number )
  .view()
```
Output
```
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
```

#### Regular expression

To filter by a regular expression you have to do is to put `~` right in front of the string literal regular expression (e.g. ~"(^[Nn]extflow)" or use slashy strings which replace the quotes with /. ~/^[Nn]extflow/).

The following example shows how to filter a channel by using a regular expression ~/^1.*/ inside a slashy string, that returns only strings that begin with 1:
```
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter(~/^1.*/)
  .view()
1
10
11
12
13
14
15
16
17
18
19
```

#### Boolean statement
A filtering condition can be defined by using a Boolean expression described by a closure `{}` and returning a boolean value. For example the following fragment shows how to combine a filter for a type qualifier Number with another filter operator using a Boolean expression to emit numbers less than 5:
```
channel
  .of( 1..22, 'X', 'Y' )
  .filter(Number)
  .filter { it < 5 }
  .view()
```
```
1
2
3
4
```

In the above example we could remove the brackets around the filter condition e.g. `filter{ it<5}`, since it specifies a closure as the operator’s argument. This is language short for `filter({ it<5})`

#### Literal value

Finally, if we only want to include elements of a specific value we can specify a literal value. In the example below we use the literal value X to filter the channel for only those elements containing the value X.
```
channel
  .of( 1..22, 'X', 'Y' )
  .filter('X')
  .view()
```
```
X
```
Add two channel filters to the Nextflow script below to view only the even numbered chromosomes.
```
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter( Number )
  .filter({ it % 2 == 0 })
  .view()
```

### Modifying the contents of a channel

If we want to modify the items in a channel, we can use transforming operators.

#### map

Applying a function to items in a channel

The map operator applies a function of your choosing to every item in a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a closure `{}` as shown in the example below:

```
chr = channel
  .of( 'chr1', 'chr2' )
  .map ({ it.replaceAll("chr","") })

chr.view()
```

Here the map function uses the groovy string function `replaceAll` to remove the chr prefix from each element.

1
2
We can also use the map operator to transform each element into a tuple.

In the example below we use the map operator to transform a channel containing fastq files to a new channel containing a tuple with the fastq file and the number of reads in the fastq file. We use the built in `countFastq` file method to count the number of records in a FASTQ formatted file.

We can change the default name of the closure parameter keyword from `it` to a more meaningful name file using `->`. When we have multiple parameters we can specify the keywords at the start of the closure, e.g. file, numreads ->.
```
fq_ch = channel
    .fromPath( 'data/yeast/reads/*.fq.gz' )
    .map ({ file -> [file, file.countFastq()] })
    .view ({ file, numreads -> "file $file contains $numreads reads" })
```
This would produce.
```
file data/yeast/reads/ref1_2.fq.gz contains 14677 reads
file data/yeast/reads/etoh60_3_2.fq.gz contains 26254 reads
file data/yeast/reads/temp33_1_2.fq.gz contains 20593 reads
file data/yeast/reads/temp33_2_1.fq.gz contains 15779 reads
file data/yeast/reads/ref2_1.fq.gz contains 20430 reads
[..truncated..]
```
We can then add a filter operator to only retain those fastq files with more than 25000 reads.
```
channel
    .fromPath( 'data/yeast/reads/*.fq.gz' )
    .map ({ file -> [file, file.countFastq()] })
    .filter({ file, numreads -> numreads > 25000})
    .view ({ file, numreads -> "file $file contains $numreads reads" })
file data/yeast/reads/etoh60_3_2.fq.gz contains 26254 reads
file data/yeast/reads/etoh60_3_1.fq.gz contains 26254 reads
```

#### Converting a list into multiple items

The flatten operator transforms a channel in such a way that every item in a list or tuple is flattened so that each single entry is emitted as a sole element by the resulting channel.

```
list1 = [1,2,3]
ch = channel
  .of(list1)
  .view()
```
```
[1, 2, 3]
```
```
ch =channel
    .of(list1)
    .flatten()
    .view()
```
The above snippet prints:
```
1
2
3
```
This is similar to the channel factory `Channel.fromList`.


#### Converting the contents of a channel to a single list item.

The reverse of the flatten operator is collect. The collect operator collects all the items emitted by a channel to a list and return the resulting object as a sole emission. This can be extremely useful when combining the results from the output of multiple processes, or a single process run multiple times.
```
ch = channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()
```
It prints a single value:
```
[1,2,3,4]
```
The result of the collect operator is a value channel and can be used multiple times.

