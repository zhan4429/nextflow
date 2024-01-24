## main.nf
### Shebang
```
#!/usr/bin/env nextflow
```
同很多Unix-like脚本一样，第一行叫做`shebang (Hash bang)`，出现在脚本第一行并以`#!`开头。它告诉系统用什么环境软件去解析这个脚本，当它存在并且脚本可执行的时候，我们可以通过直接调用该脚本来运行程序。以下为示例：

通过Interpreter来调用脚本（shebang不存在时也可使用）：
```
nextflow run ./main.nf
```

直接调用（需要shebang）：
```
chmod u+x ./main.nf ## 添加可执行权限
./main.nf           ## 系统自动使用nextflow运行
```

### 使用DSL2语言
```
nextflow.preview.dsl=2
```
注：DSL2是新功能，后续语法可能会调整

### 流程参数
```
// parameters
params.help = false
params.read_path  = "${workflow.projectDir}/data"
```
Nextflow通过`params`这个字典来允许执行时直接传入参数。上述的两个参数help，和read_path，在命令行中可通过`--help`，`--read_path /PATH/TO/READS`，来更改。

```
./main.nf --help ## 将params.help的值更改为true
./main.nf --read_path /data/reads/ ## 将params.read_path的值更改为~/reads
```

### Groovy原生支持
```
def helpMessage() {
    log.info"""
    =================================================================
    Usage: ${workflow.projectDir}/main.nf --read_path PATH/OF/READS
    =================================================================
    """.stripIndent()
}
if (params.help){
    helpMessage()
    exit 0
}
```

Nextflow基于Groovy语言，可在流程中直接使用。以上部分定义了函数helppMessage，**在接到--help输入时会输出帮助文档并退出执行**。workflow为Nextflow定义的特殊字典，`${workflow.projectDir}`对应了当前脚本(main.nf)的目录。

```
$ ./main.nf --help
N E X T F L O W  ~  version 19.09.0-edge
Launching `./main.nf` [distraught_meitner] - revision: bb88634790
WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
=================================================================
Usage: /home/ubuntu/shotgunmetagenomics-nf/tutorial/main.nf --read_path PATH/OF/READS
=================================================================
```

### Channel
```
ch_reads = Channel
    .fromFilePairs(params.read_path + '/**{1,2}.f*q*', flat: true)

ch_reads.view()
```
Netflow定义了Channel，通过Channel来在不同过程(Process)之间传递文件。Channel保证了基于不同文件的运算独立、并行进行。Channel对应的文件会被“拷贝”(经常以symbolic link的形式)到work/目录下。这个定义其实不是很容易说得清楚，建议阅读Nextflow帮助文档。

- `fromFilePair`专门为fastq打造，可以直接将不同的样本以列表形式分组
- `'/**{1,2}.f*q*'`定义了文件匹配方式，`**`表示递归地检索文件，`{1,2}.f*q*`跟bash本身的文件匹配一致，这里会匹配结尾为1.fastq.gz，2.fastq，1.fq，2.fq.gz的文件
- `flat: ture`使返回的双端序列和匹配ID存储为一个列表中。下面是一个例子：
```
$ ls /data/reads/
SRR1950772_1.fastq.gz  SRR1950772_2.fastq.gz  SRR1950773_1.fastq.gz  SRR1950773_2.fastq.gz

$ ./main.nf --read_path /data/reads/
N E X T F L O W  ~  version 19.09.0-edge
Launching `./main.nf` [elegant_volta] - revision: bb88634790
WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
[SRR1950773, /data/reads/SRR1950773_1.fastq.gz, /data/reads/SRR1950773_2.fastq.gz]
[SRR1950772, /data/reads/SRR1950772_1.fastq.gz, /data/reads/SRR1950772_2.fastq.gz]
```
如果默认情况 `(flat: false)`使用``.fromFilePairs(params.read_path + '/**{1,2}.f*q*')`，输出结果会将Read1和Read2合并在一个列表中：
```
[SRR1950773, [/data/reads/SRR1950773_1.fastq.gz, /data/reads/SRR1950773_2.fastq.gz]]
[SRR1950772, [/data/reads/SRR1950772_1.fastq.gz, /data/reads/SRR1950772_2.fastq.gz]]
```

## decont.nf

### 模块参数
```
params.index = 'hg19.fa'
params.outdir = './'
```
和上一篇介绍的流程参数一致，定义了在该模块的参数（我推测这些参数的scope是local，即仅限于该模块，有待确定）。

### 过程定义
```
process DECONT {
    input:
    ...
    output:
    ...
    script:
    ...
}
```
以上为Nextflow过程(process)的骨架结构，主要由输入(input)，输出(output)，脚本(script)来构成。

### 过程指令
```
tag "${prefix}"
cpus 8
publishDir params.outdir, mode: 'copy'
```
- tag：给每一个过程执行命名，方便在执行日志中查看
- cpus：此过程运行时的CPU数量
- publishDir：结果发布路径，运行完成后将最终的结果（由output定义）拷贝（'copy'）到该路径

### 输入、输出
```
input:
    file index_path
    tuple prefix, file(reads1), file(reads2)

output:
    tuple prefix, file("${prefix}*1.fastq.gz"), file("${prefix}*2.fastq.gz")
    tuple file("${prefix}.html"), file("${prefix}.json")
```
由关键字input:或output:(注意冒号)构成，后面接内容(全都为Channel类型)。内容可包含：
- 值(val)：可省略
- 文件(file)
- 列表(tuple)：形式为[element1, element2, elemnet3, ...]

### 脚本
```
script:
    """
    fastp -i $reads1 -I $reads2 --stdout -j ${prefix}.json -h ${prefix}.html | \\
    bwa mem -p -t $task.cpus ${index_path}/${params.index} - | \\
    samtools fastq -f12 -F256  -1  ${prefix}_fastpdecont_1.fastq.gz -2 ${prefix}_fastpdecont_2.fastq.gz -
    """
```
Nextflow过程的脚本部分为Shell脚本，注意Shell中变量要使用`\$`开头。脚本部分涉及的文件(或路径)必须是Channel。

### 配置DECONT过程参数
```
// parameters decont                                                                       // ***
params.decont_refpath = '/data/nucleotide/'                                                // ***
params.decont_index   = 'hg19.fa'                                                          // ***
params.decont_outdir  = './pipeline_output/decont_out'                                     // ***
ch_bwa_idx = file(params.decont_refpath)   
```                                              
我们知道DECONT部分需要四个参数：

- BWA索引路径
- BWA索引的前缀
- DECONT输出路径
- 输入fastq文件：[prefix, r1.fastq, r2.fastq] (tuple prefix, file(reads1), file(reads2))
注意这里BWA索引路径需要以Channel的形式传入。

### 引用模块
```
include './decont' params(index: "$params.decont_index", outdir: "$params.decont_outdir")  // ***
```
非文件或路径的参数通过`params(param1:, param2:)`传入。

### 流程定义
```
workflow{                                                                                 // ***
    DECONT(ch_bwa_idx, ch_reads)                                                          // ***
}
```
执行流程，process所需的输入文件通过上面的方式传入。

## kraken2.nf
```
params.outdir = './'

process KRAKEN2 {
    tag "${prefix}"
    cpus 8
    publishDir params.outdir, mode: 'copy'

    input:
    file index_path
    tuple prefix, file(reads1), file(reads2)

    output:
    tuple prefix, file("${prefix}.kraken2.report")
    tuple prefix, file("${prefix}.kraken2.tax")
    file "${prefix}.kraken2.out"

    script:
    """
    kraken2 \\
    --db $index_path \\
    --paired \\
    --threads $task.cpus \\
    --output ${prefix}.kraken2.out \\
    --report ${prefix}.kraken2.tax \\
    $reads1 $reads2 \\
    --use-mpa-style \\

    ### run again for bracken
    kraken2 \\
    --db $index_path \\
    --paired \\
    --threads $task.cpu \\
    --report ${prefix}.kraken2.report \\
    $reads1 $reads2 \\
    --output -
    """
}

process BRACKEN {
    tag "${prefix}"
    publishDir params.outdir, mode: 'copy'

    input:
    file index_path  // This must have the bracken database
    tuple prefix, file(kraken2_report)
    each tax

    output:
    file "${prefix}*.tsv"

    script:
    """
    TAX=$tax; \\

    bracken -d $index_path \\
    -i $kraken2_report \\
    -o ${prefix}.bracken.${tax} \\
    -l \${TAX^^}; \\

    sed 's/ /_/g' ${prefix}.bracken.${tax} | \\
    tail -n+2 | \\
    cut -f 1,7 > ${prefix}.bracken.${tax}.tsv
    """
}
```
其中`each`相当于上面脚本的for循环，对于不同的分类s（species），g（genus）分别运行bracken。
### 配置参数
```
// parameters kraken2                                              // ***
params.kraken2_refpath = '/data/minikraken2_v2_8GB_201904_UPDATE/' // ***
params.kraken2_outdir = './pipeline_output/kraken2_out'            // ***
ch_kraken_idx = file(params.kraken2_refpath)                       // ***
```
此处于DECONT过程的配置类似，这里假设Bracken的索引和Kraken2的索引在同一个目录下。

引用模块
```
include './kraken2' params(outdir: "$params.kraken2_outdir")      // ***
```

### 流程定义
```
workflow{
    DECONT(ch_bwa_idx, ch_reads)
    KRAKEN2(ch_kraken_idx, DECONT.out[0])                          // ***
    BRACKEN(ch_kraken_idx, KRAKEN2.out[0], Channel.from('s', 'g')) // ***
}
```
在上一篇博文，DECONT的输出被定义为两个tuple：
```
output:
    tuple prefix, file("${prefix}*1.fastq.gz"), file("${prefix}*2.fastq.gz")
    tuple file("${prefix}.html"), file("${prefix}.json")
```
第一个tuple为去宿主DNA之后的fastq文件，第二个为fastp的质控文件，运行KRAKEN过程的时候，我们以`DECONT.out[0]`来使用DECONT过程结果的第一个tuple。同样的，在BRACKEN过程调用时，我们以`KRAKEN2.out[0]`来使用KRAKEN2过程结果的第一个tuple.

### 执行

```
$ ls data
$ ./main.nf
N E X T F L O W  ~  version 19.09.0-edge
Launching `./main.nf` [jolly_heyrovsky] - revision: b913fda572
WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
executor >  local (8)
[6f/7c66ef] process > DECONT (SRR1950773)  [100%] 2 of 2 ✔
[51/9f4f92] process > KRAKEN2 (SRR1950773) [100%] 2 of 2 ✔
[c7/603e5e] process > BRACKEN (SRR1950772) [100%] 4 of 4 ✔
Completed at: 17-Oct-2019 12:57:14
Duration    : 3m 50s
CPU hours   : 0.5
Succeeded   : 8
```

### 查看结果

```
$ ls pipeline_output/kraken2_out/
SRR1950772.bracken.g.tsv  SRR1950772.kraken2.out     SRR1950772.kraken2.tax    SRR1950773.bracken.s.tsv  SRR1950773.kraken2.report
SRR1950772.bracken.s.tsv  SRR1950772.kraken2.report  SRR1950773.bracken.g.tsv  SRR1950773.kraken2.out    SRR1950773.kraken2.tax
```


