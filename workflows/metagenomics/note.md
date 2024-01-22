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