# 强大的并行工具 GNU Parallel

你是否曾经计算过一个非常大的数据(几百GB)？或是在里面搜索，或其它操作（一些无法并行的操作比如`grep`, `bzip2`, `wc`, `awk`, `sed`），然而这些不能并行运算的软件都只能使用一个CPU内核，不能发挥出计算机的真正实力，所以往往会造成一核有难，N核围观的局面。

![](http://oo3g995ih.bkt.clouddn.com/blog/180813/EbE2a76K7i.png?imageslim)

本期就要强烈安利一个**真·并行**工具 —— `GNU Parallel`，它是一个shell工具，可以在一台或多台计算机上并行的执行计算任务，一个计算任务可以是一条shell命令或者一个以每一行做为输入的脚本程序。`GNU Parallel`会把输入分块，然后通过管道并行执行。

如果你会使用`xargs`和`tee`命令，你会发现`GNU Parallel`非常易于使用，因为`GNU Parallel`具有与`xargs`一样的选项。`GNU Parallel`可以替代大部分的shell循环，并且用并行的方式更快的完成计算任务。

## 安装 GNU Parallel

```
(wget -O - pi.dk/3 || curl pi.dk/3/) | bash
```

这条命令同时也会安装最新版的指南

```
man parallel_tutorial
```

当然也可以用`conda`来安装，更多请参考[官方文档](https://www.gnu.org/software/parallel/)。

## 简单的栗子

下面用一个简单的例子介绍一下`Parallel`的用法：在这个例子中，有10个蛋白质序列文件，我们想要用相同的选项运行`blastp`。使用下面的`parallel`并行命令，我们将一次对两个文件运行`blastp`，并为每个文件分配1个线程。

### 下载测试数据（45 mb）

```
wget http://oo3g995ih.bkt.clouddn.com/parallel_example.tar.gz
```

解压：

```
tar xfvz parallel_example.tar.gz
```

进入文件夹：

```
cd parallel_example
```

运行命令：

```
mkdir blastp_outfiles

parallel --eta -j 2 --load 80% --noswap 'blastp -db pdb_blast_db_example/pdb_seqres.txt -query {} -out blastp_outfiles/{.}.out -evalue 0.0001 -word_size 7 -outfmt "6 std stitle staxids sscinames" -max_target_seqs 10 -num_threads 1' ::: test_seq*.fas
```

你会发现其实只有四个序列匹配到了数据库中。

单引号的内容是`parallel`的一些选项。单引号内的命令包含仅用于`blastp`的选项（需要注意的是`-num_threads`，这是每个`blastp`命令使用的线程数）。

### 必须的选项

- 文件名：`{}`
- 删除了扩展名的文件名：`{.}`
例如：`test.fa`
- 表示从命令行读入以下所有内容：`:::`
例如：`parallel gzip ::: *`的意思是解压当前目录下的所有文件，但要是只输入`parallel gzip *`，这行命令并不会起作用。你需要加入`:::`

### 一些可选的选项

- `--eta`：显示任务完成的预计剩余时间。
- `-j 2`或`-jobs 2`：同时运行的命令数，在本例中设置为2。
- `--load 80%`：最大cpu负载。在上面的命令中，我们指定最多可以运行80％的CPU。
- `--noswap`：如果服务器处于大量内存负载下，则不会启动新作业，以便在存储新信息之前从内存中删除信息。

### 使用管道输入文件列表

在之前的例子中，我们在并行命令的末尾使用`:::`来指示输入文件。你也可以使用管道`|`来进行输入：

```
mkdir blastp_outfiles2

ls test_seq*.fas | parallel --eta -j 2 --load 80% --noswap 'blastp -db pdb_blast_db_example/pdb_seqres.txt -query {} -out blastp_outfiles2/{.}.out -evalue 0.0001 -word_size 7 -outfmt "6 std stitle staxids sscinames" -max_target_seqs 10 -num_threads 1'
```

### 从文件中读取管道命令

我们也可以将上面的命令写入一个文本文件中，然后通过管道输入到`parallel`中运行。

```
# 创建一个新文件夹
mkdir blastp_outfiles3
```

写一个`shell`脚本来产生命令：

```
for f in test_seq*.fas
do
out=${f/.fas/.out};
echo "blastp -db pdb_blast_db_example/pdb_seqres.txt -query $f -out blastp_outfiles3/$out -evalue 0.0001 -word_size 7 -outfmt \"6 std stitle staxids sscinames\" -max_target_seqs 10 -num_threads 1" >> blastp_cmds.txt
done
```

读入包含命令的文件并传递到`parallel`中运行：

```
cat blastp_cmds.txt | parallel --eta -j 2 --load 80% --noswap '{}'
```

在这种情况下，因为我们输入的是整个命令，所以可以用`{}`来引用输入本身。

## 改写for loop

经过上面的例子，我们可以把许多需要写循环的地方用`parallel`进行替换，以进行并行处理。

比如，当我们需要下载一些SRA数据，原来的脚本是这样滴：

```
for i in `seq 48 62`;
do
    prefetch SRR35899${i}
done
```

现在就可以改写成：

```
seq 48 62 | parallel --eta -j 10 'prefetch SRR35899{}'
```

> 整理自：https://github.com/LangilleLab/microbiome_helper/wiki/Quick-Introduction-to-GNU-Parallel
