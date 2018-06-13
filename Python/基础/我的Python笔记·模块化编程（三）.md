# 我的Python笔记·模块化编程（三）

构建程序流程

## 本文内容

- 如何使多个程序共同运行
- 如何用Python运行其他程序
- 如何将参数传递到Python程序中
- 如何用Python对文件目录进行浏览

## 构建NGS流程

对于很多研究来说，一个程序往往不足以解决问题，所以我们经常需要将多个程序组合起来共同使用，NGS数据分析就是一个十分典型的例子：

- 下机产生`.fastq`文件；
- 首先用`TopHat`读取读长`.fastq`，将其回帖到参考基因组；
- `TopHat`输出的`.bam`会被传递到下一个程序`Cufflinks`进行转录本的拼接，产生`.gtf`；
- 对于不同条件的样本或重复实验，从中获得的转录组文件`.gtf`可以先进行筛选，然后相互比较，按照唯一的参考转录组进行拼接，这一步可以使用`Cufflink`工具包中的`Cuffcompare`；`Cufflinks`由三个程序集合而成：转录组装配`Cufflinks`、转录组比较`Cuffcompare`、检测调控和表达差异`Cuffdiff`。

注：目前已经用Hisat2取代了TopHat流程


![pipeline](https://www.melbournebioinformatics.org.au/tutorials/tutorials/rna_seq_dge_basic/media/image01.png)

```python
import os

tophat_output_dir = '/home/RNA-seq/tophat'
tophat_output_file = 'accepted_hits.bam'
bowtie_index_dir = '/home/RNA-seq/index'
cufflinks_output_dir = '/home/RNA-seq/cufflinks'
cufflinks_output_file = 'transcripts.gtf'
illumina_output_file = 'sample.fastq'

tophat_command = 'tophat -o %s %s %s' %\
                 (tophat_output_dir, bowtie_index_dir,\
                  illumina_output_file)
os.system(tophat_command)

cufflinks_command = 'cufflinks -o %s %s%s%s' %\
                    (cufflinks_output_dir, tophat_output_dir, os.sep,\
                     tophat_output_file)
os.system(cufflinks_command)
```

## 程序是如何工作的

- `os`模块就是对操作系统进行操作，使用该模块必须先导入模块`import os`，例子中使用`os`模块中的`system()`函数来运行tophat和cufflinks程序；
- 导入模块后的6行都是在进行变量赋值；在8行和第10行，定义两条命令并将其作为字符串储存在变量中；
- `os.system()`的参数是一个单一字符串，包含了需要运行的shell命令；
- cufflinks将tophat的输出作为输入。第11行的`os.sep`的值是用于连接路径和文件名的字符。

### 一些技巧

#### 在程序中交换文件名和数据

当我们连接几个相互调用的程序时，他们需要交换数据，第二个程序需要知道第一个程序的输出文件名。下面是程序间传递信息的4种方法：

1. 用`raw.input()`输入文件名。这样每次运行程序时都会要求输入文件名或参数，但当使用程序好多次后，这样重复快速的输入相同的东西总是很烦人；
2. 使用Python变量。可以将代码中的文件名、参数甚至数据直接储存为字符串。这个方法的好处在于，所有的赋值都在同一块区域，但是当想要将流程应用于不同文件时都必须修改这个程序；
3. 将文件名储存在另一个文件中。我们还可以把文件名卸载一个独立的文件中，用程序打开并读取它。那么每次使用不同的文件时只要修改这个文件即可，如果程序很庞大，这种方法就会使代码易于模块化；
4. 使用命令行参数。可以使用`sys`模块中的`sys.argv`来传递命令行参数。

#### 使用`if`语句验证路径是否存在

```python
tophat_dir = '/home/RNA-seq/tophat'
index_dir = '/home/RNA-seq/index'

if os.path.exists(tophat_dir) and os.path.exists(index_dir):
    os.system('tophat -o ' + tophat_dir + ' ' + index_dir + \
              sample.fastq)
else:
    print 'You have to create tophat and/or index directories before running your wrapper'
```

如果路径真实存在，则`os.path.exists()`方法将返回`TRUE`，也可以吧所有路径变量放在一个单独的模块里，命名为`pathvariables.py`，这样当我们更改文件名或路径时，不必去查看整个程序。最常见的运行错误发生在文件位置和目录更改时忘记给路径变量赋值。要使用路径变量就必须将这个模块导入程序：

```python
from pathvariables import tophat_dir, index_dir
if os.path.exists(tophat_dir) and os.path.exists(index_dir):
    os.system('tophat -o ' + tophat_dir + ' ' + index_dir + \
              sample.fastq)
else:
    print 'You have to create tophat and/or index directories before running your wrapper'
```

`pathvariables.py`文件的内容：

```python
tophat_dir = '/home/RNA-seq/tophat'
index_dir = '/home/RNA-seq/index'
cufflinks_dir = '/home/RNA-seq/cufflinks'
```

#### 关闭文件时的延迟

通常，使用流程的一个问题就是：我们当我们需要在第一个程序运行完毕后再调用第二的程序，尤其是如果第一个程序的输出是第二个程序的输入时，但有时候会出现一些问题，比如书写输入文件出现了延迟，此时虽然下一步调用的输入文件还未写好，程序也会开始下一个子过程。

这里有一个小技巧可以避免这样的情况：可以在之后的调用前插入一个动作来确保上一步系统调用的完成。例如，可以打开一个傀儡文件（dummy）写点东西进去，然后关掉，接着用`os.path.exists()`来检查文件是否被创建：

```python
import sys
import os
sys.path.append('home/RNA-seq/')
from pathvariables import tophat_dir, index_dir, cufflinks_dir
os.system('tophat -o ' + tophat_dir + ' ' + index_dir + \
              sample.fastq)
# here we don't know whether the tophat output file is completed and available
# we open and close a dummy file, so the operating system catches up
lag_file = open('dummy.txt', 'w')
lag_file.write('tophat completed')
lag_file.close()
# read the output file
if os.path.exists('home/RNA-seq/dummy.txt'):
    os.system('cufflinks -o ' + cufflinks_dir + ' '\
        + tophat_dir + '/accepted_hits.bam')
```

#### 使用命令行参数

从UNIX命令行向Python程序传递参数，比如：

```bash
python ngs_pipeline.py dataset_one.fastq
python ngs_pipeline.py dataset_two.fastq
```

文件名`dataset_one.fastq`和`dataset_two.fastq`是如何传递到Python程序中的呢？

参数会自动保存在名为`sys.argv`的变量中，想要知道他是如何作用的，可以试试下面的代码：

```python
import sys
print(sys.argv)
```

现在试着在终端中使用不同的参数调用这个脚本：

```bash
python arguments.py
python arguments.py Hello
python arguments.py 1 2 3
python arguments.py -o sample.fastq
```

可以看到，输入的部分都变成了字符串储存在列表中，列表中第一个元素总是所调用程序的名称，所以第一个参数应该是`sys.argv[1]`。使用`sys.argv`列表中储存的信息就能将命令行参数传递至程序中，如果程序需要很多不同的命令行选项，则使用`optparse`模块。

#### 处理文件和路径

- `os.listdir(dirtectory)`将路径下所有文件读取至一个列表中；
- `os.chdir(path)`改变路径，因为有些程序要求在特定文件夹下启动；
- `os.mkdir()`新建文件夹，`os.rmdir()`移除文件夹；
- 使用`tempfile`模块来创建临时文件。
