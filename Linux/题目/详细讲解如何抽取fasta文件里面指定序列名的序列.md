# 详细讲解如何抽取fasta文件里面指定序列名的序列

## 题目

要从fasta序列里面抽取指定序列，那么必须要首先了解什么是fasta格式：

序列文件的第一行是由大于号">"或分号";"打头的任意文字说明（习惯常用">"作为起始），用于序列标记。从第二行开始为序列本身，只允许使用既定的核苷酸或氨基酸编码符号（参见下表）。通常核苷酸符号大小写均可，而氨基酸常用大写字母。使用时应注意有些程序对大小写有明确要求。文件每行的字母一般不应超过80个字符。
下面是FASTA格式的氨基酸序列实例：

```
>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken 
ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID 
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA 
DIDGDGQVNYEEFVQMMTAK*
```

比如有下面这个fasta文件

```
>ID_1
abcd
abcd
abcd
>ID_2
efgh
efgh
efgh
>ID_3
ijkl
ijkl
ijkl
```

如果我想提取序列名为 ID_2的序列，那么我肯定需要写代码。用什么编程语言无所谓啦，我这里就用perl给大家讲解一下。我们先来回顾一下生信菜鸟团群里朋友的问题：

perl问题：大家好，我有个fasta文件，我想从里边把chrM匹配并把匹配的行打印出来 ：试写脚本如下

```perl
#!usr/bin/perl -w
open B,"fasta" or die "$!";
while (my $line = <B>) {
        if($line=~/chrM/i)
                {print $line}
}
```

没结果输出是何问题呢？？？？

他的代码错误在打印出来的就是chrM这个序列名而已，而不是序列。要做到提取序列名+序列，首先需要用编程的思想来解析一下整个任务！

1. 首先，要提取文件里面的序列，必须要用代码读取文件。
2. 然后，我们要根据某个特定的序列名来提取，所以我们碰到了>开头的序列名，就需要判断一下，是不是我们想要的。
3. 最后，我们对每一行都需要判断，是否需要打印出来，所以肯定有一个控制变量来决定是否要打印！

我最讨厌的就是写正规代码和注释，但是为了讲解清楚这个最基础代码，我破例一次，就在上面问问题的同学代码基础上面修改一下：

```perl
#!usr/bin/perl -w
open B,"fasta" or die "$!"; ##  要保证你的fasta文件的文件名就是fasta
my $sign=0; ## 这个是控制是否打印的变量
while (my $line = <B>) {  ## 一行行的读取
         chomp $line;  ## 去掉末尾的换行符
        if($line =~ /^>/){ # 判断是否是fasta序列名
                if ($line eq '>ID_2' ){  ## 如果是你需要的序列名
                        $sign=1;  ## 就把控制变量设为 1 ， 意思是可以输出啦
                }else{
                        $sign=0;  ## 如果不是你需要的序列名，就不输出
                }
        }
        print "$line\n" if $sign==1; ## 根据控制变量的情况来决定是否打印输出
}
```

如果是我自己写代码，我一般是这样写

```perl
$sign=0;
while(<>){
        chomp;
        if(/^>/){
                $_ eq '>ID_2' ? ($sign=1):($sign=0);
        }
        print "$_\n" if $sign==1;
}
```

## 我的答案

```
awk 'BEGIN{RS=">"}{if($0 ~ /ID_2/)print ">"$0}' test.fa 
```

---
:yum: 我还是一个小小白，平时的笔记一定有不少疏漏之处，恳请各位大佬批评指正！