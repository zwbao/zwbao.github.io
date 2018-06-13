# 生信编程直播第六题-下载最新版的KEGG信息，并且解析好

## 题目
[打开官网](http://www.genome.jp/kegg-bin/get_htext?hsa00001+3101) 可以在里面找到下载链接

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180128/ciIl850G8A.png?imageslim)

下载得到文本文件，可以看到里面的结构层次非常清楚：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180128/FG648Bmi8G.png?imageslim)

C开头的就是kegg的pathway的ID所在行，D开头的就是属于它的kegg的所有的基因；A,B是kegg的分类，总共是6个大类，42个小类。

```
grep ^A hsa00001.keg 
```

```
A<b>Metabolism</b>
A<b>Genetic Information Processing</b>
A<b>Environmental Information Processing</b>
A<b>Cellular Processes</b>
A<b>Organismal Systems</b>
A<b>Human Diseases</b>
```

也可以看到，到目前为止（2018年1月27日23:10:45），共有364个kegg的pathway信息啦~

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180128/0Kfab8b6c7.png?imageslim)

## 我的答案

1. 解析KEGG信息

写一个脚本：首先把`grep`得到含有pathwayID的行，然后用`awk`的切分功能把第二列提取出来，即pathway的ID，这样就得到了pathwayID的列表，然后根据这个列表`for`循环；在循环中，`grep -A<显示列数>` 显示除了显示符合范本样式的那一行之外，并显示该行之后的内容，之所以是27990是因为我之前`wc -l`看了下hsa00001.keg总共是27990行，这样得到这个pathwayID以下的内容，再用`awk`来取到基因，其中 `-v`用于传递变量，又因为`grep -A<显示列数>` 包括匹配到的行（即有pathwayID的行），所以要去掉，故`NR>1`，因为基因都是D开头，所以遇到D开头的就打印，一旦遇到不是D开头的，`exit`马上退出`awk`。

```
#!/bin/bash
#
for i in `grep ^C hsa00001.keg | awk  '{print $2}'` ;do
    grep -A27990 "$i" hsa00001.keg | awk -v kegg="$i" 'NR>1{if($0 ~ /^D/)printf kegg"\t"$2"\n";else exit}'
done
```

2.  运行脚本

```
[root@VM_180_248_centos d6t]# ./tmp.sh >kegg2gene.txt
[root@VM_180_248_centos d6t]# head kegg2gene.txt 
00010	3101
00010	3098
00010	3099
00010	80201
00010	2645
00010	2821
00010	5213
00010	5214
00010	5211
00010	2203
[root@VM_180_248_centos d6t]# cut -f1 kegg2gene.txt | sort -u | wc -l
382
[root@VM_180_248_centos d6t]# cut -f2 kegg2gene.txt | sort -u | wc -l
5591
```

3. 翻车的地方：这次的作业，很快就有了思路，马上就写好了脚本，如下：

```
#!/bin/bash
#
for i in `grep ^C hsa00001.keg | awk  '{print $2}'` ;do
    grep -A27990 "$i" hsa00001.keg | awk -v kegg="$i" 'NR>1{if($0 ~ /^D/)printf $kegg"\t"$2"\n";else exit}'
done
```

但是运行起来却不是我想象的那样，而是出现了如下结果：

```
[root@VM_180_248_centos d6t]# ./tmp2.sh >mistake.txt
[root@VM_180_248_centos d6t]# head -n20 mistake.txt 
	3101
	3098
	3099
hexokinase	80201
	2645
[EC:5.3.1.9]	2821
[EC:2.7.1.11]	5213
[EC:2.7.1.11]	5214
1	5211
[EC:3.1.3.11]	2203
[EC:3.1.3.11]	8789
aldolase,	230
aldolase,	226
aldolase,	229
isomerase	7167
dehydrogenase	2597
dehydrogenase,	26330
kinase	5232
kinase	5230
phosphoglycerate	5223
```

那么错误到底出在哪里呢？我苦思冥想，重新整理逻辑，重新写脚本，重新去翻以前的笔记……

犯的错误就是：

[SPOILER]打印变量不用加`$`！！！[/SPOILER]

---

:yum: 我还是一个小小白，平时的笔记一定有不少疏漏之处，恳请各位大佬批评指正！
