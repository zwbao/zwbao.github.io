# [生信人的linux考试](http://www.bio-info-trainee.com/2900.html)

一、在任意文件夹下面创建形如 1/2/3/4/5/6/7/8/9 格式的文件夹系列。
```shell
mkdir -p 1/2/3/4/5/6/7/8/9
# -p或--parents 若所要建立目录的上层目录目前尚未建立，则会一并建立上层目录；
```

二、在创建好的文件夹下面，比如我的是 /Users/jimmy/tmp/1/2/3/4/5/6/7/8/9 ，里面创建文本文件 me.txt

```shell
cat >me1.txt
Go to: http://www.biotrainee.com/
I love bioinfomatics.
And you ?
```

三、在文本文件 me.txt 里面输入内容:

```
Go to: http://www.biotrainee.com/
I love bioinfomatics.
And you ?
```
前三题效果如下：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180122/3aBmFafEim.png?imageslim)

四、删除上面创建的文件夹 1/2/3/4/5/6/7/8/9 及文本文件 me.txt
```shell
rm -r 1
# -r或-R：递归处理，将指定目录下的所有文件与子目录一并处理；
```

五、在任意文件夹下面创建 folder1~5这5个文件夹，然后每个文件夹下面继续创建 folder1~5这5个文件夹：

```shell
mkdir -p folder_{1..5}/folder_{1..5}
```

效果如下：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180122/464KdmD3mh.png?imageslim)

六、在第五题创建的每一个文件夹下面都创建第二题文本文件 me.txt ，内容也要一样。

```shell
cat >me1.txt
Go to: http://www.biotrainee.com/
I love bioinfomatics.
And you ?
echo folder_{1..5}/folder_{1..5} | xargs -n 1 cp ./me1.txt
# echo 通过管道传递到 xargs 命令中
# -n 1 - 告诉 xargs 命令每个命令行最多使用一个参数，并发送到 cp 命令中。
```


七，再次删除掉前面几个步骤建立的文件夹及文件
```shell
rm -r folder_{1..5}
```
八、下载 http://www.biotrainee.com/jmzeng/igv/test.bed 文件，后在里面选择含有 H3K4me3 的那一行是第几行，该文件总共有几行。
```shell
wget http://www.biotrainee.com/jmzeng/igv/test.bed
cat -n test.bed | grep 'H3K4me3' | cut -f1 #第八行
cat test.bed | wc -l #共10行
```

九、下载 http://www.biotrainee.com/jmzeng/rmDuplicate.zip 文件，并且解压，查看里面的文件夹结构

```shell
wget http://www.biotrainee.com/jmzeng/rmDuplicate.zip
unzip rmDuplicate.zip
cd rmDuplicate
tree
```

十、打开第九题解压的文件，进入 rmDuplicate/samtools/single 文件夹里面，查看后缀为 .sam 的文件，搞清楚生物信息学里面的[SAM](https://vip.biotrainee.com/d/162-sam)/[BAM](https://vip.biotrainee.com/d/163-bam) 定义是什么。

- bam和sam的区别与一致：
    - sam是带有比对信息的序列文件（即告诉你这个reads在染色体上的位置等），用于储存序列数据（SAM  format is a generic format for storing large nucleotide sequence alignments. ）。
    - BAM is the compressed binary version of the Sequence Alignment/Map (SAM) format. 生物信息中的二进制文件主要是为了节约空间，计算机机可读。可以用samtools工具实现sam和bam文件之间的转化。
    - 二者都是fastq文件经过序列比对或者mapping后输出的格式（其储存的信息都是一致的）

十一、安装 samtools 软件
参考[centos上安装samtool具体步骤](http://blog.csdn.net/u011808596/article/details/78176535)
[samtools常用命令详解](https://www.plob.org/article/7112.html)

十二、打开后缀为BAM 的文件，找到产生该文件的命令。 提示一下命令是：

```
/home/jianmingzeng/biosoft/bowtie/bowtie2-2.2.9/bowtie2-align-s --wrapper basic-0 -p 20 -x /home/jianmingzeng/reference/index/bowtie/hg38 -S /home/jianmingzeng/data/public/allMouse/alignment/WT_rep2_Input.sam -U /tmp/41440.unp
```

 查看bam文件的header信息：

 ```
 samtools view -H tmp.rmdup.bam
 ```

十三题、根据上面的命令，找到我使用的参考基因组 /home/jianmingzeng/reference/index/bowtie/hg38 具体有多少条染色体。

```
samtools view -H tmp.rmdup.bam | cut -f2 |sort -u
```
有25条染色体，chr1..22，chrM，chrX，chrY

十四题、上面的后缀为BAM 的文件的第二列，只有 0 和 16 两个数字，用 cut/sort/uniq等命令统计它们的个数。

```
samtools view tmp.rmdup.bam | cut -f2 | grep '0' |wc -l
16
samtools view tmp.rmdup.bam | cut -f2 | grep '16' |wc -l
12
```
第二列：sum of flags，每个flag用数字来表示，分别为：
```
1 read是pair中的一条（read表示本条read，mate表示pair中的另一条read）
2 pair一正一负完美的比对上
4 这条read没有比对上
8 mate没有比对上
16 这条read反向比对
32 mate反向比对
64 这条read是read1
128 这条read是read2
256 第二次比对
512 比对质量不合格
1024 read是PCR或光学副本产生
2048 辅助比对结果
通过这个和可以直接推断出匹配的情况。假如说标记不是以上列举出的数字，比如说83=（64+16+2+1），就是这几种情况值和。
```
十五题、重新打开 rmDuplicate/samtools/paired 文件夹下面的后缀为BAM 的文件，再次查看第二列，并且统计

```
samtools view tmp.rmdup.bam | cut -f2 | grep '0' |wc -l
0
samtools view tmp.rmdup.bam | cut -f2 | grep '16' |wc -l
2
```

十六题、下载 http://www.biotrainee.com/jmzeng/sickle/sickle-results.zip 文件，并且解压，查看里面的文件夹结构， 这个文件有2.3M，注意留心下载时间及下载速度。

```
wget http://www.biotrainee.com/jmzeng/sickle/sickle-results.zip 
unzip sickle-results.zip 
cd sickle-results
tree
```

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180123/6gClI5GB00.png?imageslim)

十七题、解压 sickle-results/single_tmp_fastqc.zip 文件，并且进入解压后的文件夹，找到 fastqc_data.txt 文件，并且搜索该文本文件以 >>开头的有多少行？

```
cat fastqc_data.txt | grep '^>>' | wc -l
24
```

十八题、下载 http://www.biotrainee.com/jmzeng/tmp/hg38.tss 文件，去NCBI找到TP53/BRCA1等自己感兴趣的基因对应的 refseq数据库 ID，然后找到它们的hg38.tss 文件的那一行。
https://www.ncbi.nlm.nih.gov/gene/7157

```
cat hg38.tss | grep 'NM_000546'
NM_000546	chr17	7685550	7689550	1
```

十九题、解析hg38.tss 文件，统计每条染色体的基因个数。

- 统计每条染色体的基因个数：

```
cut -f2 hg38.tss | cut -d"_" -f1 | sort | uniq -c | sort -rn
```

```
   6157 chr1
   5880 chr19
   5782 chr6
   4090 chr2
   3794 chr17
   3577 chr11
   3395 chr3
   3014 chr12
   2838 chr10
   2821 chr5
   2785 chr7
   2696 chr16
   2561 chrX
   2377 chr15
   2310 chr9
   2277 chr4
   2221 chr8
   1982 chr14
   1692 chr20
   1410 chr22
   1133 chr13
    895 chr21
    883 chr18
    414 chrY
     32 chrUn
      2 chrM
```

- chrUN 指的是不能准确地放到特定染色体上的基因组DNA克隆重叠群

二十题、解析hg38.tss 文件，统计NM和NR开头的序列，了解NM和NR开头的含义。

```
[root@VM_180_248_centos test]# echo `cat hg38.tss | grep "^NM" | sort -u |wc -l`
51064
[root@VM_180_248_centos test]# echo `cat hg38.tss | grep "^NR" | sort -u |wc -l`
15954

```

- NM：有意义的，可以转录成蛋白质的基因
- NR：不能转录成蛋白质的基因


