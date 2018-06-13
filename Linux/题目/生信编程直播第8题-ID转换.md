# 生信编程直播第8题-ID转换

## 题目

- 基础知识参见：
    - [NCBI的基因entrez ID相关文件介绍](http://www.bio-info-trainee.com/75.html)
    - [ID转换大全](http://www.biotrainee.com/thread-862-1-1.html)
    - [生信人必须了解的各种ID表示方式](http://www.biotrainee.com/thread-43-1-1.html)

不管是什么ID转换，都是找到对应关系，然后match一下即可！
然后是练习题：

首先是affymetrix芯片的探针：运行里面的R代码，得到的变量my_probe就是我们需要转换的探针ID！

```R
rm(list=ls())
 
library("hgu95av2.db")
ls('package:hgu95av2.db')
probe2entrezID=toTable(hgu95av2ENTREZID)
probe2symbol=toTable(hgu95av2SYMBOL)
probe2genename=toTable(hgu95av2GENENAME)
 
my_probe = sample(unique(mappedLkeys(hgu95av2ENTREZID)),30)
 
tmp1 = probe2symbol[match(my_probe,probe2symbol$probe_id),]
tmp2 = probe2entrezID[match(my_probe,probe2entrezID$probe_id),]
tmp3 = probe2genename[match(my_probe,probe2genename$probe_id),]
 
write.table(my_probe,'my_probe.txt',quote = F,col.names = F,row.names =F)
write.table(tmp1$symbol,'my_symbol.txt',quote = F,col.names = F,row.names =F)
write.table(tmp2$gene_id,'my_geneID.txt',quote = F,col.names = F,row.names =F)
```

如果用perl和python来写，就把对应关系保存到本地文件，或者去其它地方下载对应关系，然后用hash或者dictionary来转换！

ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/  下载对应关系

然后是illumina芯片的 探针：

```R
library("illuminaHumanv4.db")
ls('package:illuminaHumanv4.db')
 
probe2entrezID=toTable(illuminaHumanv4ENTREZID)
probe2symbol=toTable(illuminaHumanv4SYMBOL)
probe2genename=toTable(illuminaHumanv4GENENAME)
 
my_probe = sample(unique(mappedLkeys(illuminaHumanv4ENTREZID)),30)
 
probe2symbol[match(my_probe,probe2symbol$probe_id),]
probe2entrezID[match(my_probe,probe2entrezID$probe_id),]
probe2genename[match(my_probe,probe2genename$probe_id),]
```

其实探针不管是什么平台，都是很简单的！

> 所有bioconductor支持的芯片平台对应关系：[通过bioconductor包来获取所有的芯片探针与gene的对应关系](http://www.bio-info-trainee.com/1399.html)
> 可以从NCBI的GPL平台下载：http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6947
> 也可以直接加载对应的包！或者直接去公司的主页下载**manifest文件**！

最后是基因的转换：
运行下面的R代码，得到的my_symbol_gene和my_entrez_gene就是我们需要转换的ID

```R
library("illuminaHumanv4.db")
ls('package:illuminaHumanv4.db')
my_entrez_gene = sample(unique(mappedRkeys(illuminaHumanv4ENTREZID)),30)
my_symbol_gene = sample(unique(mappedRkeys(illuminaHumanv4SYMBOL)),30)
 
library("org.Hs.eg.db")
ls('package:org.Hs.eg.db')
 
entrezID2symbol <- toTable(org.Hs.egSYMBOL)
 
entrezID2symbol[match(my_entrez_gene,entrezID2symbol$gene_id),]
entrezID2symbol[match(my_symbol_gene,entrezID2symbol$symbol),]
```

## 我的答案

1. 下载对应关系文件：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180130/dHLh0iD2D2.png?imageslim)

我就先下载了最小的 gene2go.gz 试试。

2. 解压看看：

```shell
gunzip gene2go.gz
```

```shell
head gene2go
```

```
#tax_id	GeneID	GO_ID	Evidence	Qualifier	GO_term	PubMed	Category
3702	814629	GO:0005634	ISM	-	nucleus	-	Component
3702	814629	GO:0008150	ND	-	biological_process	-	Process
3702	814630	GO:0003677	IEA	-	DNA binding	-	Function
3702	814630	GO:0003700	ISS	-	DNA binding transcription factor activity	11118137	Function
3702	814630	GO:0005634	IEA	-	nucleus	-	Component
3702	814630	GO:0005634	ISM	-	nucleus	-	Component
3702	814630	GO:0006351	IEA	-	transcription, DNA-templated	-	Process
3702	814630	GO:0006355	TAS	-	regulation of transcription, DNA-templated	11118137	Process
3702	814636	GO:0003674	ND	-	molecular_function	-	Function
```

3. 写个脚本：

既然是要对应第三列的GO_ID与第四列的gene_name，只要写个循环依次`grep` 到基因名所在行，然后`cut` 到第三列的GO_ID并输出即可~

```shell
#!/bin/bash
#
for i in `cat genelist.txt` ;do
    grep "\b$i\b" gene2go | cut -f3,4 >>tmp
done
```

4. 看了论坛上前辈的答案，学到了还可以用`grep -w -f`选项。

```shell
grep -w -f genelist.txt gene2go | cut f3,4 >>tmp
```

```
-w：只显示全字符合的列。
-f<范本文件>：指定范本文件，其内容有一个或多个范本样式，让grep查找符合范本条件的文件内容，格式为每一列的范本样式。
```

---

:yum: 我还是一个小小白，平时的笔记一定有不少疏漏之处，恳请各位大佬批评指正！

