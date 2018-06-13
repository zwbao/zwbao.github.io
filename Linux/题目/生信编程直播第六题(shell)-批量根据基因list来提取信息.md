# 生信编程直播第六题(shell)-批量根据基因list来提取信息

## 题目


如果只有一个基因list，那么从大的文件里面提取信息，很容易。比如：[我的答案·详细讲解如何抽取fasta文件里面指定序列名的序列](https://vip.biotrainee.com/d/277 )

如果有多个多个gene list，就需要写循环啦，这个时候就用shell啦~

```shell
ls *txt |while read id
do
cat $id ~/annotation/CHIPseq/mm10/ucsc.refseq.bed |perl -alne '{$h{$F[0]}=1;print if exists $h{$F[3]} }' >${id%%.*}.bed
done
```

自行去UCSC下载refseq的bed，然后自己找几个gene list去这个refseq.bed里面提取自己的gene list的bed文件。

## 我的答案

1. 下载 bed 文件


- UCSC主界面如图所示，找到Table Browser点击进入；

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180129/K71mGaBl7f.png?imageslim)

- 在Table Browser里，选定人的基因组，采用最新的GRCh38版本，然后再选择Gene and Gene Predictions里的NCBI RefSeq作为想要导出的本地数据库。在导出格式里，我们选择了比较常用的BED格式，然后点击get output。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180129/dbmG7j3gb4.png?imageslim)

- 在Create one BED record per下面有一些选项，比如这里默认是Whole Gene，当然我们也可以选择启动子区域、外显子加周边区域、5' UTR区域、3' UTR区域等生成我们想要的BED文件。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180129/a5G384IHaf.png?imageslim)

- 参见：
    - [得到一个物种所有基因的TSS(转录起始位点)区域的bed文件](http://www.bio-info-trainee.com/2494.html)
    - [NGS数据格式之bed](https://mp.weixin.qq.com/s?__biz=MzUzMTEwODk0Ng==&mid=2247484331&idx=1&sn=3d2b2e644e1423e282318ef1d3bf3cec&scene=21#wechat_redirect)

2. 写个脚本：

首先用`cat`载入基因名，再`grep` 得到 bed 文件中匹配到的行，输出即可。

```
#!/bin/bash
#
for i in `cat genelist.txt` ;do
    grep "$i" test.bed >>genelist.bed
done
```

---

:yum: 我还是一个小小白，平时的笔记一定有不少疏漏之处，恳请各位大佬批评指正！

