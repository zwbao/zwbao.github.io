# 生信编程直播第7题(shell)-批量从NCBI下载数据

## 题目

首先进入https://www.ncbi.nlm.nih.gov/genome/genomes/13563 可以找到Mycobacterium相关的180个记录，现在需要批量下载每一个记录的4个数据：

https://www.ncbi.nlm.nih.gov/genome/13563?genome_assembly_id=305879
![mark](http://oo3g995ih.bkt.clouddn.com/blog/180129/gG2dikcKaE.png?imageslim)

回到搜索界面，找到那个download table的控件点击即可！这样就可以拿到里面的ftp列表！再把自己的180个匹配出来，cut有ftp地址的那一列，拿到ftp地址。

然后写一个shell脚本即可~

## 答案

1. 获取ftp地址：

```shell
cut -f19,20 genomes_proks.txt | sed '/^R/d' >ftp.list
```

2. 下载：

- 脚本如下：

```shell
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.fna.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.gff.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.gbff.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.faa.gz $id;done
```

- `wget`的一些选项：

```shell
-A<后缀名>：指定要下载文件的后缀名，多个后缀名之间使用逗号进行分隔；
-c：断点续传，当文件特别大或者网络特别慢的时候，往往一个文件还没有下载完，连接就已经被切断，此时就需要断点续传；
-r -A：下载指定格式文件；
-np：不遍历父目录；
-k：让下载得到的 HTML 或 CSS 中的链接指向本地文件；
-L：仅顺着关联的连接；
-p：下载所有用于显示 HTML 页面的图片之类的元素；
-nd 表示不在本机重新创建目录结构；
-i<文件>：从指定文件获取要下载的URL地址；
```
