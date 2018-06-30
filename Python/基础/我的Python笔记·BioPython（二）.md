# 我的Python笔记·BioPython（二）

![](https://biopython.org/assets/images/biopython_logo_white.png)

## 本文内容


## 问题

前一部分主要讲了如何用`Biopython`处理数据文件（如`FATSA`和`GenBank`），现在我们将学习用`Biopython`访问NCBI数据库（如`PubMed`和`GenBank`），以及`Expasy`资源（如`UniProt`），检索和解析他们的内容。

下面的例子是如何找到关于`PyCogent`的文献（`PyCogent`是与`Biopython`互补的库）：

```python
from Bio import Entrez
from Bio import Medline
keyword = "PyCogent"
# search publications in PubMed
Entrez.email = "my_email@address.com"
handle = Entrez.esearch(db="pubmed", term=keyword)
record = Entrez.read(handle)
pmids = record['IdList']
print(pmids)
# retrieve Medline entries from PubMed
handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
medline_records = Medline.parse(handle)
records = list(medline_records)
n = 1
for record in records:
    if keyword in record["TI"]:
        print(n, ')', record["TI"])
        n += 1
```

代码的输出：

```python
['22479120', '18230758', '17708774']
1 ) Abstractions, algorithms and data structures for structural bioinformatics in PyCogent.
2 ) PyCogent: a toolkit for making sense from sequence.
```

## 程序是如何工作的

### `Entrez`模块

`Biopython`访问NCBI资源的模块为`Entrez`。`Entrez`模块利用了 Entrez Programming Utilities（也称作 EUtils），包含八个工具，详情请见NCBI 的网站：http://www.ncbi.nlm.nih.gov/entrez/utils/. 每个工具都能在 Python 的`Bio.Entrez` 模块中找到对应函数。

可以输入以下命令来查看可用的方法和属性：

```python
from Bio import Entrez
print(dir(Entrez))
```

在这个输出中可以看到例子中用到的`email`属性和`esearch()`，`efetch()`函数。

- 使用 email 参数，这样如果遇到什么问题，NCBI 可以通过邮件联系到你。你可以在每次请求 `Entrez`的时候明确的设置这个参数（例如，在参数列表中包含 `email="my_email@address.com"` ），也可以设置一个全局的 email 地址：

```python
>>> from Bio import Entrez
>>> Entrez.email = "my_email@address.com"
```

#### `esearch()`

我们可以使用`esearch()`来搜索任意的数据库。这个函数有两个强制性参数：
- `db`：搜索的数据库（默认是`pubmed`）；
- `term`：查询的文本。

在例子中我们查询的关键词是`PyCogent`，如果要搜索超过一个关键词也可以使用`AND`和`OR`，也能用关键词类别，如`[Year]`，`[Organism]`，`[Gene]`等，比如`handle = Entrez.esearch(db="nucleotide",term="Cypripedioideae[Orgn] AND matK[Gene]")`。`Bio.Entrez.esearch()`的结果可以用`Entrez.read()`读取，返回一个字典。

可选参数`retmax`来设置返回多少条数据。其他参数包括`datatype`，`reldate`，`mindate`和`maxdate`后两个为`YYYY/MM/DD`格式。`datatype`用于选择日期类型（`mdat`：修改日期；`pdat`：发表日期；`edat`：Entrez日期）；`reldate`必须是一个整数，表示那些匹配`datatype`在`n`天之内的数据记录ID。`mindate`和`maxdate`是日期范围。

> 更多内容请参见 ESearch 帮助页面：https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch

搜索到结果后，我们可以通过`efetch()`来获得关于每条记录的详细内容。

#### `efetch()`

在上面的例子中当我们获得了`pmids`后，我们就可以用`efetch()`来得到详细信息，试着加入一句`print(handle.read())`：

```python
from Bio import Entrez
from Bio import Medline
keyword = "PyCogent"
# search publications in PubMed
Entrez.email = "my_email@address.com"
handle = Entrez.esearch(db="pubmed", term=keyword)
record = Entrez.read(handle)
pmids = record['IdList']
# print(pmids)
# retrieve Medline entries from PubMed
handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
print(handle.read())
# medline_records = Medline.parse(handle)
# records = list(medline_records)
# n = 1
# for record in records:
#     if keyword in record["TI"]:
#         print(n, ')', record["TI"])
#         n += 1
```

`efetch()`函数有许多可选参数：`retmode`指出检索记录格式，`rettype`指出记录类型，这取决于要访问的数据库。`PubMed`的`rettype`值可以是`abstract`，`citation`，`medline`等，对于`UniProt`可以设置为`fasta`来检索蛋白质的序列；`retmax`是返回记录的总数。

再举一个例子，用`efetch()`下载序列数据，保存到本地，然后并用`SeqIO`来解析数据：

```python
import os
from Bio import SeqIO
from Bio import Entrez
Entrez.email = "my_email@address.com" # Always tell NCBI who you are
filename = "gi_186972394.gbk"
if not os.path.isfile(filename):
# Downloading...
    net_handle = Entrez.efetch(db="nucleotide",id="186972394",rettype="gb", retmode="text")
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print("Saved")
print("Parsing...")
record = SeqIO.read(filename, "genbank")
print(record)
```

> 更多内容请参见 EFetch 帮助页面：https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

### `Medline`模块

要解析`efetch()`下载的`PubMed`数据，就需要导入`BioPython`的`Medline`模块，他提供了`Medline.parse()`函数，这个函数结果可以方便地转换成一个列表。这个列表包含了字典（在上面的例子中检索到了三条记录嘛），插入一句`print(records[0].keys())`，看看有哪些可用的键：

```python
dict_keys(['PMID', 'OWN', 'STAT', 'LR', 'IS', 'VI', 'IP', 'DP', 'TI', 'PG', 'AB', 'FAU', 'AU', 'LA', 'GR', 'PT', 'DEP', 'PL', 'TA', 'JT', 'JID', 'PMC', 'EDAT', 'MHDA', 'CRDT', 'PHST', 'AID', 'PST', 'SO'])
```

最常用的键是`TI`（title），`PMID`，`PG`（页面），`AB`（摘要）和`AU`（作者）
，当然也不是每一个字典都包括所有的键。

注意，如果需要解析单条记录，则可以使用`Medline.read()`函数而不是`Medline.parse()`。

### 获取Entrez数据库的信息

`Entrez.einfo()`函数为每个 NCBI 的数据库提供了条目索引，最近更新的时间以及可用的链接。此外，你可以很容易地使用 EInfo 通过 Entrez 获取所有数据库名字的列表：

```python
from Bio import Entrez
Entrez.email = "my_email@address.com" # Always tell NCBI who you are
handle = Entrez.einfo()
info = handle.read()
print(info)
```

输出结果如下：

```python
<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD einfo 20130322//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20130322/einfo.dtd">
<eInfoResult>
<DbList>

	<DbName>pubmed</DbName>
	<DbName>protein</DbName>
	<DbName>nuccore</DbName>
	<DbName>ipg</DbName>
	<DbName>nucleotide</DbName>
	<DbName>nucgss</DbName>
	<DbName>nucest</DbName>
	<DbName>structure</DbName>
	<DbName>sparcle</DbName>
	<DbName>genome</DbName>
	<DbName>annotinfo</DbName>
	<DbName>assembly</DbName>
	<DbName>bioproject</DbName>
	<DbName>biosample</DbName>
	<DbName>blastdbinfo</DbName>
	<DbName>books</DbName>
	<DbName>cdd</DbName>
	<DbName>clinvar</DbName>
	<DbName>clone</DbName>
	<DbName>gap</DbName>
	<DbName>gapplus</DbName>
	<DbName>grasp</DbName>
	<DbName>dbvar</DbName>
	<DbName>gene</DbName>
	<DbName>gds</DbName>
	<DbName>geoprofiles</DbName>
	<DbName>homologene</DbName>
	<DbName>medgen</DbName>
	<DbName>mesh</DbName>
	<DbName>ncbisearch</DbName>
	<DbName>nlmcatalog</DbName>
	<DbName>omim</DbName>
	<DbName>orgtrack</DbName>
	<DbName>pmc</DbName>
	<DbName>popset</DbName>
	<DbName>probe</DbName>
	<DbName>proteinclusters</DbName>
	<DbName>pcassay</DbName>
	<DbName>biosystems</DbName>
	<DbName>pccompound</DbName>
	<DbName>pcsubstance</DbName>
	<DbName>pubmedhealth</DbName>
	<DbName>seqannot</DbName>
	<DbName>snp</DbName>
	<DbName>sra</DbName>
	<DbName>taxonomy</DbName>
	<DbName>biocollections</DbName>
	<DbName>unigene</DbName>
	<DbName>gencoll</DbName>
	<DbName>gtr</DbName>
</DbList>

</eInfoResult>
```

我们可以把这个XML读入一个`Python`变量中：

```python
from Bio import Entrez
Entrez.email = "my_email@address.com" # Always tell NCBI who you are
handle = Entrez.einfo()
record = Entrez.read(handle)
print(record.keys())
```

现在`record`是有一个确定键值的字典：

```python
dict_keys(['DbList'])
```

这个键对应的值储存在上面XML文件里面包含的数据库名字的列表：

```python
['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'nucgss', 'nucest', 'structure', 'sparcle', 'genome', 'annotinfo', 'assembly', 'bioproject', 'biosample', 'blastdbinfo', 'books', 'cdd', 'clinvar', 'clone', 'gap', 'gapplus', 'grasp', 'dbvar', 'gene', 'gds', 'geoprofiles', 'homologene', 'medgen', 'mesh', 'ncbisearch', 'nlmcatalog', 'omim', 'orgtrack', 'pmc', 'popset', 'probe', 'proteinclusters', 'pcassay', 'biosystems', 'pccompound', 'pcsubstance', 'pubmedhealth', 'seqannot', 'snp', 'sra', 'taxonomy', 'biocollections', 'unigene', 'gencoll', 'gtr']
```

对于这些数据库可以用`EInfo`来获取更多的信息：

```python
from Bio import Entrez
Entrez.email = "my_email@address.com" # Always tell NCBI who you are
handle = Entrez.einfo(db="pubmed")
record = Entrez.read(handle)
print(record["DbInfo"]["Description"])
print(record["DbInfo"]["Count"])
print(record["DbInfo"]["LastUpdate"])
```

输出结果如下：

```python
PubMed bibliographic record
28579404
2018/06/23 15:58
```

通过`record["DbInfo"].keys()`可以获取储存在这个记录里的其他信息。这里最有用的就是`ESearch`可用的搜索值列表：

```python
from Bio import Entrez
Entrez.email = "my_email@address.com" # Always tell NCBI who you are
handle = Entrez.einfo(db="pubmed")
record = Entrez.read(handle)
for field in record["DbInfo"]["FieldList"]:
    print("%(Name)s, %(FullName)s, %(Description)s" % field)
```

输出的结果是一个很长的列表，但是这会告诉你在使用 PubMed 的时候，你可以通过`Jones[AUTH] `搜索作者，或者通过`Sanger[AFFL]` 将作者范围限制在 Sanger Centre。这会非常方便，特别是当你对某个数据库不太熟悉的时候。

> 更多内容请参见 Biopython Tutorial and Cookbook - Chapter 9  Accessing NCBI’s Entrez databases：http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc120
