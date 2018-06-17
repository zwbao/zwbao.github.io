# 我的Python笔记·BioPython（一）

## 简介

`Biopython`（https://biopython.org/）是一个计算分子生物学模块的集合，用它可以实现许多生物信息学项目中的基本任务，比如：
- 解析文件格式，这些信息包括如基因和蛋白质序列、蛋白质结构、PubMed记录等等；
- 从资源库下载文件，这些资源包括NCBI，ExPASy等；
- 运行（本地或远程）常用的生物信息学算法，如BLAST，ClustalW等；
- 运行`Biopython`实现的算法，进行聚类、机器学习、数据分析和可视化。

我们也不仅可以用`Biopython`建立一个研究流程，也可以为特定的任务写新的代码，而让`Biopython`进行更多标准操作，你还可以修改`Biopython`的开源代码以更好地适应你的需求。

## 让我们开始吧

### 安装`Biopython`

直接使用`pip`安装：

```python
pip install biopython
```

升级旧版本：

```python
pip install biopython --upgrade
```

### 一个简单的示例

在前面的笔记中我们已经学习了如何用基本的字符串操作来解析序列文件，现在让我们用`Biopython`来试试看。`Biopython`提供了快速处理序列文件的工具，使利用不同文件格式、注释序列记录和把他们写入文件等操作变得非常简单。下面实例我们会用到4个从`Bio`导入的模块：`Seq`用来创建序列对象；`IUPAC`用来定义一个序列对象用的字符集（如DNA或蛋白质）；`SeqRecord`允许创建一个包含ID、注释、描述等的序列记录数据；`SeqIO`用来读写格式化的序列文件。

```python
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
# read the input sequence
dna = open("hemoglobin-gene.txt").read().strip()
dna = Seq.Seq(dna, IUPAC.unambiguous_dna)
# transcribe and translate
mrna = dna.transcribe()
protein = mrna.translate()

# write the protein sequence to a file
protein_record = SeqRecord(protein, id='sp|P69905.2|HBA_HUMAN', description="Hemoglobin subunit alpha, human")
outfile = open("HBA_HUMAN.fasta", "w")
SeqIO.write(protein_record, outfile, "fasta")
outfile.close()
```

`hemoglobin-gene.txt`文件内容如下：

```
ATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAA
```

`HBA_HUMAN.fasta`输出文件如下：

```
>sp|P69905.2|HBA_HUMAN Hemoglobin subunit alpha, human
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
AVHASLDKFLASVSTVLTSKYR*
```

### 程序是如何工作的

程序首先从文本文件中读取序列，然后将它转录成mRNA序列，翻译成肽段序列，最后写入`HBA_HUMAN.fasta`输出文件。让我们逐步分析导入Python会话3中的对象们。

#### `Seq`对象

`Seq.Seq`类创建了一个**序列**对象，用户在创建`Seq`对象时可以指定（或不指定）其字符集。

```python
>>> from Bio.Seq import Seq
>>> my_seq = Seq("AGTACACTGGT")
>>> my_seq
Seq('AGTACACTGGT', Alphabet())
>>> my_seq.alphabet
Alphabet()
```

当我们指定字符集时：

```python
>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import IUPAC
>>> my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
>>> my_seq
Seq('AGTACACTGGT', IUPACUnambiguousDNA())
>>> my_seq.alphabet
IUPACUnambiguousDNA()
```

当然，这也许是个氨基酸序列：

```python
>>> from Bio.Seq import Seq
>>> from Bio.Alphabet import IUPAC
>>> my_seq = Seq("AGTACACTGGT", IUPAC.unambiguous_dna)
>>> my_seq
Seq('AGTACACTGGT', IUPACUnambiguousDNA())
>>> my_seq.alphabet
IUPACUnambiguousDNA()
```

`Biopython`包含了一个预编译好的字符集，覆盖了所有生物序列类型。最频繁使用的是IUPAC定义的字符集(http://www. chem. qmw. ac. uk/iupac)。如果用户需要使用字符集，就必须从`Bio. Alphabet`模块导入IUPAC模块。它包括字符集IUPACUnamiguousDNA(基本的ACTG字母)，IUPACAmbiguousDNA(包含二义字母)，ExtendedIUPACDNA(包含修饰的碱基)，IUPACUnamiguousRNA,  IUPACAmbiguousRNA,  ExtendedIUPACRNA，IUPACProtein(IUPAC标准氨基酸)和ExtendedIUPACProtein(包括硒代半胱氨酸，X等)。在上面的实例中定义的dna变量是一个以IUPAC.unambiguous_dna字符集为特征的序列对象。

##### 转录和翻译序列

`Seq`对象可以用`transcribe()`方法和`translate()`方法来进行转录和翻译。`transcribe()`方法只是把所有的`T`替换成`U`，同时把字符集设置成RNA。我们也可以用`reverse_complement()`方法来得到反向互补序列，如：

```python
>>> from Bio import Seq
>>> my_seq = Seq.Seq("AGTACACTGGT")
>>> cdna = my_seq.reverse_complement()
>>> cdna
Seq('ACCAGTGTACT', Alphabet())
```

#### 把序列当成字符串工作

在`Biopython`中，我们可以像处理字符串一样处理序列对象。例如，可以索引、切片、分割、转换序列大小写，计算出现字符个数等等：

```python
>>> from Bio import Seq
>>> my_seq = Seq.Seq("AGTACACTGGT")
>>> my_seq[0]
'A'
>>> my_seq[0:3]
Seq('AGT', Alphabet())
>>> my_seq.split('T')
[Seq('AG', Alphabet()), Seq('ACAC', Alphabet()), Seq('GG', Alphabet()), Seq('', Alphabet())]
>>> my_seq.count('A')
3
>>> my_seq.count('A')/len(my_seq)
0.2727272727272727
```

需要注意的是将`Seq`对象分割后返回的是多个`Seq`对象，这些`Seq`对象可以用`+`连接。

除此之外，我们还可以用`find()`方法搜索序列中的子串，如果没有找到会返回`-1`，找到了则返回目标序列的最左端匹配字符的位置。当然，也可以配合`Python`的`re`模块或`Biopython`的`Bio.motif`模块用正则表达式进行搜索。

#### `SeqRecord`对象

`SeqRecord`类提供序列及其注释的容器，在上面的实例中我们将翻译得到的`protein`变量转换成`SeqRecord`对象。

```python
protein_record = SeqRecord(protein, id='sp|P69905.2|HBA_HUMAN', description="Hemoglobin subunit alpha, human")
```

`SeqRecord` 类包括下列属性:

- **.seq**

  – 序列自身（即 `Seq` 对象）。

- **.id**

  – 序列ID。通常类同于accession number。

- **.name**

  – 序列名/id 。可以是accession number, 也可是clone名（类似GenBank record中的LOCUS id）。

- **.description**

  – 序列描述。

- **.letter_annotations**

  – 对照序列的每个字母逐字注释（per-letter-annotations），以信息名为键（keys），信息内容为值（value）所构成的字典。值与序列等长，用Python列表、元组或字符串表示。`.letter_annotations`可用于质量分数或二级结构信息 (如 Stockholm/PFAM 比对文件)等数据的存储。

- **.annotations**

  – 用于储存附加信息的字典。信息名为键（keys），信息内容为值（value）。用于保存序列的零散信息（如unstructured information）。

- **.features**

  – `SeqFeature` 对象列表，储存序列的结构化信息（structured information），如：基因位置, 蛋白结构域。

- **.dbxrefs**

  – 储存数据库交叉引用信息（cross-references）的字符串列表。

#### `SeqIO`模块

`Biopython`的`SeqIO`模块提供了多种常用文件格式的解析器。这些解析器从一个输入文件（从本地或数据库中）提取信息，而后自动转换成`SeqRecord`对象。`SeqIO`模块也提供了一种方法把`SeqRecord`对象写入到格式化的文件中。

##### 解析文件

序列文件的解析有两种方式：`SeqIO.parse()`和`SeqIO.read()`：这两种方法都有两个必须的参数和一个可选参数：

1. 输入文件，它指定从哪里读取数据；
2. 数据格式，如`fasta`或`genbank`，完整的支持格式列表参见：http:// biopython.org/wiki/SeqIO
3. 参数指定序列数据的字符集（可选）。

`SeqIO.parse()`和`SeqIO.read()`的区别在于：`SeqIO.parse()`返回的是一个迭代器，即从输入文件中产生几个`SeqRecord`对象，你可以用`for`或`while`循环来遍历他们；而当文件中只包含一条记录时，就像上面的例子，则必须用`SeqIO.read()`。换言之，`SeqIO.parse()`能处理输入文件中任意数目的记录，而`SeqIO.read()`只能处理一条记录的文件，后者会先检查文件中是否只有一条记录，否则将会产生错误。

###### 解析大文件

对于大量记录，可以使用`SeqIO.index()`方法，它需要两个参数：记录文件和文件格式。`SeqIO.index()`方法返回一个字典式对象，用它可以访问所有记录而不用把他们都读取到内存中。字典的键是记录的ID，值包含整个记录，后者也能用属性来访问，如`id`，`description`等。需要注意的是，这些类字典对象是只读的，也就是说创建后不能删除或插入。

###### 写文件

`SeqIO.write()`方法可将一个或多个`SeqRecord`对象写入指定格式的文件中。用这个方法需要三个参数：一个或多个`SeqRecord`对象，输出文件以及输出的格式。


## 示例

### 用`SeqIO`模块来解析一个多序列FASTA文件

代码会输出这些序列的标识、序列和长度。

```python
from Bio import SeqIO
fasta_file = open("Uniprot.fasta","r")
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    print(seq_record.id) 
    print(repr(seq_record.seq)) 
    print(len(seq_record)) 
fasta_file.close()
```

### 用`SeqIO`模块来解析文件并将其内容储存到列表或字典中

下面的代码将文件解析并储存到列表中，输出为第一条记录的标识和序列。

```python
from Bio import SeqIO
uniprot_iterator = SeqIO.parse("Uniprot.fasta", "fasta")
records = list(uniprot_iterator)
print(records[0].id)
print(records[0].seq)
```

此外可以使用字典，键是记录的ID，值包含记录的信息：

```python
from Bio import SeqIO
uniprot_iterator = SeqIO.parse("Uniprot.fasta", "fasta")
records = SeqIO.index("Uniprot.fasta","fasta")
print(len(records['sp|P03372|ESR1_HUMAN'].seq))
```

### 序列文件格式的转换

可以使用`SeqIO.parse()`和`SeqIO.write()`来转换序列格式，下面的脚本把一个`Genbank`文件转换成`FASTA`文件：

```python
from Bio import SeqIO
genbank_file = open ("AY810830.gb", "r")
output_file = open("AY810830.fasta", "w")
records = SeqIO.parse(genbank_file, "genbank")
SeqIO.write(records, output_file, "fasta")
output_file.close()
```
