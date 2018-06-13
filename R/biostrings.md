## FASTA文件中的DNA/氨基酸序列

### FastA格式

是一种基于文本用于表示核苷酸序列或氨基酸序列的格式。在这种格式中碱基对或氨基酸用单个字母来编码，且允许在序列前添加序列名及注释。序列文件的第一行是由大于号">"或分号";"打头的任意文字说明（习惯常用">"作为起始），用于序列标记。从第二行开始为序列本身，只允许使用既定的核苷酸或氨基酸编码符号。通常核苷酸符号大小写均可，而氨基酸常用大写字母。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180401/l0KE1fmc70.png?imageslim)

`biostrings`包可以用于处理生物序列比如DNA，RNA和蛋白质。还可以用于将fasta序列读写输出。此外，该包具有用于模式匹配（短读取对齐）的功能，以及实现Smith-Waterman局部比对和经典序列比对中使用的Needleman-Wunsch全局比对的成对比对函数。

`biostrings`包的功能十分强大，以下内容摘自 Biostrings Quick Overview

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180327/gJ65IK73Ba.png?imageslim)
![mark](http://oo3g995ih.bkt.clouddn.com/blog/180327/aDgF8k4aic.png?imageslim)

### 下载并加载biostrings

```R
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library(Biostrings)
```

### 表示序列

```R
# 直接输入序列
dna1 <- DNAString("ACGT-N")
dna1
dna2 <- DNAStringSet(c("ACGT", "GTCA", "GCTA"))
dna2
# 从fasta读入
dna3 = readDNAStringSet("./test/sequence.fasta",format="fasta")
dna3
# 详细参数
readBStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
readDNAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
readRNAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
readAAStringSet(filepath, format="fasta",
               nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
# nrec:最多从文件读入多少记录，负数则表示该项被忽略；
# skip:一个非负整数，表示在开始读入记录之前要跳过的记录数；
# seek.first.rec:TRUE或者FALSE（默认）。 如果为TRUE，则一定是以“>”（对于FASTA）或“@”（对于FASTQ）作为第一行。 如果找不到这样的行，就会出现错误。值得一提的是，TRUE还可以用于解析带有FASTA数据的GFF3文件。
# use.names:返回的向量是否应该被命名，对于FASTA，名称取自描述行。对于FASTQ，它们取自记录序列ID。 删除名称可以帮助减少内存占用，例如：对于包含数百万个记录的FASTQ文件。
```

### Biostrings定义的常量

包括DNA、RNA、AA、密码子的符号，载入Biostrings包后这些常量就可以直接使用：

```R
> DNA_BASES
[1] "A" "C" "G" "T"
> DNA_ALPHABET
 [1] "A" "C" "G" "T" "M" "R" "W" "S" "Y" "K" "V" "H" "D" "B" "N" "-" "+" "."
> RNA_BASES
[1] "A" "C" "G" "U"
> RNA_ALPHABET
 [1] "A" "C" "G" "U" "M" "R" "W" "S" "Y" "K" "V" "H" "D" "B" "N" "-" "+" "."
> RNA_GENETIC_CODE
UUU UUC UUA UUG UCU UCC UCA UCG UAU UAC UAA UAG UGU UGC UGA UGG CUU CUC CUA CUG CCU CCC CCA CCG CAU CAC CAA CAG CGU 
"F" "F" "L" "L" "S" "S" "S" "S" "Y" "Y" "*" "*" "C" "C" "*" "W" "L" "L" "L" "L" "P" "P" "P" "P" "H" "H" "Q" "Q" "R" 
CGC CGA CGG AUU AUC AUA AUG ACU ACC ACA ACG AAU AAC AAA AAG AGU AGC AGA AGG GUU GUC GUA GUG GCU GCC GCA GCG GAU GAC 
"R" "R" "R" "I" "I" "I" "M" "T" "T" "T" "T" "N" "N" "K" "K" "S" "S" "R" "R" "V" "V" "V" "V" "A" "A" "A" "A" "D" "D" 
GAA GAG GGU GGC GGA GGG 
"E" "E" "G" "G" "G" "G" 
attr(,"alt_init_codons")
[1] "UUG" "CUG"
> AA_ALPHABET
 [1] "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T" "W" "Y" "V" "U" "O" "B" "J" "Z" "X" "*" "-"
[29] "+" "."
> AMINO_ACID_CODE
    A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y 
"Ala" "Arg" "Asn" "Asp" "Cys" "Gln" "Glu" "Gly" "His" "Ile" "Leu" "Lys" "Met" "Phe" "Pro" "Ser" "Thr" "Trp" "Tyr" 
    V     U     O     B     J     Z     X 
"Val" "Sec" "Pyl" "Asx" "Xle" "Glx" "Xaa" 
> GENETIC_CODE
TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT 
"F" "F" "L" "L" "S" "S" "S" "S" "Y" "Y" "*" "*" "C" "C" "*" "W" "L" "L" "L" "L" "P" "P" "P" "P" "H" "H" "Q" "Q" "R" 
CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC 
"R" "R" "R" "I" "I" "I" "M" "T" "T" "T" "T" "N" "N" "K" "K" "S" "S" "R" "R" "V" "V" "V" "V" "A" "A" "A" "A" "D" "D" 
GAA GAG GGT GGC GGA GGG 
"E" "E" "G" "G" "G" "G" 
attr(,"alt_init_codons")
[1] "TTG" "CTG"
> IUPAC_CODE_MAP
     A      C      G      T      M      R      W      S      Y      K      V      H      D      B      N 
   "A"    "C"    "G"    "T"   "AC"   "AG"   "AT"   "CG"   "CT"   "GT"  "ACG"  "ACT"  "AGT"  "CGT" "ACGT" 
```

### 获取序列基本信息

```R
# 检索一个子序列（类似于标准的R函数substr）
dna1[2:4]
dna2[2:3]
# [[可以从DNAStringSet取出一个元素作为一个DNAString，跟list中的用法一致
dna2[[2]]
# DNAStringSet可以设置名字
names(dna2) <- paste0("seq", 1:3)
dna2
# rev运用在DNAStringSet只是反向元素的顺序，而在 DNAString是整个string的反转。reverse是反转序列
dna1
rev(dna1)
reverse(dna1)
dna2
rev(dna2)
reverse(dna2)
# width查看序列长度
width(dna1)
width(dna2)
# letterFrequency用于统计指定字符的频率或比例
letterFrequency(dna1, DNA_BASES)
letterFrequency(dna1, DNA_ALPHABET)
letterFrequency(dna1, DNA_BASES, as.prob = TRUE)
# 计算GC含量
letterFrequency(dna1, "GC", as.prob = TRUE)
# 分别统计2核苷酸组合和3核苷酸组合
dinucleotideFrequency(dna1)
trinucleotideFrequency(dna1)
```

### 此外

`treeio`包里的`read.fasta`，也可以读入序列，它有个好处是你直接`as.character(aa[[i]])`就是个字母的向量，方便比较。

```R
require(treeio)
fasta <- system.file("./test/sequence.fasta", package="ggtree")
aa <- read.fasta(fasta)
```

### 利用biostrings包来实现：生信编程直播第二题-hg19基因组序列的一些探究

- 使用测试数据：

```
>chr_1
ATCGTCGaaAATGAANccNNttGTA
AGGTCTNAAccAAttGggG
>chr_2
ATCGAATGATCGANNNGccTA
AGGTCTNAAAAGG
>chr_3
ATCGTCGANNNGTAATggGA
AGGTCTNAAAAGG
>chr_4
ATCGTCaaaGANNAATGANGgggTA
```

- 代码如下：

```R
dna3 = readDNAStringSet("./test/sequence.fasta",format="fasta")
DNA_count <- letterFrequency(dna3, DNA_BASES)
GC_percentage <- letterFrequency(dna3, "GC", as.prob = TRUE)
Length <- width(dna3)
results <- cbind(DNA_count,GC_percentage, Length)
results
```

```
      A C  G  T       G|C Length
[1,] 13 7 10 10 0.3863636     44
[2,] 11 5  8  6 0.3823529     34
[3,] 10 3 10  6 0.3939394     33
[4,]  9 2  7  4 0.3600000     25
```

## FASTQ文件中的Reads

### FASTQ格式

FASTQ文件中每个序列通常有四行：
第一行：必须以“@”开头，后面跟着唯一的序列ID标识符，然后跟着可选的序列描述内容，标识符与描述内容用空格分开；
第二行：序列字符（核酸为[AGCTN]+，蛋白为氨基酸字符）；
第三行：必须以“+”开头，后面跟着可选的ID标识符和可选的描述内容，如果“+”后面有内容，该内容必须与第一行“@”后的内容相同；
第四行：碱基质量字符，每个字符对应第二行相应位置碱基或氨基酸的质量，该字符可以按一定规则转换为碱基质量得分，碱基质量得分可以反映该碱基的错误率。这一行的字符数与第二行中的字符数必须相同。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180330/8l16g91aBH.png?imageslim)

具体信息可参见：浅谈FastQ和FastA格式

`ShortRead`包用于处理fastq文件。以下是主要功能：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180401/mGEfDl9A5H.png?imageslim)

### 下载并加载ShortRead

```R
source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
library(ShortRead)
```

### 处理FASTQ文件

```R
> fq <- readFastq("./test/sequence.fastq")
> fq
class: ShortReadQ
length: 1 reads; width: 59 cycles
> head(sread(fq),1)
  A DNAStringSet instance of length 1
    width seq
[1]    59 GATTTGGGGTTCAAAGC...TGTTCAACTCACAGTT
> head(quality(fq),1)
class: FastqQuality
quality:
  A BStringSet instance of length 1
    width seq
[1]    59 !''*((((***+))%%%...CF>>>>>>CCCCCCC6
> encoding(quality(fq))
 !  "  #  $  %  &  '  (  )  *  +  ,  -  .  / 
 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 
 0  1  2  3  4  5  6  7  8  9  :  ;  <  =  > 
15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
 ?  @  A  B  C  D  E  F  G  H  I  J  K  L  M 
30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 
 N  O  P  Q  R  S  T  U  V  W  X  Y  Z  [ \\ 
45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 
 ]  ^  _  `  a  b  c  d  e  f  g  h  i  j  k 
60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 
 l  m  n  o  p  q  r  s  t  u  v  w  x  y  z 
75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 
 {  |  }  ~ 
90 91 92 93
```

同样地，用`ShortRead`也可以读入fasta文件

```R
> fa <- readFasta("./test/sequence.fasta")
> fa
class: ShortRead
length: 4 reads; width: 25..44 cycles
> head(sread(fa),1)
  A DNAStringSet instance of length 1
    width seq
[1]    44 ATCGTCGAAAATGAANC...TCTNAACCAATTGGGG
```


## 参考资料

1. Introduction to Bioconductor for Sequence Data

http://www.bioconductor.org/help/workflows/sequencing/

2. Biostrings的bioconductor主页

http://www.bioconductor.org/packages/release/bioc/html/Biostrings.html

3. R/BioC序列处理之二：Biostrings序列的基本操作

https://www.plob.org/article/7568.html

4. Biostrings in R

https://web.stanford.edu/class/bios221/labs/biostrings/lab_1_biostrings.html

5. ShortRead的bioconductor主页

http://www.bioconductor.org/packages/release/bioc/html/ShortRead.html


---

aaa