# bioconductor注释：probe+cdf芯片

不同芯片有不同的包，bioconductor包可以做芯片ID转换。

## 用affy包读取affymetix的基因表达芯片数据-CEL格式数据

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180405/31fBa3dmjC.png?imageslim)

`affy`包是R语言的bioconductor系列包的一个，就一个功能，读取affymetix的基因表达芯片数据-CEL格式数据，处理成表达矩阵。但是，这个`affy`包支持的芯片平台是有限的，一般是hgu 95系列和133系列~~

首先，使用`GEOquery`包下载GEO数据库里的CEL文件，比如GSE1428，测了10 young (19-25 years old) 和12 older (70-80 years old) male的样品，下载的原始数据为79.4M。

```R
source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
getGEOSuppFiles("GSE1428")
```

GSE1428_RAW.tar 文件是压缩打包的 CEL（Affymetrix array 数据原始格式）文件。先解包数据，再解压数据。

```R
untar("GSE1428/GSE1428_RAW.tar", exdir="test")
cels <- list.files("test", pattern = "[gz]")
sapply(paste("test", cels, sep="/"), gunzip)
```

之后就可以用`affy`包来读取它了~

```R
biocLite("affy")
library(affy)
dir_cels='test'
affy_data = ReadAffy(celfile.path=dir_cels)
eset.mas5 = mas5(affy_data)
exprSet.nologs = exprs(eset.mas5)
exprSet = log(exprSet.nologs, 2)  #transform to Log_2 if needed
```

读取的过程较长，也可以选择rma函数而不是mas5函数对表达数据进行normalization。

```R
done.
22283 ids to be processed
|                    |
|####################|
```

读取之后的表达矩阵如图所示：

![mark](_v_images/_1522849047_1214.png)

## 用oligo包来读取affymetix的基因表达芯片数据-CEL格式数据

上一节中，`affy`包处理的芯片平台是有限的，一般是hgu 95系列和133系列，但无法读取和分析一些新版Affy芯片。而`oligo`包的处理方法以解决这些问题。而且，除了用于Affy芯片处理外，`oligo`包还可处理NimbleGen芯片。

同样地，首先下载并解压数据，以GSE48452为例，下载的原始数据为318.4 M。

```R
getGEOSuppFiles("GSE48452")
untar("GSE48452/GSE48452_RAW.tar", exdir="test2")
cels <- list.files("test2", pattern = "[gz]")
sapply(paste("test2", cels, sep="/"), gunzip)
```

使用`oligo`包：

```R
biocLite("oligo")
library(oligo)
geneCELs=list.celfiles('./test2',listGzipped=T,full.name=T)
# 读取cel文件
affyGeneFS <- read.celfiles(geneCELs)
# 这一步是normalization，会比较耗时
eset <- rma(affyGeneFS)
# 用exprs函数可以提取其中的表达量矩阵
data.exprs <- exprs(eset)
```

`read.celfiles`函数读取数据的过程中会检查系统中是否已安装相应的芯片注释软件包，没有安装的话将有警告。如果在读取的过程中没有自动安装相应注释软件包，请手动安装后再重新读取数据。


[使用oligo软件包处理芯片数据](https://blog.csdn.net/u014801157/article/details/66974577)

## 用lumi包来处理illumina的bead系列表达芯片

它封装好了一个函数，lumiExpresso可以直接处理LumiBatch对象，这个函数结合了,N,T,B,Q(normalization,transformation,backgroud correction,qulity control)四个步骤，其中Q这个步骤又包括8种统计学图片。在该包的文章有详细说明：http://bioinformatics.oxfordjournals.org/content/24/13/1547.full 


而 LumiBatch 对象是通过 lumiR.batch 读取的芯片文件被Illumina Bead Studio toolkit 处理的结果，也就是通常我们从公司或者GEO下载的数据( level 3 的 process data)

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180405/L3beH6BI29.png?imageslim)

以这个文件为例：ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30669/suppl/GSE30669_HEK_Sample_Probe_Profile.txt.gz ，就可以直接用`lumi`包的`lumiR.batch`函数读取文件成为`LumiBatch`对象，然后被`lumiExpresso`函数直接处理，然后被`exprs`函数提取表达矩阵。

```R
biocLite("lumi")
library(lumi)
## 首先是从illumina的芯片结果文件，用R的lumi包来获取表达矩阵。
setwd('./lumitest')
fileName <- 'GSE30669_HEK_Sample_Probe_Profile.txt' # Not Run
x.lumi <- lumiR.batch(fileName) ##, sampleInfoFile='sampleInfo.txt')
pData(phenoData(x.lumi))
## Do all the default preprocessing in one step
lumi.N.Q <- lumiExpresso(x.lumi)
### retrieve normalized data
dataMatrix <- exprs(lumi.N.Q)

## 下面是从GEO里面下载表达矩阵
rm(list=ls())
library(GEOquery)
library(limma)
GSE30669 <- getGEO('GSE30669', destdir=".",getGPL = F)
exprSet=exprs(GSE30669[[1]])
GSE30669[[1]]
pdata=pData(GSE30669[[1]])
exprSet=exprs(GSE30669[[1]])
## 很明显可以看到前面得到的dataMatrix 和后面得到的 exprSet 都是我们想要的表达矩阵
```






