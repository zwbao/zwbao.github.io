# phyloseq | 用 R 分析微生物组数据及可视化

## 导入数据

phyloseq 包是一个集OTU 数据导入，存储，分析和图形可视化于一体的工具。它不但利用了 R 中许多经典的工具进行生态学和系统发育分析（例如：vegan，ade4，ape， picante），同时还结合 ggplot2 以轻松生成发表级别的可视化结果。phyloseq 使用的S4类将一个研究所有相关的测序数据及元数据存储为单个对象，从而更容易共享数据并重复结果。

- GitHub 地址：https://github.com/joey711/phyloseq

### 安装

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq", version = "3.8")
```

### 数据导入

一个 phyloseq 类，通常由以下几个部分组成：

- `otu_table` ：一个数字矩阵 `matrix`，包含了 OTU 在每个样本中的丰度信息；
- `sample_data` ：一个`data.frame`，包含了所有样本的表型信息，行名必须匹配`otu_table` 中的样本名；
- `tax_table` ：一个字符矩阵 `matrix`，包含了 OTU 的物种信息，行名必须匹配`otu_table` 中的 OTU 名。

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g0pt78c5eij234l1jpu0x.jpg)

在下面的例子中，我用 R 随机创造了一些模拟数据，要是有真实数据的话，也可以直接导入。

```R
# Create a pretend OTU table that you read from a file, called otumat
otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
# Create a pretend taxonomy table
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
```

看一下我们导入的 OTU 丰度信息和 OTU 分类信息：

```R
otumat
```

```R
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       96      50      36      35      59      80      83      63
## OTU2       52      67      39      39      37      57      20      15
## OTU3       94      18      15      11      14      75       1      12
## OTU4       27      88      98     100      59      27      30      30
## OTU5       26      66      93      85      41      30     100      41
## OTU6       17      16      97      86      18      25      94      31
## OTU7       63      19      16      43      89      25      17      63
## OTU8       31      92      22      14      58       1      45       2
## OTU9      100      33      19      77      43       1      14      69
## OTU10      13      35      80      43      34      45      24      47
##       Sample9 Sample10
## OTU1       38       35
## OTU2       64       94
## OTU3       42       58
## OTU4       94       78
## OTU5       92      100
## OTU6       62       37
## OTU7       15       82
## OTU8       25       35
## OTU9       42       18
## OTU10      71       72
```

```R
taxmat
```

```R
##       Domain Phylum Class Order Family Genus Species
## OTU1  "x"    "d"    "q"   "v"   "l"    "k"   "i"    
## OTU2  "a"    "d"    "x"   "a"   "k"    "o"   "r"    
## OTU3  "h"    "a"    "h"   "c"   "d"    "j"   "k"    
## OTU4  "t"    "f"    "j"   "e"   "n"    "y"   "o"    
## OTU5  "o"    "q"    "s"   "w"   "d"    "y"   "j"    
## OTU6  "e"    "r"    "p"   "k"   "b"    "v"   "t"    
## OTU7  "m"    "l"    "y"   "u"   "b"    "y"   "q"    
## OTU8  "d"    "o"    "w"   "g"   "p"    "w"   "v"    
## OTU9  "f"    "o"    "a"   "n"   "l"    "u"   "e"    
## OTU10 "h"    "r"    "d"   "j"   "u"    "f"   "a"
```

接下来，我们需要将他们组合成一个 phyloseq 对象：

```R
library("phyloseq")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
```

```R
## OTU Table:          [10 taxa and 10 samples]
##                      taxa are rows
##       Sample1 Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8
## OTU1       96      50      36      35      59      80      83      63
## OTU2       52      67      39      39      37      57      20      15
## OTU3       94      18      15      11      14      75       1      12
## OTU4       27      88      98     100      59      27      30      30
## OTU5       26      66      93      85      41      30     100      41
## OTU6       17      16      97      86      18      25      94      31
## OTU7       63      19      16      43      89      25      17      63
## OTU8       31      92      22      14      58       1      45       2
## OTU9      100      33      19      77      43       1      14      69
## OTU10      13      35      80      43      34      45      24      47
##       Sample9 Sample10
## OTU1       38       35
## OTU2       64       94
## OTU3       42       58
## OTU4       94       78
## OTU5       92      100
## OTU6       62       37
## OTU7       15       82
## OTU8       25       35
## OTU9       42       18
## OTU10      71       72
```

```R
TAX
```

```R
## Taxonomy Table:     [10 taxa by 7 taxonomic ranks]:
##       Domain Phylum Class Order Family Genus Species
## OTU1  "x"    "d"    "q"   "v"   "l"    "k"   "i"    
## OTU2  "a"    "d"    "x"   "a"   "k"    "o"   "r"    
## OTU3  "h"    "a"    "h"   "c"   "d"    "j"   "k"    
## OTU4  "t"    "f"    "j"   "e"   "n"    "y"   "o"    
## OTU5  "o"    "q"    "s"   "w"   "d"    "y"   "j"    
## OTU6  "e"    "r"    "p"   "k"   "b"    "v"   "t"    
## OTU7  "m"    "l"    "y"   "u"   "b"    "y"   "q"    
## OTU8  "d"    "o"    "w"   "g"   "p"    "w"   "v"    
## OTU9  "f"    "o"    "a"   "n"   "l"    "u"   "e"    
## OTU10 "h"    "r"    "d"   "j"   "u"    "f"   "a"
```

```R
physeq = phyloseq(OTU, TAX)
physeq
```

```R
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10 taxa and 10 samples ]
## tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]
```

```R
plot_bar(physeq, fill = "Family")
```

![](https://joey711.github.io/phyloseq/import-data_files/figure-html/unnamed-chunk-6-1.png)

下面导入样本的表型信息：

```R
# Create random sample data
sampledata = sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  Depth = sample(50:1000, size=nsamples(physeq), replace=TRUE),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
sampledata
```

```R
##          Location Depth
## Sample1         D   337
## Sample2         B    74
## Sample3         D    68
## Sample4         C   397
## Sample5         B   142
## Sample6         D   970
## Sample7         D    69
## Sample8         C   253
## Sample9         A   497
## Sample10        D   237
```

使用`ape`包建立 OTU 系统发育树并导入 phyloseq 对象：

```R
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)
```

![](https://joey711.github.io/phyloseq/import-data_files/figure-html/unnamed-chunk-8-1.png)

现在，我们有了`otu_table`，`sample_data`，`tax_table`，`phy_tree`这四类数据，以下两种方法都可以将他们合并为一个 phyloseq 对象：1. 使用`merge_phyloseq`函数在之前创建的`physeq`对象中加入`sample_data`和`phy_tree`数据；2. 使用`physeq`函数重新创建一个`physeq`对象。

```R
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
physeq1
```

```R
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 2 sample variables ]
## tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 10 tips and 9 internal nodes ]
```

使用`physeq`函数重新创建一个`physeq`对象：

```R
physeq2 = phyloseq(OTU, TAX, sampledata, random_tree)
physeq2
```

```R
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 2 sample variables ]
## tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 10 tips and 9 internal nodes ]
```

看看用这两种方法创建的`physeq`对象是否一样：

```R
identical(physeq1, physeq2)
```

```R
## [1] TRUE
```

尝试对这些数据进行一些可视化：

```R
plot_tree(physeq1, color="Location", label.tips="taxa_names", ladderize="left", plot.margin=0.3)
```

![](https://joey711.github.io/phyloseq/import-data_files/figure-html/treeplot-1.png)

```R
plot_tree(physeq1, color="Depth", shape="Location", label.tips="taxa_names", ladderize="right", plot.margin=0.3)
```

![](https://joey711.github.io/phyloseq/import-data_files/figure-html/treeplot-2.png)

再做一些热图：

```R
plot_heatmap(physeq1)
```

![](https://joey711.github.io/phyloseq/import-data_files/figure-html/heatmap-random-1.png)

```R
plot_heatmap(physeq1, taxa.label="Phylum")
```

![](https://joey711.github.io/phyloseq/import-data_files/figure-html/heatmap-random-2.png)

**Reference**

- https://joey711.github.io/phyloseq/import-data.html

## 分析实战


首先载入`phyloseq`和`ggplot`：

```R
library("phyloseq")
library("ggplot2")
# 设置 ggplot2 主题
theme_set(theme_bw())
```

在这次的教程中以`phyloseq`自带的`GlobalPatterns`数据集和`enterotype`数据集作为例子。`GlobalPatterns`数据来自一篇2011年的 PNAS 文章 （Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample），比较了 25 个环境样本和 3 个已知的“模拟微生物群落”，包含了 9 种样本类型。`enterotype`数据集即人类肠型数据集来自一篇2011 年的 Nature 文章（Enterotypes of the human gut microbiome），该文章使用宏基因组鸟枪法测序比较了22名受试者的粪便微生物群落，作者进一步将这些微生物群落与来自其他研究的受试者的粪便群落进行比较。

### 物种多样性分析

导入`GlobalPatterns`数据集：

```R
data(GlobalPatterns)
```

首先对这些数据做一些预处理：

```R
# prune OTUs that are not present in at least one sample
GP <- prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)
# Define a human-associated versus non-human categorical variable:
human <- get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
# Add new human variable to sample data:
sample_data(GP)$human <- factor(human)
```

我们可以使用`plot_richness`函数轻松创建一个复杂的图形，用于比较`GlobalPatterns`数据集中不同环境类型的样本的物种多样性。

```R
alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
(p <- plot_richness(GP, "human", "SampleType", measures=alpha_meas))
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g0vl6upuilj21w310cn0d.jpg)

还可以往上面加`ggplot`图层：

```R
p + geom_boxplot(data=p$data, aes(x=human, y=value, color=NULL), alpha=0.1)
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g0vl8js8khj21w310c77q.jpg)

以上就是`GlobalPatterns`数据集中样本的 Alpha 多样性分析。每幅图都展示了不同的物种多样性指数，不同的颜色表示不同的样本类型。在每个小组中，将样本进一步分为与人类相关的（TRUE）或不相关的（FALSE），并且在这两个组的顶部叠加箱图，说明这些与人类相关样本不如环境样本物种丰富。

### 绘制进化树

`phyloseq`还可以很容易地绘制带注释的系统发育树，并加入特定分类的信息，以及观察到的个体数量。

这一节中，我们先取只包含`Chlamydiae`门的数据：

```R
GP.chl <- subset_taxa(GP, Phylum=="Chlamydiae")
```

现在我们将创建这个`GlobalPatterns`数据集的子集的进化树并用`SampleType`变量进行着色，该变量表示微生物组样本的来源。以下命令还可以标记每个分类中每个样本（如果有的话）中观察到的菌群数量，菌群丰度越高，形状就越大。

```R
plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="Abundance")
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g0vlookuilj226m10eaex.jpg)

### 条形图

导入`enterotype`数据集：

```
data(enterotype)
```

在`enterotype`数据集中，可用的数据是菌群的相对丰度，而不是原始的数量，所以就不能用来做物种多样性分析。而且 OTU 仅有属级别的注释信息。 

我们先从一个简单的秩丰度条形图开始，计算数据集中每个OTU的累积分数丰度。在这个条形图中，我们通过样本总数（280）进一步标准化。

```R
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
N <- 30
barplot(sort(taxa_sums(enterotype), TRUE)[1:N]/nsamples(enterotype), las=2)
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g0vm4hfy7mj211k0pymxr.jpg)

上面的例子，我们用到了`R`最基础的`barplot`函数和`phyloseq`提供的`taxa_sums`和`nsamples`函数。你可以看到大约在排名第十的菌属后菌群的相对丰度出现了剧烈下滑。因此我们可以尝试在最多的十个菌属中来分析肠型数据。

```R
TopNOTUs <- names(sort(taxa_sums(enterotype), TRUE)[1:10]) 
ent10   <- prune_taxa(TopNOTUs, enterotype)
print(ent10)
```

```R
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 10 taxa and 280 samples ]
## sample_data() Sample Data:       [ 280 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 10 taxa by 1 taxonomic ranks ]
```

这个数据集中共有 280 个样本。接着，看看原始的`metadata`中包含了哪些表型信息：

```R
sample_variables(ent10)
```

```R
## [1] "Enterotype"     "Sample_ID"      "SeqTech"        "SampleID"      
## [5] "Project"        "Nationality"    "Gender"         "Age"           
## [9] "ClinicalStatus"
```

在下面的例子中，我们根据每个排名前十的菌属来分面（`facet`），在每个子图中以三种不同的测序方法分成三组，并以颜色表示每个菌属来源的样品的肠型标签。可以发现，不同的测序技术对检测到哪个属和菌属的丰度水平具有很大的影响。

```R
plot_bar(ent10, "SeqTech", fill="Enterotype", facet_grid=~Genus)
```

![](http://ww1.sinaimg.cn/mw690/c5d7b0ebly1g0vmouslrij21g40v6acs.jpg)

从上图还可以看出，肠型1 中有较高丰度的 Bacteroides， 肠型2 中有较高丰度的 Prevotella；而对于肠型3 可以观察到高丰度的 Blautia 仅在454-焦磷酸测序的数据中出现，在 Illumina 或 Sanger 的数据中丰度都很低，说明高丰度的 Blautia 可能是引入的误差。

### 热图

下面以`GlobalPatterns`数据集的子集（仅包含`Crenarchaeota`门的数据）来展示热图的功能：

```R
data("GlobalPatterns")
gpac <- subset_taxa(GlobalPatterns, Phylum=="Crenarchaeota")
(p <- plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family"))
```

![](http://ww1.sinaimg.cn/mw690/c5d7b0ebly1g0vn18z2ejj21as10igoo.jpg)

如果需要自定义`X`轴和`Y`轴的标签，只要这样就好啦：

```R
p$scales$scales[[1]]$name <- "My X-Axis"
p$scales$scales[[2]]$name <- "My Y-Axis"
print(p)
```

除此之外，你还可以很方便地选择不同的距离矩阵（例子中使用了`Bray Curtis`距离矩阵）来进行排序，是不是很简单~

### 微生物网络图

这一节中，我们使用`enterotype`数据集，来创建网络图：

```R
data(enterotype)
plot_net(enterotype, maxdist=0.4, color="SeqTech", shape="Enterotype")
```

![](http://ww1.sinaimg.cn/mw690/c5d7b0ebly1g0vn9hjxlzj21jg0xdq9e.jpg)

## 排序分析

排序分析（Ordination analysis）是探索复杂的系统发育测序数据的得力工具。简而言之，排序（ordination）的过程就是在一个可视化的低维空间重新排列这些样方，使得样方之间的距离最大程度地反映出平面散点图内样方之间的关系信息。`phyloseq`分别通过`ordinate()`和`plot_ordination()`函数，来进行排序分析及其可视化。虽然排序分析有很多方法和参数，但最基本的代码就是这样：

```R
my.physeq <- import("Biom", BIOMfilename="myBiomFile.biom")
my.ord    <- ordinate(my.physeq)
plot_ordination(my.physeq, my.ord, color="myFavoriteVarible")
```

### PCoA 主坐标分析

PCoA（principal co-ordinates analysis）主坐标分析法是一种与 PCA 类似的降维排序方法。PCoA与PCA的区别在于PCA是基于原始的物种组成矩阵所做的分析，仅仅比较的是物种丰度的不同，而 PCoA 首先根据不同的距离算法计算样品之间的距离，然后对距离矩阵进行处理，使图中点间的距离正好等于原来的差异数据，实现定性数据的定量转换。

下面以 **Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample** 这篇文章的 **Figure 5** 为例，作者展示了 unweighted-UniFrac 距离矩阵 PCoA 的前三个轴。

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g14fwcrakkj20gt0graea.jpg)

下面将在完整的数据集上再现 unweighted-UniFrac 距离矩阵的计算。由于 OTU 数量众多，所以可能需要很长的时间，建议对大型数据集进行并行化运算，有关并行化的详细信息，请参阅`UniFrac()`的文档及示例。

```R
data(GlobalPatterns)
# GPUF <- UniFrac(GlobalPatterns)
# Load the pre-computed distance matrix, GPUF
load(system.file("doc", "Unweighted_UniFrac.RData", package="phyloseq"))
# Calculate the PCoA on this distance matrix, GPUF
GloPa.pcoa = ordinate(GlobalPatterns, method="PCoA", distance=GPUF)
```

在查看结果之前，我们可以用碎石图来查看排名前几的轴的贡献率：

```R
plot_scree(GloPa.pcoa, "Scree plot for Global Patterns, UniFrac/PCoA")
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g14gg4pltxj20v80kcdg3.jpg)

通过碎石图，可以看出前三个轴代表距离总变化的43％，第四轴代表另外9％，因此也可能需要进行探索。碎石图是任何排序方法的重要工具，因为轴的相对贡献率在不同的数据集之间可能变化很大。

下面开始重复 **Figure 5**，但我们将原来三维的图，转化为两张图：

```R
(p12 <- plot_ordination(GlobalPatterns, GloPa.pcoa, "samples", color="SampleType") + 
  geom_point(size=5) + geom_path() + scale_colour_hue(guide = FALSE) )
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g14go1zni1j21b810ign3.jpg)

```R
(p13 <- plot_ordination(GlobalPatterns, GloPa.pcoa, "samples", axes=c(1, 3),
  color="SampleType") + geom_line() + geom_point(size=5) )
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g14gpvutfij21ar10ijtp.jpg)

### NMDS 非度量多维尺度分析

非度量多维尺度分析是一种**将多维空间的研究对象**（样本或变量）**简化到低维空间进行定位、分析和归类，同时又保留对象间原始关系的数据分析方法**。适用于无法获得研究对象间精确的相似性或相异性数据，仅能得到他们之间等级关系数据的情形。其基本特征是将对象间的相似性或相异性数据看成点间距离的单调函数，在保持原始数据次序关系的基础上，用新的相同次序的数据列替换原始数据进行度量型多维尺度分析。

其特点是根据样品中包含的物种信息，以点的形式反映在多维空间上，而对不同样品间的差异程度，则是通过点与点间的距离体现的，最终获得样品的空间定位点图。

```R
# (Re)load UniFrac distance matrix and GlobalPatterns data
data(GlobalPatterns)
load(system.file("doc", "Unweighted_UniFrac.RData", package="phyloseq"))
# perform NMDS, set to 2 axes
GP.NMDS <- ordinate(GlobalPatterns, "NMDS", GPUF)
```

```R
## Run 0 stress 0.1432774 
## Run 1 stress 0.1432625 
## ... New best solution
## ... Procrustes: rmse 0.003499558  max resid 0.01265589 
## Run 2 stress 0.2467057 
## Run 3 stress 0.2234687 
## Run 4 stress 0.1432774 
## ... Procrustes: rmse 0.003483686  max resid 0.01260292 
## Run 5 stress 0.1432774 
## ... Procrustes: rmse 0.003487121  max resid 0.01261713 
## Run 6 stress 0.1856306 
## Run 7 stress 0.1432625 
## ... New best solution
## ... Procrustes: rmse 1.405265e-05  max resid 3.151096e-05 
## ... Similar to previous best
## Run 8 stress 0.2574296 
## Run 9 stress 0.1809282 
## Run 10 stress 0.1809285 
## Run 11 stress 0.1432625 
## ... New best solution
## ... Procrustes: rmse 4.788746e-06  max resid 1.24343e-05 
## ... Similar to previous best
## Run 12 stress 0.1432774 
## ... Procrustes: rmse 0.003470474  max resid 0.01255052 
## Run 13 stress 0.1432625 
## ... Procrustes: rmse 1.152937e-05  max resid 3.31644e-05 
## ... Similar to previous best
## Run 14 stress 0.1432625 
## ... New best solution
## ... Procrustes: rmse 2.462092e-06  max resid 6.050872e-06 
## ... Similar to previous best
## Run 15 stress 0.1856315 
## Run 16 stress 0.1432625 
## ... Procrustes: rmse 9.826556e-06  max resid 3.260556e-05 
## ... Similar to previous best
## Run 17 stress 0.1809284 
## Run 18 stress 0.1856306 
## Run 19 stress 0.1432774 
## ... Procrustes: rmse 0.003490316  max resid 0.01261738 
## Run 20 stress 0.1432774 
## ... Procrustes: rmse 0.003471665  max resid 0.01255283 
## *** Solution reached
```

```R
(p <- plot_ordination(GlobalPatterns, GP.NMDS, "samples", color="SampleType") +
  geom_line() + geom_point(size=5) )
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g14gzm4g7mj21as10djth.jpg)

使用非度量多维尺度分析对`GlobalPatterns`数据集的样本之间的 unweighted-UniFrac 距离矩阵进行排序。每个样本点根据环境类型着色，如果它们属于同一类型，则由一条线连接。与`GlobalPatterns`文章中的图5进行比较。该图很好地显示了来自不同栖息地的微生物群落之间的相对差异。但是，它未能表明群落之间的差异。对于提供解释样本（或样本组）之间差异的群落信息的排序方法，我们使用对应分析（CA）。

### CA 对应分析

在这一节中，我们使用 Correspondence Analysis 的排序方法的各种功能继续探索`GlobalPatterns`数据集。

让我们首先进行对应分析并查看碎石图。碎石图表明`GlobalPatterns`数据集具有相当高的维度，前两个 CA 轴占总变异的不到17％，但随着轴数的增加，贡献度并没有急剧下降，每个轴仅比前一个下降一点点。

首先，为了节约运行时间，我们只取`GlobalPatterns`数据集的子集：

```R
data(GlobalPatterns)
# Take a subset of the GP dataset, top 200 species
topsp <- names(sort(taxa_sums(GlobalPatterns), TRUE)[1:200])
GP    <- prune_taxa(topsp, GlobalPatterns)
# Subset further to top 5 phyla, among the top 200 OTUs.
top5ph <- sort(tapply(taxa_sums(GP), tax_table(GP)[, "Phylum"], sum), decreasing=TRUE)[1:5]
GP     <- subset_taxa(GP, Phylum %in% names(top5ph))
# Define a human-associated versus non-human categorical variable:
human <- get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
# Add new human variable to sample data:
sample_data(GP)$human <- factor(human)
```

进行 CA 对应分析：

```R
# Now perform a unconstrained correspondence analysis
gpca  <- ordinate(GP, "CCA")
# Scree plot
plot_scree(gpca, "Scree Plot for Global Patterns Correspondence Analysis")
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g14hxgn7ssj215w0po74k.jpg)

```R
(p12 <- plot_ordination(GP, gpca, "samples", color="SampleType") + 
  geom_line() + geom_point(size=5) )
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g14hyr35ovj21as10iwgf.jpg)

```R
(p34 <- plot_ordination(GP, gpca, "samples", axes=c(3, 4), color="SampleType") + 
  geom_line() + geom_point(size=5) )
```

![](http://ww1.sinaimg.cn/large/c5d7b0ebly1g14i08jm0kj21as10itag.jpg)

第一张图说明粪便和模拟数据紧密地聚集在一起，远离 CA1 上的其他样本。在 CA2 上，皮肤和舌头样本也类似地分开。总而言之，前两个轴已经解释了人类相关“环境”与数据集中其他非人类环境的分离，以及舌头和皮肤样本与粪便的分离。

我们现在将进一步研究这种数据结构，使我们能够比较单个分类群在相同图形空间上的贡献度：

```R
p1  <- plot_ordination(GP, gpca, "species", color="Phylum")
(p1 <- ggplot(p1$data, p1$mapping) + geom_point(size=5, alpha=0.5) + 
  facet_wrap(~Phylum) +  scale_colour_hue(guide = FALSE) )
```

![](https://static.xmt.cn/8845d9c8deb24407b261f0ee1428142e.png)

我们还可以将物种点总结为 2D 密度图：

```R
(p3 <- ggplot(p1$data, p1$mapping) + geom_density2d() +
  facet_wrap(~Phylum) +  scale_colour_hue(guide = FALSE) )
```

![](https://static.xmt.cn/de33ce3e38104373b17b87d3a18d9091.png)

上面的图表展示了一些有用的模式和有趣的异常值，但如果我们想要具体了解这些门在多大程度上影响每个轴呢？下面的代码使用了箱线图来展示这个结果。

```R
library("reshape2")
# Melt the species-data.frame, DF, to facet each CA axis separately
mdf <- melt(p1$data[, c("CA1", "CA2", "Phylum", "Family", "Genus")], 
            id=c("Phylum", "Family", "Genus") )
# Select some special outliers for labelling
LF <- subset(mdf, variable=="CA2" & value < -1.0)
# build plot: boxplot summaries of each CA-axis, with labels
p <- ggplot(mdf, aes(Phylum, value, color=Phylum)) + 
  geom_boxplot() + 
  facet_wrap(~variable, 2) + 
  scale_colour_hue(guide = FALSE) +
  theme_bw() + 
  theme( axis.text.x = element_text(angle = -90, vjust = 0.5) )
# Add the text label layer, and render ggplot graphic
(p <- p + geom_text(data=subset(LF, !is.na(Family)),
  mapping = aes(Phylum, value+0.1, color=Phylum, label=Family), 
  vjust=0,
  size=2))
```

![](https://static.xmt.cn/ac6cc012e7ee4848bc3e207b6185d7b3.png)

通过这种方法，更容易看到相对于其他门的异常聚集的特定物种，例如上图中的拟杆菌（Bacteroidetes）属的 Prevotellaceae，在舌/皮肤样品中负 CA2 方向贡献最多。

另一种方法是直接展示物种丰度来说明菌群的差异变化，使用前面提及的`plot_bar`函数：

```R
plot_bar(GP, x="human", fill="SampleType", facet_grid= ~ Phylum)
```

![](https://static.xmt.cn/c1d9f448ff014d64acec625f1bfe1c55.png)

从上面这个图中也可以看出一些模式：

- Cyanobacteria 和 Actinobacteria 在人类样品较少；
- Firmicutes 在人类样本中较多；
- Acidobacteria 和 Verrucomicrobia 在粪便样本中较多。

根据之前的研究，这些结果并不是非常令人惊讶，但希望它可以应用于你所拥有的其他数据集中。

### DPCoA 双主坐标分析

下面的示例展示了使用了`phyloseq`中的`plot_ordination()`函数来进行双主坐标分析 ：

```R
# Perform ordination
GP.dpcoa <- ordinate(GP, "DPCoA") 
# Generate default ordination bi-plot
pdpcoa <- 
  plot_ordination(
    physeq = GP, 
    ordination = GP.dpcoa, 
    type="biplot",
    color="SampleType", 
    shape="Phylum")
# Adjust the shape scale manually 
# to make taxa hollow and samples filled (advanced)
shape.fac <- pdpcoa$data$Phylum
man.shapes <- c(19, 21:25)
names(man.shapes) <- c("Samples", levels(shape.fac)[levels(shape.fac)!="Samples"])
p2dpcoa <- pdpcoa + scale_shape_manual(values=man.shapes)
p2dpcoa
```

![](https://static.xmt.cn/c2d73180712344398c028b53b7139ba4.png)

```R
# Show just Samples or just Taxa
plot_ordination(GP, GP.dpcoa, type="taxa", shape="Phylum")
```

![](https://static.xmt.cn/b686e1696f5749dc8a3c5b963b584962.png)

```R
plot_ordination(GP, GP.dpcoa, type="samples", color="SampleType")
```

![](https://static.xmt.cn/7dfab912689b41a69b5522123d252334.png)

```R
# Split
plot_ordination(GP, GP.dpcoa, type="split",
                color="SampleType", shape="Phylum") +
  ggplot2::scale_colour_discrete()
```

![](https://static.xmt.cn/50f3ad80a4fd4cc29dd6d903f0d8bd87.png)

## 聚类分析

层级聚类（例如，`hclust`）是可视化样本距离矩阵的另一种流行的方式。在下面的示例中，我们使用了 unweighted UniFrac  距离矩阵和 UPGMA 方法（`hclust`参数： `method="average"`）：

```R
# (Re)load UniFrac distance matrix and GlobalPatterns data
data(GlobalPatterns)
load(system.file("doc", "Unweighted_UniFrac.RData", package="phyloseq"))
# Manually define color-shading vector based on sample type.
colorScale    <- rainbow(length(levels(get_variable(GlobalPatterns, "SampleType"))))
cols          <- colorScale[get_variable(GlobalPatterns, "SampleType")] 
GP.tip.labels <- as(get_variable(GlobalPatterns, "SampleType"), "character")
# This is the actual hierarchical clustering call, specifying average-link clustering
GP.hclust     <- hclust(GPUF, method="average")
plot(GP.hclust, col=cols)
```

![](https://static.xmt.cn/be206676fbae484bb61621855ca81445.png)

> 本文翻译整理自：http://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html
