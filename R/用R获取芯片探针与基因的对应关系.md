# Bioconductor注释专题：用R获取芯片探针与基因的对应关系

今天要讲的是如何用R的bioconductor包来得到芯片探针与基因的对应关系~

一般重要的芯片在R的bioconductor里面都是有包的，不同的芯片对应不同的包，常见的如下：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180407/jmdKEfa801.png?imageslim)

先安装AnnotationDbi

```R
source("http://bioconductor.org/biocLite.R")
biocLite("AnnotationDbi")
```

以`hgu95av2.db`为例，下载对应的数据库：

```R
biocLite("hgu95av2.db")
```

然后载入这两个包

```R
library(AnnotationDbi)
library(hgu95av2.db)
```

看下数据库的信息~

```R
> hgu95av2.db
ChipDb object:
| DBSCHEMAVERSION: 2.1
| Db type: ChipDb
| Supporting package: AnnotationDbi
| DBSCHEMA: HUMANCHIP_DB
| ORGANISM: Homo sapiens
| SPECIES: Human
| MANUFACTURER: Affymetrix
| CHIPNAME: Human Genome U95 Set
| MANUFACTURERURL: http://www.affymetrix.com/support/technical/byproduct.affx?product=hgu95
| EGSOURCEDATE: 2015-Sep27
| EGSOURCENAME: Entrez Gene
| EGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
| CENTRALID: ENTREZID
| TAXID: 9606
| GOSOURCENAME: Gene Ontology
| GOSOURCEURL: ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/
| GOSOURCEDATE: 20150919
| GOEGSOURCEDATE: 2015-Sep27
| GOEGSOURCENAME: Entrez Gene
| GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
| KEGGSOURCENAME: KEGG GENOME
| KEGGSOURCEURL: ftp://ftp.genome.jp/pub/kegg/genomes
| KEGGSOURCEDATE: 2011-Mar15
| GPSOURCENAME: UCSC Genome Bioinformatics (Homo sapiens)
| GPSOURCEURL: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19
| GPSOURCEDATE: 2010-Mar22
| ENSOURCEDATE: 2015-Jul16
| ENSOURCENAME: Ensembl
| ENSOURCEURL: ftp://ftp.ensembl.org/pub/current_fasta
| UPSOURCENAME: Uniprot
| UPSOURCEURL: http://www.uniprot.org/
| UPSOURCEDATE: Thu Oct  1 23:31:58 2015

Please see: help('select') for usage information
```

库所包含的内容及可以作为检索键的列分别可以用`columns`命令和`keytypes`命令查看： 

```R
> columns(hgu95av2.db)
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"     
 [4] "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
 [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
[10] "GENENAME"     "GO"           "GOALL"       
[13] "IPI"          "MAP"          "OMIM"        
[16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"        
[19] "PFAM"         "PMID"         "PROBEID"     
[22] "PROSITE"      "REFSEQ"       "SYMBOL"      
[25] "UCSCKG"       "UNIGENE"      "UNIPROT"     
> keytypes(hgu95av2.db)
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"     
 [4] "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
 [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
[10] "GENENAME"     "GO"           "GOALL"       
[13] "IPI"          "MAP"          "OMIM"        
[16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"        
[19] "PFAM"         "PMID"         "PROBEID"     
[22] "PROSITE"      "REFSEQ"       "SYMBOL"      
[25] "UCSCKG"       "UNIGENE"      "UNIPROT"   
```

要是想要看上面的具体某一列可以用`key`

```R
> head(keys(hgu95av2.db, keytype="SYMBOL"))
[1] "A1BG"  "A2M"   "A2MP1" "NAT1"  "NAT2"  "NATP" 
```

最后，要是我们有一些PROBEID需要转换成SYMBOL，可以这么做：

```R
> # 模拟一些PROBEID
> k <- head(keys(hgu95av2.db,keytype="PROBEID"))
> # 使用select进行选择
> select(hgu95av2.db, keys=k, columns=c("SYMBOL"), keytype="PROBEID")
'select()' returned 1:1 mapping between keys and
columns
    PROBEID  SYMBOL
1   1000_at   MAPK3
2   1001_at    TIE1
3 1002_f_at CYP2C19
4 1003_s_at   CXCR5
5   1004_at   CXCR5
6   1005_at   DUSP1
```

`select`需要先指定你要使用的数据库，这里就是hgu95av2.db，接下来的keys参数是要检索的key，可以是所有id，或是感兴趣的id的列表。columns给出的是你要检索的对应内容，我们这里是根据id来检索基因的symbol，因此在columns参数中只指定了这一项。

- 参考：
    - http://bioconductor.org/packages/2.12/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
    - http://blog.chinaunix.net/uid-12084847-id-3851353.html
    - http://www.bio-info-trainee.com/1399.html