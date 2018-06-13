# 生信编程直播第12题-json格式数据的格式化

## 题目


json数据大家统一用我给的测试数据，自己下载： http://biotrainee.com/jbrowse/JBrowse-1.12.1/sample_data/json/modencode/modencodeMetaData.json 范例如下：

```
{
   "types" : {
      "data set" : {
         "pluralLabel" : "data sets"
      }
   },
   "items" : [
      {
         "technique" : "ChIP-chip",
         "factor" : "BEAF-32",
         "target" : "Non TF Chromatin binding factor",
         "principal_investigator" : "White, K.",
         "Tracks" : [
            "fly/White_INSULATORS_WIG/BEAF32"
         ],
         "submission" : "21",
         "label" : "BEAF-32;Embryos 0-12 hr;ChIP-chip",
         "category" : "Other chromatin binding sites",
         "type" : "data set",
         "Developmental-Stage" : "Embryos 0-12 hr",
         "organism" : "D. melanogaster"
      },
      {
         "technique" : "ChIP-chip",
         "factor" : "CP190",
         "target" : "Non TF Chromatin binding factor",
         "principal_investigator" : "White, K.",
         "Tracks" : [
            "fly/White_INSULATORS_WIG/CP190"
         ],
         "submission" : "22",
         "label" : "CP190;Embryos 0-12 hr;ChIP-chip",
         "category" : "Other chromatin binding sites",
         "type" : "data set",
         "Developmental-Stage" : "Embryos 0-12 hr",
         "organism" : "D. melanogaster"
      },
```

因为帖子长度有限，我就只截取了一部分，请自己下载查看，如果是完整的json，可以用在线工具查看结构： http://json.parser.online.fr/

如果不懂json格式的，请自行搜索哈，现在TCGA在GDC的metadata信息，就是json格式的。

我们需要从这个json文件里面提取：technique factor target principal_investigator submission label category type Developmental-Stage organism key  这几列信息，当然，是可以用正则表达式做的。
完成之后应该是：http://biotrainee.com/jbrowse/JBrowse-1.12.1/sample_data/json/modencode/modencodeMetaData.csv

```
"technique","factor","target","principal_investigator","submission","label","category","type","Developmental-Stage","organism","key"
"ChIP-chip","BEAF-32","Non TF Chromatin binding factor <b>with some bold text for testing</b>","White, K.","21","fly/White_INSULATORS_WIG/BEAF32","Other chromatin binding sites","data set","Embryos 0-12 hr","D. melanogaster","BEAF-32;Embryos 0-12 hr;ChIP-chip"
"ChIP-chip","CP190","Non TF Chromatin binding factor","White, K.","22","fly/White_INSULATORS_WIG/CP190","Other chromatin binding sites","data set","Embryos 0-12 hr","D. melanogaster","CP190;Embryos 0-12 hr;ChIP-chip"
"ChIP-chip","GAF","Non TF Chromatin binding factor","White, K.","23","fly/White_INSULATORS_WIG/GAF","Other chromatin binding sites","data set","Embryos 0-12 hr","D. melanogaster","GAF;Embryos 0-12 hr;ChIP-chip"
"ChIP-chip","mod(mdg4)","Non TF Chromatin binding factor","White, K.","24","fly/White_INSULATORS_WIG/MDG4","Other chromatin binding sites","data set","Embryos 0-12 hr","D. melanogaster","mod(mdg4);Embryos 0-12 hr;ChIP-chip"
"ChIP-chip","Su(Hw)","Non TF Chromatin binding factor","White, K.","27","fly/White_INSULATORS_WIG/SuHw","Other chromatin binding sites","data set","Embryos 0-12 hr","D. melanogaster","Su(Hw);Embryos 0-12 hr;ChIP-chip"
"RNA-tiling-array","total-RNA","mRNA","Celniker, S.","40","fly/CEL_CLINES/1182-4H","RNA expression profiling","data set","Late Embryonic stage","D. melanogaster","total-RNA;1182-4H;Late Embryonic stage;embryo-derived cell-line;RNA-tiling-array"
"ChIP-chip","HTZ-1","Chromatin Structure","Lieb, J.","43","worm/LIEB_WIG_CHIPCHIP_HVARS","Chromatin structure","data set","Mixed Embryos","C. elegans","HTZ-1;Mixed Embryos;20 degree celsius;ChIP-chip"
"ChIP-chip","pol2","Non TF Chromatin binding factor","Lieb, J.","44","worm/LIEB_WIG_CHIPCHIP_POL2/8WG16_N2_MXEMB","Other chromatin binding sites","data set","Mixed Embryos","C. elegans","pol2;Mixed Embryos;20 degree celsius;ChIP-chip"
"RNA-tiling-array","total-RNA","mRNA","Celniker, S.","48","fly/CEL_CLINES/CME-L1","RNA expression profiling","data set","Larvae 3rd instar","D. melanogaster","total-RNA;CME-L1;Larvae 3rd instar;ventral prothoracic disc;RNA-tiling-array"
```

我就不多做介绍了，主要难点在于理解json，本次作业，推荐大家用已有的包，正则表达式虽然可以做，但是太麻烦了~

给一个perl代码如下；

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';
use 5.10.0;
 
use JSON 2;
 
my $data = from_json( do { local $/; open my $f, '<', $ARGV[0]; scalar <$f> } );
 
my @fields = qw( technique factor target principal_investigator submission label category type Developmental-Stage organism key );
 
say join ',', map "\"$_\"", @fields;
 
for my $item ( @{$data->{items}} ) {
    $item->{key} = $item->{label};
    no warnings 'uninitialized';
    for my $track ( @{$item->{Tracks}} ) {
        $item->{label} = $track;
        say join ',', map "\"$_\"", @{$item}{@fields};
    }
}
```

希望有同学可以推陈出新，不要局限于我们的作业，可以自己下载TCGA的metadata信息，自己尝试提取，格式化。

## 我的答案

1. 我发现用简单的`sed`即可完成，先预处理一下：把前七行没用的信息去掉，再去掉换行符；把多出来的空格去掉，因为每一个item之间由`},`分割，故将其替换为换行符。

```
sed '1,7d' modencodeMetaData.json | tr '\n' ' ' | sed 's/{//g;s/\s//g;s/},/\n/g' >tmp.txt
```

```
head -2 tmp.txt
```

```
"technique":"ChIP-chip","factor":"BEAF-32","target":"NonTFChromatinbindingfactor","principal_investigator":"White,K.","Tracks":["fly/White_INSULATORS_WIG/BEAF32"],"submission":"21","label":"BEAF-32;Embryos0-12hr;ChIP-chip","category":"Otherchromatinbindingsites","type":"dataset","Developmental-Stage":"Embryos0-12hr","organism":"D.melanogaster"
"technique":"ChIP-chip","factor":"CP190","target":"NonTFChromatinbindingfactor","principal_investigator":"White,K.","Tracks":["fly/White_INSULATORS_WIG/CP190"],"submission":"22","label":"CP190;Embryos0-12hr;ChIP-chip","category":"Otherchromatinbindingsites","type":"dataset","Developmental-Stage":"Embryos0-12hr","organism":"D.melanogaster"
```

2. 将标签去掉，写个`for`循环即可

```
#!/bin/bash
#
for i in technique factor target principal_investigator submission label category type Developmental-Stage organism Tracks ;do
sed -i "s/\"$i\"\://g" tmp.txt
done
```

```
head -2 tmp.txt
```

```
"ChIP-chip","BEAF-32","NonTFChromatinbindingfactor","White,K.",["fly/White_INSULATORS_WIG/BEAF32"],"21","BEAF-32;Embryos0-12hr;ChIP-chip","Otherchromatinbindingsites","dataset","Embryos0-12hr","D.melanogaster"
"ChIP-chip","CP190","NonTFChromatinbindingfactor","White,K.",["fly/White_INSULATORS_WIG/CP190"],"22","CP190;Embryos0-12hr;ChIP-chip","Otherchromatinbindingsites","dataset","Embryos0-12hr","D.melanogaster"

```

3. 关于为什么`sed`不能替换换行符：是因为sed是按行处理文本数据的，每次处理一行数据后，都会在行尾自动添加trailing newline，其实就是行的分隔符即换行符。

- 详见：[linux sed命令，如何替换换行符“\n”](http://blog.csdn.net/u011729865/article/details/71773840)

---

:yum: 我还是一个小小白，平时的笔记一定有不少疏漏之处，恳请各位大佬批评指正！

