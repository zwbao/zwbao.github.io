# 生信编程直播第五题-根据GTF画基因的多个转录本结构

```
[root@VM_180_248_centos d6t]# grep ^C hsa00001.keg | grep 'hsa' |wc -l
364
```

```
awk '{if($1 == C)path=$2}'


#!/bin/bash
#


awk '{if($1 == C){path=$2;printf $path"\t";if($1 ==D )gene = $2;printf $gene"\t"}}' hsa00001.keg |head

```