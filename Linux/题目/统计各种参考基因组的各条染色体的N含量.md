# 统计各种参考基因组的各条染色体的N含量

## 题目

染色体里面的碱基一般只是ATCG这样正常的碱基，但是参考基因组毕竟不完美，比如人类的hg19参考基因组里面的chrY就高达36.38M的N碱基，而它全长还不到60M。如果计算覆盖度的时候没有能考虑到N这样的碱基，就会到底一些非常大的误差。因此统计各种参考基因组的各条染色体的N含量。

## 我的答案

```
awk 'BEGIN{RS=">";FS="\n"} NR>1{seq="";for(i=2;i<=NF;i++) seq=seq""$i;len =length(seq);count_N=0;for(i=0;i<=len;i++){tmp = substr(seq,i+1,1);if(tmp ~ /[N]/) count_N++};printf $1"\t"count_N}' test.fa
```

