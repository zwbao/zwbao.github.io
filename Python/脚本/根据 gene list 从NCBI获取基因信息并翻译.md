# 根据 gene list 从NCBI获取基因信息并翻译

## Code

[get_gene_summary](/Python/脚本/get_gene_summary.py ':include :type=code')

```python
from Bio import Entrez # pip install biopython
from translate_api.translate_api import api # pip install translate_api
import re

rm_pattern = re.compile('\[.*?\]')

for line in open('gene.txt'): # 输入去重的基因列表
    gene = str(line.strip())
    gene_term = "(" + gene +"[Gene Name]) AND Homo sapiens[Organism]"
    Entrez.email = "shinningbzw@foxmail.com"
    handle = Entrez.esearch(db="gene", term=gene_term)
    gene_id = Entrez.read(handle)['IdList'][0]
    sum_handle = Entrez.esummary(db="gene", id=gene_id)
    sum_record = Entrez.read(sum_handle)
    r_gene_sum = sum_record['DocumentSummarySet']['DocumentSummary'][0]['Summary']
    gene_sum = rm_pattern.sub('', r_gene_sum)
    translation = api(gene_sum)
    print(gene + "\n" + gene_sum + "\n" + translation + "\n")
```

## Result

```
ZEB2
The protein encoded by this gene is a member of the Zfh1 family of 2-handed zinc finger/homeodomain proteins. It is located in the nucleus and functions as a DNA-binding transcriptional repressor that interacts with activated SMADs. Mutations in this gene are associated with Hirschsprung disease/Mowat-Wilson syndrome. Alternatively spliced transcript variants have been found for this gene.
由该基因编码的蛋白质是Zfh1家族的双手锌指/同源域蛋白质的成员。它位于细胞核中，起着与受激活的SMAD相互作用的DNA结合转录抑制因子的作用。该基因的突变与先天性巨结肠病/ Mowat-Wilson综合征有关。已经发现该基因的剪接转录物变体。

DMD
This gene spans a genomic range of greater than 2 Mb and encodes a large protein containing an N-terminal actin-binding domain and multiple spectrin repeats. The encoded protein forms a component of the dystrophin-glycoprotein complex (DGC), which bridges the inner cytoskeleton and the extracellular matrix. Deletions, duplications, and point mutations at this gene locus may cause Duchenne muscular dystrophy (DMD), Becker muscular dystrophy (BMD), or cardiomyopathy. Alternative promoter usage and alternative splicing result in numerous distinct transcript variants and protein isoforms for this gene. 
该基因跨越大于2Mb的基因组范围，并编码含有N-末端肌动蛋白结合结构域和多个血影蛋白重复序列​​的大蛋白质。编码的蛋白质形成肌营养不良蛋白 - 糖蛋白复合物（DGC）的组分，其结合内部细胞骨架和细胞外基质。该基因位点的缺失，重复和点突变可能导致杜氏肌营养不良症（DMD），贝克尔肌营养不良症（BMD）或心肌病。替代的启动子用法和可变剪接导致该基因的许多不同的转录物变体和蛋白质同种型。

GJB3
This gene is a member of the connexin gene family. The encoded protein is a component of gap junctions, which are composed of arrays of intercellular channels that provide a route for the diffusion of low molecular weight materials from cell to cell. Mutations in this gene can cause non-syndromic deafness or erythrokeratodermia variabilis, a skin disorder. Alternative splicing results in multiple transcript variants encoding the same protein. 
该基因是连接蛋白基因家族的成员。编码的蛋白质是间隙连接的组分，其由细胞间通道阵列组成，其提供了从细胞到细胞的低分子量材料扩散的途径。该基因的突变可引起非综合征性耳聋或变性红皮肤角化症，即皮肤病。选择性剪接导致编码相同蛋白质的多种转录物变体。
```