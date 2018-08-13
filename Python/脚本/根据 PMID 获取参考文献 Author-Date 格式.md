# 根据 PMID 获取参考文献 Author-Date 格式


## Code

```python
from Bio import Entrez # pip install biopython
Entrez.email = "shinningbzw@foxmail.com"

def citation(id):
    handle = Entrez.esummary(db="pubmed", id = id)
    record = Entrez.read(handle)

    title = record[0]['Title']
    pub_year = record[0]['PubDate'].split()[0]
    source = record[0]['Source']
    so = record[0]['SO'].split(";")[1].split(":")[0]
    pages = record[0]['Pages']

    if len(record[0]['AuthorList']) > 3: # 最多三个作者
        cit = record[0]['AuthorList'][0] + ", " + record[0]['AuthorList'][1] + ", "+ record[0]['AuthorList'][2] + ", et al. "+"("+ pub_year + "). " + "\"" + title + "\" " + source + " " + so + ": " + pages + "."
    else:
        cit = record[0]['AuthorList'][0] + ", " + record[0]['AuthorList'][1] + ", "+ record[0]['AuthorList'][2] +"("+ pub_year + "). " + "\"" + title + "\" " + source + " " + so + ": " + pages + "."
    
    return cit

for line in open('zwbao.txt'):
    pmid = line.split("\t")[18]
    if pmid != "citation_id":
        cit_list = pmid.strip().split(";")
        cit_sum_list = []

        if cit_list[0] != "":
            if len(cit_list) <= 3:
                for i in range(len(cit_list)):
                    cit_sum_list.append(cit_list[i])
                    cit_sum_list.append(citation(cit_list[i]))
            else: # 最多插入三篇参考文献
                for i in range(3):
                    cit_sum_list.append(cit_list[i])
                    cit_sum_list.append(citation(cit_list[i])) 
            print("\t".join(cit_sum_list))
        else:
            print("")
```

## Result

```

11279515	Wakamatsu N, Yamada Y, Yamada K, et al. (2001). "Mutations in SIP1, encoding Smad interacting protein-1, cause a form of Hirschsprung disease." Nat Genet 27(4): 369-70.
11279515	Wakamatsu N, Yamada Y, Yamada K, et al. (2001). "Mutations in SIP1, encoding Smad interacting protein-1, cause a form of Hirschsprung disease." Nat Genet 27(4): 369-70.	11592033	Yamada K, Yamada Y, Nomura N, et al. (2001). "Nonsense and frameshift mutations in ZFHX1B, encoding Smad-interacting protein 1, cause a complex developmental disorder with a great variety of clinical features." Am J Hum Genet 69(6): 1178-85.	18414213	Richards CS, Bale S, Bellissimo DB, et al. (2008). "ACMG recommendations for standards for interpretation and reporting of sequence variations: Revisions 2007." Genet Med 10(4): 294-300.
11891681	Zweier C, Albrecht B, Mitulla B, et al. (2002). ""Mowat-Wilson" syndrome with and without Hirschsprung disease is a distinct, recognizable multiple congenital anomalies-mental retardation syndrome caused by mutations in the zinc finger homeo box 1B gene." Am J Med Genet 108(3): 177-81.
16688751	Heinritz W, Zweier C, Froster UG, et al. (2006). "A missense mutation in the ZFHX1B gene associated with an atypical Mowat-Wilson syndrome phenotype." Am J Med Genet A 140(11): 1223-7.
```