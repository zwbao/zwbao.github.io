'''
根据基因列表从NCBI获得基因功能并通过谷歌翻译获得中文释义 
@author: zwbao
'''


from Bio import Entrez # pip install biopython
from translate_api.translate_api import api # pip install translate_api
import re

Entrez.email = "shinningbzw@foxmail.com" # email
output_file = r'gene_sum.txt' # 输出文件
input_file = r'genelist.txt'# 输入文件：去重后的基因列表 （将基因列保存为 txt，uniq *.txt>gene_list.txt ）

gene_list = []
line_c = []
count = len(open(input_file, 'r').readlines())
print("Waiting...")

# get gene list
for line in open(input_file):
    if line != "基因":
        gene_list.append(line)

gene_list.remove(gene_list[0])
rm_pattern = re.compile('\[.*?\]')

with open(output_file, 'a+', encoding='utf-8') as f:
    for line in gene_list:
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
        f.write(gene + "\n" + gene_sum + "\n" + translation + "\n")
        line_c.append("b")
        if count % len(line_c) == 0:
            perc = (len(line_c) / count) * 100
            print("Completed " + str(int(perc)) + "%")