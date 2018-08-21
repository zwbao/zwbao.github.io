'''
根据pmid从pubmed批量获得参考文献 Author-Date 格式
@author: zwbao
'''

from Bio import Entrez # pip install biopython
Entrez.email = "shinningbzw@foxmail.com" # email
output_file = r'ljj2citation.txt' # 输出文件
input_file = r'ljj.txt'# 输入文件：把excel中citation那一列保存为此文件

def citation(id):
    handle = Entrez.esummary(db="pubmed", id = id)
    record = Entrez.read(handle)

    title = record[0]['Title']
    pub_year = record[0]['PubDate'].split()[0]
    source = record[0]['Source']
    so = record[0]['SO'].split(";")[1].split(":")[0]
    pages = record[0]['Pages']
    cit = ""

    if len(record[0]['AuthorList']) > 3:
        cit = record[0]['AuthorList'][0] + ", " + record[0]['AuthorList'][1] + ", "+ record[0]['AuthorList'][2] + ", et al. "+"("+ pub_year + "). " + "\"" + title + "\" " + source + " " + so + ": " + pages + "."
    elif len(record[0]['AuthorList']) == 3:
        cit = record[0]['AuthorList'][0] + ", " + record[0]['AuthorList'][1] + " and "+ record[0]['AuthorList'][2] +" " +"("+ pub_year + "). " + "\"" + title + "\" " + source + " " + so + ": " + pages + "."
    elif len(record[0]['AuthorList']) == 2:
        cit = record[0]['AuthorList'][0] + " and " + record[0]['AuthorList'][1] +" " + "("+ pub_year + "). " + "\"" + title + "\" " + source + " " + so + ": " + pages + "."
    elif len(record[0]['AuthorList']) == 1:
        cit = record[0]['AuthorList'][0] + " " + "("+ pub_year + "). " + "\"" + title + "\" " + source + " " + so + ": " + pages + "."
    return cit

count = len(open(input_file, 'r').readlines())
line_c = []
print("Waiting...")

with open(output_file, 'a+', encoding='utf-8') as f:
    for line in open(input_file):
        new_cit_list = []
        if line.strip() != "citation_id":
            cit_list = line.strip().split(";")
            cit_sum_list = []

            if cit_list[0] != "":
                for i in cit_list:
                    if not str(i).startswith("NBK"):
                        new_cit_list.append(i)
                if len(new_cit_list) <= 3:
                    for i in range(len(new_cit_list)):
                        cit_sum_list.append(new_cit_list[i])
                        cit_sum_list.append(citation(new_cit_list[i]))
                else:
                    for i in range(3):
                        cit_sum_list.append(new_cit_list[i])
                        cit_sum_list.append(citation(new_cit_list[i])) 
                f.write("\t".join(cit_sum_list) + "\n")
            else:
                f.write("\n")
        line_c.append("b")
        if count % len(line_c) == 0:
            perc = (len(line_c) / count) * 100
            print("Completed " + str(int(perc)) + "%")