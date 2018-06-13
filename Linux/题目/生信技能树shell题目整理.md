# 生信技能树shell题目整理
## 批量根据基因list来提取信息
根据多个基因list来提取fasta文件里面指定序列名的序列
```shell
ls *txt |while read id
do
cat $id ~/annotation/CHIPseq/mm10/ucsc.refseq.bed |perl -alne '{$h{$F[0]}=1;print if exists $h{$F[3]} }' >${id%%.*}.bed
done
```
## 批量从NCBI下载数据
首先进入https://www.ncbi.nlm.nih.gov/genome/genomes/13563 可以找到Mycobacterium相关的180个记录，现在需要批量下载每一个记录的4个数据，回到搜索界面，点击download table，下载文档，里面有所有菌株的列表信息，`cut`有ftp地址的那一列，拿到ftp地址如下：
```
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/015/405/GCA_000015405.1_ASM1540v1
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/620/625/GCA_000620625.1_ASM62062v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 620625.1_ASM62062v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/972/925/GCA_000972925.1_ASM97292v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 972925.1_ASM97292v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/021/385/GCA_001021385.1_ASM102138v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 21385.1_ASM102138v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/328/565/GCA_000328565.1_ASM32856v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 328565.1_ASM32856v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/426/545/GCA_001426545.1_Root135]ftp://ftp.ncbi.nlm.nih.gov/genom ... 001426545.1_Root135[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/428/895/GCA_001428895.1_Root265]ftp://ftp.ncbi.nlm.nih.gov/genom ... 001428895.1_Root265[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/494/595/GCA_001494595.1_Mycobacterium_massilipolynesiensis]ftp://ftp.ncbi.nlm.nih.gov/genom ... assilipolynesiensis[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/855/GCA_001499855.1_ASM149985v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99855.1_ASM149985v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/915/GCA_001499915.1_ASM149991v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99915.1_ASM149991v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/373/905/GCA_000373905.1_ASM37390v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 373905.1_ASM37390v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/014/165/GCA_000014165.1_ASM1416v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 0014165.1_ASM1416v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/016/005/GCA_000016005.1_ASM1600v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 0016005.1_ASM1600v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/262/165/GCA_000262165.1_ASM26216v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 262165.1_ASM26216v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/416/365/GCA_000416365.2_ASM41636v2]ftp://ftp.ncbi.nlm.nih.gov/genom ... 416365.2_ASM41636v2[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/580/405/GCA_001580405.1_ASM158040v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 80405.1_ASM158040v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/611/855/GCA_001611855.1_ASM161185v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 11855.1_ASM161185v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/644/575/GCA_001644575.1_ASM164457v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 44575.1_ASM164457v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/695/755/GCA_001695755.1_ASM169575v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 95755.1_ASM169575v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/886/515/GCA_001886515.1_ASM188651v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 86515.1_ASM188651v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/984/215/GCA_001984215.1_ASM198421v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 84215.1_ASM198421v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/364/405/GCA_000364405.1_ASM36440v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 364405.1_ASM36440v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/382/405/GCA_000382405.1_ASM38240v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 382405.1_ASM38240v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/383/495/GCA_000383495.1_ASM38349v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 383495.1_ASM38349v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/426/065/GCA_000426065.1_ASM42606v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 426065.1_ASM42606v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/686/745/GCA_000686745.1_ASM68674v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 686745.1_ASM68674v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/709/305/GCA_000709305.1_Myco_sp_TKK-01-0059_V2]ftp://ftp.ncbi.nlm.nih.gov/genom ... o_sp_TKK-01-0059_V2[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/746/215/GCA_000746215.1_ASM74621v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 746215.1_ASM74621v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/428/285/GCA_001428285.1_Soil538]ftp://ftp.ncbi.nlm.nih.gov/genom ... 001428285.1_Soil538[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/157/375/GCA_900157375.1_PRJEB19165]ftp://ftp.ncbi.nlm.nih.gov/genom ... 157375.1_PRJEB19165[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/157/385/GCA_900157385.1_PRJEB19184]ftp://ftp.ncbi.nlm.nih.gov/genom ... 157385.1_PRJEB19184[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/416/385/GCA_000416385.1_M6May]ftp://ftp.ncbi.nlm.nih.gov/genom ... A_000416385.1_M6May[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/419/295/GCA_000419295.1_SP]ftp://ftp.ncbi.nlm.nih.gov/genom ... /GCA_000419295.1_SP[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/455/125/GCA_000455125.1_Mycobacterium_sp_UM_WGJ]ftp://ftp.ncbi.nlm.nih.gov/genom ... bacterium_sp_UM_WGJ[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/455/185/GCA_000455185.1_Mycobacterium_sp_UM_RHS]ftp://ftp.ncbi.nlm.nih.gov/genom ... bacterium_sp_UM_RHS[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/455/205/GCA_000455205.1_Mycobacterium_sp_UM_CSW]ftp://ftp.ncbi.nlm.nih.gov/genom ... bacterium_sp_UM_CSW[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/523/635/GCA_000523635.1_ASM52363v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 523635.1_ASM52363v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/523/975/GCA_000523975.1_ASM52397v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 523975.1_ASM52397v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/972/905/GCA_000972905.1_ASM97290v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 972905.1_ASM97290v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/972/915/GCA_000972915.1_ASM97291v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 972915.1_ASM97291v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/972/935/GCA_000972935.1_ASM97293v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 972935.1_ASM97293v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/414/075/GCA_001414075.1_ASM141407v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 14075.1_ASM141407v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/414/095/GCA_001414095.1_ASM141409v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 14095.1_ASM141409v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/005/GCA_001440005.1_ASM144000v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40005.1_ASM144000v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/085/GCA_001440085.1_ASM144008v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40085.1_ASM144008v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/105/GCA_001440105.1_ASM144010v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40105.1_ASM144010v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/125/GCA_001440125.1_ASM144012v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40125.1_ASM144012v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/135/GCA_001440135.1_ASM144013v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40135.1_ASM144013v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/155/GCA_001440155.1_ASM144015v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40155.1_ASM144015v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/185/GCA_001440185.1_ASM144018v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40185.1_ASM144018v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/205/GCA_001440205.1_ASM144020v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40205.1_ASM144020v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/225/GCA_001440225.1_ASM144022v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40225.1_ASM144022v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/245/GCA_001440245.1_ASM144024v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40245.1_ASM144024v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/265/GCA_001440265.1_ASM144026v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40265.1_ASM144026v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/275/GCA_001440275.1_ASM144027v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40275.1_ASM144027v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/440/305/GCA_001440305.1_ASM144030v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 40305.1_ASM144030v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/825/GCA_001499825.1_ASM149982v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99825.1_ASM149982v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/835/GCA_001499835.1_ASM149983v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99835.1_ASM149983v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/845/GCA_001499845.1_ASM149984v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99845.1_ASM149984v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/905/GCA_001499905.1_ASM149990v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99905.1_ASM149990v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/925/GCA_001499925.1_ASM149992v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99925.1_ASM149992v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/965/GCA_001499965.1_ASM149996v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99965.1_ASM149996v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/985/GCA_001499985.1_ASM149998v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99985.1_ASM149998v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/499/995/GCA_001499995.1_ASM149999v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 99995.1_ASM149999v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/500/025/GCA_001500025.1_ASM150002v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 00025.1_ASM150002v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/500/045/GCA_001500045.1_ASM150004v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 00045.1_ASM150004v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/500/065/GCA_001500065.1_ASM150006v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 00065.1_ASM150006v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/500/085/GCA_001500085.1_ASM150008v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 00085.1_ASM150008v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/500/105/GCA_001500105.1_ASM150010v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 00105.1_ASM150010v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/500/125/GCA_001500125.1_ASM150012v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 00125.1_ASM150012v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/500/145/GCA_001500145.1_ASM150014v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 00145.1_ASM150014v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/545/925/GCA_001545925.1_ASM154592v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 45925.1_ASM154592v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/235/GCA_001665235.1_ASM166523v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65235.1_ASM166523v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/255/GCA_001665255.1_ASM166525v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65255.1_ASM166525v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/295/GCA_001665295.1_ASM166529v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65295.1_ASM166529v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/365/GCA_001665365.1_ASM166536v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65365.1_ASM166536v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/395/GCA_001665395.1_ASM166539v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65395.1_ASM166539v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/535/GCA_001665535.1_ASM166553v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65535.1_ASM166553v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/575/GCA_001665575.1_ASM166557v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65575.1_ASM166557v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/605/GCA_001665605.1_ASM166560v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65605.1_ASM166560v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/615/GCA_001665615.1_ASM166561v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65615.1_ASM166561v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/645/GCA_001665645.1_ASM166564v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65645.1_ASM166564v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/685/GCA_001665685.1_ASM166568v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65685.1_ASM166568v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/755/GCA_001665755.1_ASM166575v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65755.1_ASM166575v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/825/GCA_001665825.1_ASM166582v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65825.1_ASM166582v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/665/875/GCA_001665875.1_ASM166587v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 65875.1_ASM166587v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/745/GCA_001666745.1_ASM166674v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66745.1_ASM166674v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/755/GCA_001666755.1_ASM166675v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66755.1_ASM166675v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/785/GCA_001666785.1_ASM166678v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66785.1_ASM166678v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/815/GCA_001666815.1_ASM166681v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66815.1_ASM166681v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/835/GCA_001666835.1_ASM166683v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66835.1_ASM166683v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/865/GCA_001666865.1_ASM166686v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66865.1_ASM166686v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/875/GCA_001666875.1_ASM166687v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66875.1_ASM166687v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/895/GCA_001666895.1_ASM166689v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66895.1_ASM166689v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/915/GCA_001666915.1_ASM166691v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66915.1_ASM166691v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/935/GCA_001666935.1_ASM166693v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66935.1_ASM166693v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/666/955/GCA_001666955.1_ASM166695v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66955.1_ASM166695v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/015/GCA_001667015.1_ASM166701v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67015.1_ASM166701v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/035/GCA_001667035.1_ASM166703v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67035.1_ASM166703v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/065/GCA_001667065.1_ASM166706v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67065.1_ASM166706v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/075/GCA_001667075.1_ASM166707v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67075.1_ASM166707v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/105/GCA_001667105.1_ASM166710v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67105.1_ASM166710v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/115/GCA_001667115.1_ASM166711v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67115.1_ASM166711v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/145/GCA_001667145.1_ASM166714v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67145.1_ASM166714v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/155/GCA_001667155.1_ASM166715v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67155.1_ASM166715v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/185/GCA_001667185.1_ASM166718v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67185.1_ASM166718v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/265/GCA_001667265.1_ASM166726v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67265.1_ASM166726v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/275/GCA_001667275.1_ASM166727v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67275.1_ASM166727v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/315/GCA_001667315.1_ASM166731v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67315.1_ASM166731v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/425/GCA_001667425.1_ASM166742v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67425.1_ASM166742v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/455/GCA_001667455.1_ASM166745v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67455.1_ASM166745v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/465/GCA_001667465.1_ASM166746v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67465.1_ASM166746v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/505/GCA_001667505.1_ASM166750v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67505.1_ASM166750v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/535/GCA_001667535.1_ASM166753v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67535.1_ASM166753v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/585/GCA_001667585.1_ASM166758v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67585.1_ASM166758v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/595/GCA_001667595.1_ASM166759v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67595.1_ASM166759v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/625/GCA_001667625.1_ASM166762v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67625.1_ASM166762v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/665/GCA_001667665.1_ASM166766v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67665.1_ASM166766v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/695/GCA_001667695.1_ASM166769v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67695.1_ASM166769v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/735/GCA_001667735.1_ASM166773v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67735.1_ASM166773v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/745/GCA_001667745.1_ASM166774v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67745.1_ASM166774v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/775/GCA_001667775.1_ASM166777v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67775.1_ASM166777v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/785/GCA_001667785.1_ASM166778v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67785.1_ASM166778v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/835/GCA_001667835.1_ASM166783v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67835.1_ASM166783v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/865/GCA_001667865.1_ASM166786v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67865.1_ASM166786v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/925/GCA_001667925.1_ASM166792v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67925.1_ASM166792v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/667/995/GCA_001667995.1_ASM166799v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 67995.1_ASM166799v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/668/575/GCA_001668575.1_ASM166857v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 68575.1_ASM166857v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/668/615/GCA_001668615.1_ASM166861v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 68615.1_ASM166861v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/668/625/GCA_001668625.1_ASM166862v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 68625.1_ASM166862v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/668/695/GCA_001668695.1_ASM166869v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 68695.1_ASM166869v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/668/725/GCA_001668725.1_ASM166872v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 68725.1_ASM166872v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/669/335/GCA_001669335.1_ASM166933v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 69335.1_ASM166933v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/665/GCA_001672665.1_ASM167266v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72665.1_ASM167266v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/675/GCA_001672675.1_ASM167267v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72675.1_ASM167267v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/685/GCA_001672685.1_ASM167268v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72685.1_ASM167268v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/745/GCA_001672745.1_ASM167274v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72745.1_ASM167274v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/815/GCA_001672815.1_ASM167281v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72815.1_ASM167281v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/895/GCA_001672895.1_ASM167289v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72895.1_ASM167289v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/915/GCA_001672915.1_ASM167291v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72915.1_ASM167291v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/935/GCA_001672935.1_ASM167293v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72935.1_ASM167293v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/975/GCA_001672975.1_ASM167297v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72975.1_ASM167297v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/672/995/GCA_001672995.1_ASM167299v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 72995.1_ASM167299v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/673/055/GCA_001673055.1_ASM167305v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 73055.1_ASM167305v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/673/155/GCA_001673155.1_ASM167315v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 73155.1_ASM167315v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/673/235/GCA_001673235.1_ASM167323v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 73235.1_ASM167323v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/673/405/GCA_001673405.1_ASM167340v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 73405.1_ASM167340v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/673/415/GCA_001673415.1_ASM167341v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 73415.1_ASM167341v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/673/535/GCA_001673535.1_ASM167353v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 73535.1_ASM167353v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/673/555/GCA_001673555.1_ASM167355v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 73555.1_ASM167355v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/673/615/GCA_001673615.1_ASM167361v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 73615.1_ASM167361v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/766/635/GCA_001766635.1_ASM176663v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 66635.1_ASM176663v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/853/525/GCA_001853525.1_ASM185352v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 53525.1_ASM185352v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/905/305/GCA_001905305.1_ASM190530v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 05305.1_ASM190530v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/905/565/GCA_001905565.1_ASM190556v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 05565.1_ASM190556v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/905/655/GCA_001905655.1_ASM190565v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 05655.1_ASM190565v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/907/615/GCA_001907615.1_ASM190761v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 07615.1_ASM190761v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/942/625/GCA_001942625.1_ASM194262v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 42625.1_ASM194262v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/953/975/GCA_001953975.1_ASM195397v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 53975.1_ASM195397v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/954/045/GCA_001954045.1_ASM195404v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 54045.1_ASM195404v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/954/135/GCA_001954135.1_ASM195413v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 54135.1_ASM195413v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/954/195/GCA_001954195.1_ASM195419v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 54195.1_ASM195419v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/954/215/GCA_001954215.1_ASM195421v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 54215.1_ASM195421v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/954/275/GCA_001954275.1_ASM195427v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 54275.1_ASM195427v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/157/365/GCA_900157365.1_PRJEB19151]ftp://ftp.ncbi.nlm.nih.gov/genom ... 157365.1_PRJEB19151[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/455/145/GCA_000455145.1_Mycobacterium_sp_UM_WWY]ftp://ftp.ncbi.nlm.nih.gov/genom ... bacterium_sp_UM_WWY[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/744/355/GCA_000744355.1_ASM74435v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 744355.1_ASM74435v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/987/455/GCA_000987455.1_ASM98745v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 987455.1_ASM98745v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/756/795/GCA_001756795.1_ASM175679v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 56795.1_ASM175679v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/835/505/GCA_001835505.1_ASM183550v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 35505.1_ASM183550v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/835/515/GCA_001835515.1_ASM183551v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 35515.1_ASM183551v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/835/525/GCA_001835525.1_ASM183552v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 35525.1_ASM183552v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/835/535/GCA_001835535.1_ASM183553v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 35535.1_ASM183553v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/854/525/GCA_001854525.1_ASM185452v1]ftp://ftp.ncbi.nlm.nih.gov/genom ... 54525.1_ASM185452v1[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/050/035/GCA_001050035.1_Mycobacterium_komanii_GPK_1020]ftp://ftp.ncbi.nlm.nih.gov/genom ... um_komanii_GPK_1020[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/100/615/GCA_900100615.1_IMG-taxon_2593339259_annotated_assembly]ftp://ftp.ncbi.nlm.nih.gov/genom ... _annotated_assembly[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/101/555/GCA_900101555.1_IMG-taxon_2639762559_annotated_assembly]ftp://ftp.ncbi.nlm.nih.gov/genom ... _annotated_assembly[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/107/555/GCA_900107555.1_IMG-taxon_2636416058_annotated_assembly]ftp://ftp.ncbi.nlm.nih.gov/genom ... _annotated_assembly[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/110/825/GCA_900110825.1_IMG-taxon_2642422555_annotated_assembly]ftp://ftp.ncbi.nlm.nih.gov/genom ... _annotated_assembly[/url]
[url=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/113/075/GCA_900113075.1_IMG-taxon_2639762525_annotated_assembly]ftp://ftp.ncbi.nlm.nih.gov/genom ... _annotated_assembly[/url]
```
然后写一个shell脚本即可：
```shell
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.fna.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.gff.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.gbff.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.faa.gz $id;done
```
- 解析：
1. 首先下载得到 genomes_proks.txt，第一行内容为：
```
#Organism/Name	Strain	CladeID	BioSample	BioProject	Group	SubGroup	Assembly	Size (Mb)	GC%	Replicons	WGS	Scaffolds	Genes	Proteins	Release Date	Modify Date	Level	RefSeq FTP	GenBank FTP
```
2. 使用以下代码切割：
```shell
cut -f19,20    #默认的字段分隔符为“TAB”；-f：显示指定字段的内容
```
得到ftp地址
3.用shell脚本批量下载：
```shell
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.fna.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.gff.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.gbff.gz $id;done
cat ftp.list |while read id; do wget -c -r -np -k -L -p  -nd -A.faa.gz $id;done
```
使用到的命令：
- read：read命令从键盘读取变量的值，通常用在shell脚本中与用户进行交互的场合。
- wget：wget命令用来从指定的URL下载文件。
```
-c   断点续传
-r   递归下载，下载指定网页某一目录下（包括子目录）的所有文件
-np 递归下载时不搜索上层目录，如wget -c -r www.xxx.org/pub/path/，没有加参数-np，就会同时下载path的上一级目录pub下的其它文件
-k 将绝对链接转为相对链接，下载整个站点后脱机浏览网页，最好加上这个参数
-L 递归时不进入其它主机
-p 下载网页所需的所有文件
-nd 递归下载时不创建一层一层的目录，把所有的文件下载到当前目录
-A<后缀名>：指定要下载文件的后缀名，多个后缀名之间使用逗号进行分隔
```
## 区分染色体分别运行scalpel软件
- Scalpel is available here: http://scalpel.sourceforge.net/ 
- 文章是： http://www.nature.com/nmeth/journal/v11/n10/full/nmeth.3069.html 
- 软件说明书写的也比较详细：http://scalpel.sourceforge.net/manual.html

他提供了3种情况的找INDELs变异，我目前需要用的就是对我的全基因组测序数据来找，所以用single模式；为了节省对计算资源的消耗，作者建议我单独对每条染色体分别处理。

1. 软件安装
```shell
## Download and install Scalpel
cd ~/biosoft
mkdir Scalpel &&  cd Scalpel
wget [url=https://downloads.sourceforge.net/project/scalpel/scalpel-0.5.3.tar.gz]https://downloads.sourceforge.ne ... calpel-0.5.3.tar.gz[/url]  
tar zxvf scalpel-0.5.3.tar.gz
cd scalpel-0.5.3
make
~/biosoft/Scalpel/scalpel-0.5.3/scalpel-discovery  --help
~/biosoft/Scalpel/scalpel-0.5.3/scalpel-export  --help
```
它需要自己指定--bed参数来选择染色体运行，而且不是给一个chr1就可以了，需要指定染色体及其起始终止坐标：single region in format chr:start-end (example: 1:31656613-31656883)
所以就考验shell编程技巧啦！

2. 制作 ~/reference/genome/hg19/hg19.chr.bed  这个文件，我就不多说了，前面我们已经讲过了！
```
chr10        1        135534747
chr11        1        135006516
chr12        1        133851895
chr13        1        115169878
chr14        1        107349540
chr15        1        102531392
chr16        1        90354753
chr17        1        81195210
chr18        1        78077248
chr19        1        59128983
chr1        1        249250621
chr20        1        63025520
chr21        1        48129895
chr22        1        51304566
chr2        1        243199373
chr3        1        198022430
chr4        1        191154276
chr5        1        180915260
chr6        1        171115067
chr7        1        159138663
chr8        1        146364022
chr9        1        141213431
```
3. 区分染色体分别运行scalpel软件代码如下：
```shell
cat ~/reference/genome/hg19/hg19.chr.bed |while read id
do
arr=($id) 
 
# arr=($a) will split the $a to $arr , ${arr[0]} ${arr[1]} ~~~, but ${arr[@]}  is the whole array .
# OLD_IFS="$IFS" 
# IFS="," 
# arr=($a) 
# IFS="$OLD_IFS" 
 
#arr=($a)用于将字符串$a分割到数组$arr ${arr[0]} ${arr[1]} ... 分别存储分割后的数组第1 2 ... 项 ，${arr[@]}存储整个数组。
#变量$IFS存储着分隔符，这里我们将其设为逗号 "," OLD_IFS用于备份默认的分隔符，使用完后将之恢复默认。
 
echo ${arr[0]}:${arr[1]}-${arr[2]}
  
 
date
start=`date +%s`
 
~/biosoft/Scalpel/scalpel-0.5.3/scalpel-discovery --single \
--bam  ~/data/project/myGenome/fastq/bamFiles/jmzeng.filter.rmdup.bam \
--ref ~/reference/genome/hg19/hg19.fa \
--bed ${arr[0]}:${arr[1]}-${arr[2]}  \
--window 600 --numprocs 5  --dir ${arr[0]}
 
end=`date +%s`
runtime=$((end-start))
echo "Runtime for ${arr[0]}:${arr[1]}-${arr[2]} was $runtime"
 
done
```
## 人类有多少基因是有GO数据库注释信息的呢？
直接从NCBI的ftp服务器里面下载文件 ftp://ftp-trace.ncbi.nih.gov/gene/DATA/gene2go 这个是有go注释的基因情况，然后写脚本统计多少基因是有GO功能的。

```shell
cat gene2go | grep -v '#' | cut -f2 | sort -u | wc -l
210607
```

- 解析：
1. 用`grep -v`选取不带`#`的行，即不选第一行
2. 用`cut -f2`选取第二列，即 GeneID 
3. 用`sort -u`排序和去掉相同的GeneID
4. 用`wc -l`统计行数，即GeneID数

