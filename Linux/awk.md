# awk 奥克

NR number of record 第几行

NF numberof field 此行有多少列，注意：空格也是输入分隔符
`$NF`表示每一行的最后一列
`$(NF-1)`表示倒数第二列

一行 record
一列 field
$0：表示整一行

awk 默认的分隔符是 空格

字符要加上双引号

如果不加文件，会你输入一行执行一次该命令

如何自定义输入分隔符
在最前面定义全局变量
`BEGIN{FS=“，”}`
FS定义分隔符，以上代码表示逗号做分隔符

输入分隔符和输出分隔符不一样
定义输出分隔符：`BEGIN{OFS=“，”}`

同时改输入分隔符和输出分隔符：`BEGIN{FS=“，”；OFS=“，”}`
用分号分隔

当输入两个文件时，第二个文件会紧接着下一行输出

内部变量：FILENAME 等于文件名

打印第一列 `awk '{print $1}' FILENAME`

```
2.1 内置变量之字段变量$
	awk每次从文件中读入一行，并按照分隔符将其拆分为多个字段，并将所有字段以`$1,$`2,`$3,...的形式表示。$`0代表整行数据。
	默认分隔符是空格，如需要修改分隔符，可以使用-F选项，如 awk -F: '{print $1}' /etc/passwd 可以打印出当前系统中所有的用户名。
	修改分隔符也可以使用FS变量，见2.2。

2.2 内置变量之记录变量
    FS 读入文本信息的字段分隔符，默认空格
    RS 读入文本信息的行分隔符，默认换行
    OFS 输出格式的字段分隔符，默认空格
    ORS 输出格式的行分隔符，默认换行

    例如，同样是查看当前系统所有的用户名，可以使用如下命令：
    awk 'BEGIN {FS=":"}{print $1}' /etc/passwd 
        这个awk有两个pattern{action}语句：
        BEGIN {FS=":"}：BEGIN语句代表在读入文件之前进行的操作，也就是设定字段分隔符为":"
        {print $1}：这里的pattern为空，代表对所有行进行处理，每行均输出第一个字段

2.3  内置变量之数据变量
    NR：awk所处理的记录数（行数），如果有多个文件，则NR代表累计的处理字段数。
    NF：代表当前处理行的字段数
    FNR：与NR不同，如果有多个文件，则每个文件的记录数单独记录
    ARGV：数组变量，保存当前执行的awk命令行语句。
    FILENAME：awk命令所处理的文件名
    ENVIRON：当前shell环境变量及其值的数组，结构类似于python的字典。

    例如，
	awk 'BEGIN {print ARGV[0],ARGV[1]}' /etc/passwd 
		输出结果awk /etc/passwd
	awk 'END {print FILENAME}' /etc/passwd 
		输出结果 /etc/passwd
	awk 'BEGIN {print ENVIRON["PATH"]}' 
		输出当前shell的PATH环境变量
		注意：在awk中，字符串需要使用双引号“”，单引号特定用于包括pattern{action}。

2.4 用户自定义变量
    有两个方法：
    2.4.1 在脚本中定义
    awk 'BEGIN {VAR="test string"; print VAR}'	输出test string
    2.4.2 在命令行中定义
    使用-v选项，如
    awk -v VAR="test string"  'BEGIN {print VAR}'		输出test string
```

```
7.1 if-else语句
    使用格式：if (condition) {then-body;} [else {else-body}]
    例子
    awk -F: '{if (`$1=="root") print $`1,"Admin"; else print $1,"CommanUser"}' /etc/passwd
    输出当前系统中的用户名，如果用户名root则同时输出Admin，否则输出CommandUser。
7.2 while
    使用格式：while （condition）{statement1;statement2,...}
    例子
    awk -F: '{i=0; while (i<3) {print $i;i++}}' /etc/passwd
7.3 do-while
	使用格式：do {statement1;statement2,...} while （condition）
7.4 for
    1. 使用格式：
    for(variable assignment; condition; iteration process)  {statement1;statement2,...}
    例子：
    awk -F: '{for (i=0; i<3; i++) print $i}' /etc/passwd
    2. 也可以遍历数组
    for (i in array) {statement1;statement2,...}
	例子
	awk -F: '$NF!~/^$/{BASH[$NF]++} END{for (i in BASH) {printf "%-20s%d\n",i,BASH[i]}}' /etc/passwd 
	其中NF为每一行的字段数，`$NF代表最后一个字段，正则表达式/^$`/代表空行，所以这个pattern的意思是最后一个字段不为空的行，后面的action新建了一个数组BASH（数组见后述），数组的下标为满足条件行的最后一个字段（其实就是shell地址），自加运算符会对所有相同的shell地址进行运算。而END语句会在所有行执行完之后运行，其结果是格式化输出数组BASH的内容。
    输出结果为：
        /bin/sync           1
        /bin/bash           2
        /bin/false          11
        /usr/sbin/nologin   17
	这个命令还可以做一下延伸，将字段按照数目排序，再将printf接一个管道命令到sort中去，如下：
	awk -F: '`$NF!~/^$/{BASH[$NF]++} END{for (i in BASH) {printf "%-20s%d\n",i,BASH[i]  |"sort -k1"}}' /etc/passwd
	其含义为，将输出结果再按照第一列进行排序，输出结果为：
        /bin/bash           2
        /bin/false          11
        /bin/sync           1
        /usr/sbin/nologin   17	
7.5 case
    使用格式：
    switch (expression) {case Value or /RegExp/:statement1,statement1,... default :statement}
7.6 break 和 continue
	配合循环或case语句使用，同C语言。
7.7 next
    提前结束对本行文本的处理，并开始处理下一行
    例子
    awk -F: '{if (NR%2==0) print NR,$1;else next}' /etc/passwd
    代表只显示文件/etc/passwd偶数行的用户名。
```
