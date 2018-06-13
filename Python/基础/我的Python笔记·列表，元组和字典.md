# 我的Python笔记·列表，元组和字典

## 本文内容

1. 列表和元组
2. 字典

## 问题

写一个程序，计算输入核苷酸序列的终止密码子和起始密码子数量。该程序需要打印两个元素：起始密码子数目和终止密码子数目。

```python
codon_table = {
    'TAA':'STOP',
    'TAG':'STOP',
    'TGA':'STOP',
    'ATG':'Start'}
rna = ''
for line in open('sequence.fasta'):
    if not line.startswith('>'):
        rna = rna + line.strip()
for frame in range(3):
    print('Reading frame'+str(frame + 1))
    count_STOP = 0
    count_START = 0
    for i in range(frame,len(rna),3):
        codon = rna[i:i + 3]
        if codon in codon_table:
            if codon_table[codon] == 'STOP':
                count_STOP += 1
            if codon_table[codon] == 'Start':
                count_START += 1
    print('START',count_START)
    print('STOP',count_STOP)
```

## 程序是如何工作的

- 首先需要建立密码子和其是否为终止密码子或起始密码子的映射关系，即字典；
- 读入`fasta`文件，将核苷酸序列赋值给变量`rna`;
- 有三种判定密码子的方式，所以需要用`for`循环进行遍历；
- 最后用`if`来判定读到的密码子是否为终止密码子或起始密码子并进行计数。

### 列表和元组

#### 列表
序列是Python中最基本的数据结构。序列中的每个元素都分配一个数字，即它的位置，或索引，第一个索引是0，第二个索引是1，依此类推。

序列都可以进行的操作包括索引，切片，加，乘，检查成员。
此外，Python已经内置确定序列的长度以及确定最大和最小的元素的方法。

创建一个列表，只要把逗号分隔的不同的数据项使用方括号括起来即可。如下所示：

```python
codon_list = ["TAA", "TAG", "TGA"]
```

##### 访问列表中的值
使用下标索引来访问列表中的值，同样也可以使用方括号的形式截取字符，如下所示：

```python
codon_list = ["TAA", "TAG", "TGA"]
print ("codon_list[0]: ", codon_list[0])
print ("codon_list[0:2]: ", codon_list[0:2])
```

##### 更新列表
可以对列表的数据项进行修改或更新，也可以使用append()方法来添加列表项，如下所示：

```python
codon_list = [] # 定义一个空列表
codon_list.append("TAA") # 使用 append() 添加元素
codon_list.append("TGA")
print(codon_list)
```

##### 删除列表元素
可以使用 del 语句来删除列表的元素，如下：

```python
codon_list = []
codon_list.append("TAA")
codon_list.append("TGA")
del codon_list[0]
print(codon_list)
```

##### 计算数值列表

Python有一系列用于数字和字符串列表的内置函数：
- `len()`返回列表长度，即列表中的项目数；
- `max()`返回列表的最大元素；
- `min()`返回最小元素；
- `sum()`将所有元素加起来。

#### 元组

Python的元组与列表类似，不同之处在于元组的元素不能修改。元组使用小括号，列表使用方括号，如：`(a, b, c)`，或简单地列出序列数据项并用逗号隔开：`a, b, c`。索引和切片的操作在元组上也适用。

### 字典
开头定义的`codon_table`对象是一个字典。字典为“键(key)：值(value)”对的对象的无序集合，用花括号括起来：`{'TAA':'STOP'}`。需要注意的是，字典是不可变对象（键）映射任意对象（值）的结构。不可变对象可以是数字，字符串和元组。列表和字典本身不能用作字典键，但可以作为值。键和他的值之间用冒号分割，键值对之间用逗号分隔。键一般是唯一的，如果重复最后的一个键值对会替换前面的，值不需要唯一。

#### 访问字典里的值
```python
codon_table = {
    'TAA':'STOP',
    'TAG':'STOP',
    'TGA':'STOP',
    'ATG':'Start'}
print ("codon_table['TAA']:", codon_table['TAA'])
```

#### 修改字典
向字典添加新内容的方法是增加新的键/值对，修改或删除已有键/值对：

```python
codon_table = {
    'TAA':'STOP',
    'TAG':'STOP',
    'TGA':'STOP'}
codon_table['ATG'] = "Start"; # Add new entry
print ("codon_table['ATG']:", codon_table['ATG'])
```

#### 删除字典元素
能删单一的元素也能清空字典。显示删除一个字典用del命令：

```python
codon_table = {
    'TAA':'STOP',
    'TAG':'STOP',
    'TGA':'STOP',
    'ATG':'Start'}

del codon_table['ATG'] # 删除键是'ATG'的条目
codon_table.clear() # 清空词典所有条目
del codon_table # 删除词典
```

#### 字典键的特性
- 不允许同一个键出现两次。创建时如果同一个键被赋值两次，后一个值会被记住；
- 键必须不可变，所以可以用数字，字符串或元组充当，所以用列表就不行。

### 字典搜索

在这次的问题中需要遍历RNA序列三次：一次从序列第一位开始，一次从第二位开始，一是第三位。代码如下：

```python
for i in range(frame,len(rna),3):
    codon = rna[i:i + 3]
    if codon in codon_table:
        if codon_table[codon] == 'STOP':
            count_STOP += 1
        if codon_table[codon] == 'Start':
            count_START += 1
```

以三个字母为单位进行扫描，每三个字母都会有从`codon_table`字典中判断是否为终止密码子或起始密码子。

## 自测题

写一个基于序列的二级结构元素预测的程序。

**提示：**使用如下偏好表，其中`pref_H`是α螺旋，`pref_E`是β折叠：

```python
pref_H = {
    'A':1.45,
    'C':0.77,
    'D':0.98,
    'E':1.53,
    'F':1.12,
    'G':0.53,
    'H':1.24,
    'I':1.00,
    'K':1.07,
    'L':1.34,
    'M':1.20,
    'N':0.73,
    'P':0.59,
    'Q':1.17,
    'R':0.79,
    'S':0.79,
    'T':0.82,
    'V':1.14,
    'W':1.14,
    'Y':0.61
}

pref_E = {
    'A':0.97,
    'C':1.30,
    'D':0.80,
    'E':0.26,
    'F':1.28,
    'G':0.81,
    'H':0.71,
    'I':1.60,
    'K':0.74,
    'L':1.22,
    'M':1.67,
    'N':0.65,
    'P':0.62,
    'Q':1.23,
    'R':0.90,
    'S':0.72,
    'T':1.20,
    'V':1.65,
    'W':1.19,
    'Y':1.29
}
```

## 我的答案

```python
pref_H = {
    'A':1.45,
    'C':0.77,
    'D':0.98,
    'E':1.53,
    'F':1.12,
    'G':0.53,
    'H':1.24,
    'I':1.00,
    'K':1.07,
    'L':1.34,
    'M':1.20,
    'N':0.73,
    'P':0.59,
    'Q':1.17,
    'R':0.79,
    'S':0.79,
    'T':0.82,
    'V':1.14,
    'W':1.14,
    'Y':0.61
}

pref_E = {
    'A':0.97,
    'C':1.30,
    'D':0.80,
    'E':0.26,
    'F':1.28,
    'G':0.81,
    'H':0.71,
    'I':1.60,
    'K':0.74,
    'L':1.22,
    'M':1.67,
    'N':0.65,
    'P':0.62,
    'Q':1.23,
    'R':0.90,
    'S':0.72,
    'T':1.20,
    'V':1.65,
    'W':1.19,
    'Y':1.29
}

rna = ''
struc_rna = ''
for line in open('sequence.fasta'):
    if not line.startswith('>'):
        rna = rna + line.strip()
    for i in range(len(rna)):
        nt = rna[i]
        if pref_H[nt] >= 1 and pref_E[nt] < pref_H[nt]:
            struc_rna = struc_rna + 'H'
        elif pref_H[nt] >= 1 and pref_E[nt] > pref_H[nt]:
            struc_rna = struc_rna + 'E'
        else:
            struc_rna = struc_rna + 'L'
print(struc_rna)
```

aaa