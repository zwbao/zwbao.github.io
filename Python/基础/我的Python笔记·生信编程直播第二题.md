# 我的Python笔记·生信编程直播第二题

## 本文内容

1. 读取和写入文件
2. 字符串处理

## 问题

这次的题目来自[生信编程直播第二题-hg19基因组序列的一些探究](https://vip.biotrainee.com/d/265 )，上一次这题我没做出来，这次正好补上自己的答案~

统计每条染色体长度，ATCGN的含量：

- 测试数据
```
>chr_1
ATCGTCGaaAATGAANccNNttGTA
AGGTCTNAAccAAttGggG
>chr_2
ATCGAATGATCGANNNGccTA
AGGTCTNAAAAGG
>chr_3
ATCGTCGANNNGTAATggGA
AGGTCTNAAAAGG
>chr_4
ATCGTCaaaGANNAATGANGgggTA
```

- 我的答案
```python
InputFile = open('test.fa','r')
seq = ''

print("chromsome\tA\tT\tC\tG\tN\tlength")

    line = line.upper()
    if line[0] == '>' and seq == '':
        header = line.strip()
    elif line[0] != '>':
        seq = seq + line.strip()
    elif line[0] == '>' and seq != '':
        print(header,end='\t')
        length = len(seq)
        for nt in "ATCGN":
            if nt == "N":
                print(str(seq.count(nt))+"\t"+str(length),end='\n')
            else:
                print(seq.count(nt),end='\t')
        seq = ''
        header = line.strip()
# 注意最后一条数据
print(header,end='\t')
length = len(seq)
for nt in "ATCGN":
    if nt == "N":
        print(str(seq.count(nt))+"\t"+str(length),end='\n')
    else:
        print(seq.count(nt),end='\t')
        
InputFile.close()
```

## 程序是如何工作的

- 首先从fasta文件中读入测试数据，通过`for`循环来逐行读取数据，因为测试数据包含大小写，所以需要用`upper`方法使读入的行都变成大写，方便之后的判断；
- fasta格式以`>`开头为第一行，注明了描述信息，换行之后是序列信息，就可以根据开头有没有`>`来判断是否是描述行，同时用变量`seq`来储存序列信息；
- 如果一行是以`>`开头，`seq`为0则说明读到新的序列，则保存描述信息；如果不是以`>`开头，则说明改行是序列信息，将这一行保存到`seq`中；如果这一行是以`>`开头，但是`seq`中已经保存了数据，说明程序已经读到下一行数据的开头，这时候就可以把上一条序列的信息输出；
- 需要注意的是，当程序读完最后一条序列的所有信息后，因为没有下一个`>`来触发进行数据处理的操作，所以需要在循环之外进行处理。

### 如何读取和写入文件

#### 读取数据
```python
InputFile = open('test.fa')
lines = InputFile.readlines()
InputFile.close()
```

该程序做了三件事：

1. 打开了一个文本文件。在单引号或双引号之间指定文件名，该文件名假设与Python程序位于同一目录，否则需要在文件名之前添加完整路径。
2. 从文件中读取数据。`readlines()`读取文件中的所有内容，按分隔符逐行储存字符串（**包括每行最后的换行符**），最后会返回一个字符串列表。而`read()`则是读取整个文件，作为一个单一的字符串。
3. 关闭文本文件。程序结束时Python会立即自动关闭文件，但如果没有关闭文件就想第二次打开，就会出现问题。

当然也可以：

```python
for line in open('test.fa'):
    line = line.strip()
```

第一行会打开文件进行读取，并且运行`for`循环遍历每一行，`strip()`函数则将删除开头或结尾处的空字符或换行符`\n`，在上面的答案中如果不使用`strip()`，在用`len()`统计长度时就会多出两个字符（就是换行符`\n`）。

#### 写入数据
```python
OutputFile = open('results.txt','w')
OutputFile.write('A'+'\t'+'13\n')
OutputFile.close()
```
这段代码执行了以下操作：
1. 打开一个即将要写入数据的文本文件，与用于读取的`open()`的区别在于字母`'w'`，表示只能用于写入（write），当多次写入时只会保留最后一次的写入（覆盖），若需要在文本中追加内容，则可以使用`'a'`。
2. 将字符串写入文件。`write()`函数只能写入字符串数据，所以需要将写入的数据转化为字符串，稍后将解释如何将字符转化为字符串，而且`write()`函数不包含换行，需要在末尾添加`\n`。
3. 关闭使用后的文件。如果忘记关闭文件，当你打开文件时会发现这个文本是**空的**，因为你要写入的东西还在计算机内存中。只有达到一定字符的字符串才会被保存到文件中。但是关闭文件，无论写入多少个字符都会被保存。

#### 将数字转化为文本

```python
print("A"+13)
```

当字符串和数字同时输出时，程序会报错，需要统一成数字或者字符，因为字符串不能转化为数字，所以往往把数字转化为字符串。

```python
print("A"+str(13))
```

但是，使用`str()`函数来转换数字有很大的弊端：文本文件中的数字是未格式化的，无法对其，尤其是浮点数。`str()`的替代方法之一是使用字符串格式化，可以在数值转化为字符串时使用百分号指示要分配给整数的位数。

```python
'Result:%3i' % (17)
```

`%3i`表示字符串被格式化为三位的整数，整数的实际值在结尾的括号中。相同地，可以插入浮点数字符串`%x.yf`，`x`是总字符数（包括小数点），`y`是小数位数。

#### elif

`elif`意为else if， 当且仅当前述的所有if/elif的语句都不成立时，用来检索条件。所以需要注意if/elif语句的顺序，只要一个条件得到满足，所有后续在相同if/elif/else块中的语句将被Python忽略。此外，只能在一个if/elif/else块钟使用一次else语句，不满足其他预设条件时将会执行它。

## 自测题

1. 读取和写入多序列FASTA文件，并将每条记录（序列+标题）写入不同的文件中。
2. 读取和过滤FASTA文件，仅当起始氨基酸为甲硫氨酸（M）且含有至少两个色氨酸（W）时，将记录写入新文件。
3. 计算在FASTA格式中多个DNA序列的核苷酸频率。

## 我的答案

```python
# 1
InputFile = open('4_protein.fa','r')
seq = ''
n = 0

for line in InputFile:
    if line[0] == '>':
        n += 1
        OutputFile = open(str(n)+'.fa','a')
        OutputFile.write(line)
        OutputFile.close()
    else:
        OutputFile = open(str(n)+'.fa','a')
        OutputFile.write(line)
        OutputFile.close()

# 2
InputFile = open('4_protein.fa','r')
seq = ''

for line in InputFile:
    if line[0] == '>' and seq == '':
        header = line
    elif line[0] != '>':
        seq = seq + line
    elif line[0] == '>' and seq != '':
        if seq[0] == 'M' and seq.count('W') >= 2:
            OutputFile = open(header[1:9].strip()+'.fa','w')
            OutputFile.write(header + seq)
            OutputFile.close()
            print(header + seq)
        seq = ''
        header = line
# take care of the very last record of the input file
if seq[0] == 'M' and seq.count('W') >= 2:
    OutputFile = open(header[1:9].strip()+'.fa','w')
    OutputFile.write(header + seq)
    OutputFile.close()
    print(header + seq)

# 3
InputFile = open('4_DNA.fa','r')
seq = ''

for line in InputFile:
    if line[0] == '>' and seq == '':
        header = line
    elif line[0] != '>':
        seq = seq + line
    elif line[0] == '>' and seq != '':
        print(header)
        length = float(len(seq))
        for nt in "ATCG":
            number = seq.count(nt)
            print(nt,number/length)
        seq = ''
        header = line
print( header)
for nt in "ATCG":
    number = seq.count(nt)
    print(nt,number/length)
```
aaa