# 第一个Python程序

早就计划着赶紧把Python学一下，无奈因为各种琐事（或是我的拖延症），一直拖啊拖。所以看到当当打折，就赶紧买了几本书压压惊。因为之前已经有一些编程基础，我也不太有耐心看一些太过入门的书，所以就从《Python生物信息学数据管理》这本书开始看，这本书每个章节都从一个生物学的实际问题出发，通过问题来逐渐学习Python，比较符合我的风格。但是这本书上有些基础内容又讲的不够详细，我就翻一翻《Python从入门到实践》，互相取长补短，或者直接Google，这样往往可以学到的更多。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180430/Fc2fKimlI6.jpg?imageslim)

从这周开始，我开始更新我的Python笔记，这个笔记主要是按照《Python生物信息学数据管理》这本书，当然内容也不完全相同，也夹杂了自己的逻辑，这本书每个章节后面都有课后题，我也会对于我认为有意义的题目给出我的答案仅供参考，也许还有许多不足的地方，请各位大神多多指教。

## 本文内容

1. 学会用`#`来注释
2. 学会将字符串储存到变量中
3. 学会使用`for`循环，和`print`打印结果

## 问题

计算胰岛素序列中的氨基酸频率：

```python
# insulin [Homo sapiens] GI:386828
insulin = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
for aa in "ACDEFGHIKLMNPQRSTVWY":
    number = insulin.count(aa)
    print(aa, number)
```

## 程序是如何工作的

- `#`开头的是注释
- 第二行`insulin = "MAL..."`则将蛋白质序列储存在名为`insulin`的变量中；
- `for aa in "ACD...":`表示遍历二十个氨基酸的循环，需要注意的是`for`语句末尾的冒号告诉Python，下一行是循环的第一行，而且在Python中用缩进（注意！是**四个空格**，若是使用编辑器比如Atom，则把Tab Length改成4）来判断代码之间的关系，在前面的示例中，因为后面两行的缩进，所以表示了这两行包括在`for`循环中；当`for`循环使用的序列是字符串时，循环内的代码会一个个遍历整个字符串；
- `number = insulin.count(aa)`这一行使用`count`计算一个字符在一段文本中出现个数；
- `print(aa, number)`将氨基酸和个数打印出来，在`Python3`中使用的是`print`函数。

### 注释

说到注释，必须得改改当初学习Perl的坏习惯了，以前写程序都懒得写注释，感觉都是一次性的东西，用完就丢。。。

前几天让师兄帮我改改Python的代码，因为没有注释，差点逻辑有些混乱了，被师兄好好教育了下。。。

### 字符串变量

在Python字符串中需要有单引号`'abc'`、双引号`"abc"`或三引号`'''abc'''，"""abc"""`进行封闭。

#### 这三者的联系和区别

单引号和双引号都可以用来表示一个字符串，比如：

```python
insulin = "MALWMR"
insulin = 'MALWMR'
```

在上面的例子中，是没有任何区别的（当然，也都可以使用`\n`、`\t`等等转义字符）。但是，当用单引号括起来的字符串中出现单引号时就需要用`\`转义了，这时候，要是用双引号括起来，里面的单引号就不用转义~

总之，就是当你用单引号`' '`定义字符串的时候，它就会认为你字符串里面的双引号`" "`是普通字符，从而不需要转义。反之当你用双引号定义字符串的时候，就会认为你字符串里面的单引号是普通字符无需转义。

3个单引号及3个双引号的用处：

通常情况下我们用单引号或者双引号定义一个字符串的时候只能把字符串连在一起写成一行，如果非要写成多行，就得在每一行后面加一个`\`表示连字符。比如：

```python
str1 = "List of name:\  
        Hua Li\  
        Chao Deng"  
```
而且即使你这样写也不能得到期望的输出，实际上输出是下面这样的：

```python
>>> str1 = "List of name:\  
...         Hua Li\  
...         Chao Deng"  
>>> print(str1)  
List of name:        Hua Li        Chao Deng  
```

那么该如何得到我们期望的一行一个名字的输出格式呢？这就是3个引号的作用了：

```python
>>> str1 = """List of name: 
... Hua Li 
... Chao Deng 
... """  
>>> print(str1)  
List of name:  
Hua Li  
Chao Deng  
```

#### 索引和切片

通过方括号中的数字所以可以方便地提取字符串中的内容，索引从0开始，负的索引则表示从末尾开始计数；当在方括号中插入`:`可以提取字符串的一部分（切片）。

```python
insulin = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
print(insulin[0],insulin[-1]) #打印第一个字符和随后一个字符
print(insulin[0:3]) #打印第一个到第三个字符
print(insulin[0:]) #打印整个字符串
```

#### 字符串算术

加号将两个字符串串联，乘号则表示重复字符串几次（只能是整数）：

```python
>>> 'protein'+' '+'degradation'
'protein degradation'
>>> 'protein' * 2
'proteinprotein'
```

#### 计算字符串长度

`len()`函数返回字符串的长度：

```python
>>> len('protein')
7
```

#### 字符计数

`s.count()`可以计算字符串中某字符的个数：

```python
>>> 'protein'.count('r') # 计算r的个数
1
```

- 这里先简单介绍下python中的函数（Function）和方法（Method），以后的文章再展开：
我的理解来说，函数就是可以直接调用来实现一些功能的，比如上面的`len()`；而方法呢就是上面的`s.count()`，同样封装了独立的功能，但是方法是需要通过对象来调用的，表示针对这个对象要做的操作。

另外，在python中，对任何一个对象，包、函数或类，都可以用`help(ObjectName)`的形式，显示其函数的帮助。用`dir(ObjectName)`的方法，可以显示对象的属性和方法列表。

## 自测题

1. 端粒蛋白质序列中氨基酸出现的频率，哪种氨基酸最频繁？
2. DNA序列中碱基的出现频率。
3. 一次一个残基地打印出氨基酸序列。

## 我的答案

```python
# 1
Telomerase_seq = "MPRAPRCRAVRSLLRSHYREVLPLATFVRRLGPQGWRLVQRGDPAAFRALVAQCLVCVPWDARPPPAAPSFRQVSCLKELVARVLQRLCERGAKNVLAFGFALLDGARGGPPEAFTTSVRSYLPNTVTDALRGSGAWGLLLRRVGDDVLVHLLARCALFVLVAPSCAYQVCGPPLYQLGAATQARPPPHASGPRRRLGCERAWNHSVREAGVPLGLPAPGARRGGSASRSLPLPKRPRRGAAPEPERTPVGQGSWAHPGRTRGPSDRGFCVVSPARPAEEATSLEGALSGTRHSHPSVGRQHHAGPPSTSRPPRPWDTPCPPVYAETKHFLYSSGDKEQLRPSFLLSSLRPSLTGARRLVETIFLGSRPWMPGTPRRLPRLPQRYWQMRPLFLELLGNHAQCPYGVLLKTHCPLRAAVTPAAGVCAREKPQGSVAAPEEEDTDPRRLVQLLRQHSSPWQVYGFVRACLRRLVPPGLWGSRHNERRFLRNTKKFISLGKHAKLSLQELTWKMSVRDCAWLRRSPGVGCVPAAEHRLREEILAKFLHWLMSVYVVELLRSFFYVTETTFQKNRLFFYRKSVWSKLQSIGIRQHLKRVQLRELSEAEVRQHREARPALLTSRLRFIPKPDGLRPIVNMDYVVGARTFRREKRAERLTSRVKALFSVLNYERARRPGLLGASVLGLDDIHRAWRTFVLRVRAQDPPPELYFVKVDVTGAYDTIPQDRLTEVIASIIKPQNTYCVRRYAVVQKAAHGHVRKAFKSHVSTLTDLQPYMRQFVAHLQETSPLRDAVVIEQSSSLNEASSGLFDVFLRFMCHHAVRIRGKSYVQCQGIPQGSILSTLLCSLCYGDMENKLFAGIRRDGLLLRLVDDFLLVTPHLTHAKTFLRTLVRGVPEYGCVVNLRKTVVNFPVEDEALGGTAFVQMPAHGLFPWCGLLLDTRTLEVQSDYSSYARTSIRASLTFNRGFKAGRNMRRKLFGVLRLKCHSLFLDLQVNSLQTVCTNIYKILLLQAYRFHACVLQLPFHQQVWKNPTFFLRVISDTASLCYSILKAKNAGMSLGAKGAGPLPSEAVQWLCHQAFLLKLTRHRVTYVPLLGSLRTAQTQLSRKLPGTTLTALEAAANPALPSDFKTILD"

count_aa = {}
max_count = 0

for aa in "ACDEFGHIKLMNPQRSTVWY":
    number = Telomerase_seq.count(aa)
    count_aa[aa] = (number)
    if number > max_count:
        max_count = number
for key,value in count_aa.items():
    if value == max_count:
        print(key)

# 2
seq = "AAAACCCGGT"

for nt in "ATCG":
    number = seq.count(nt)
    print(nt,number)

# 3
insulin = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
for i in insulin:
    print(i)

insulin = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
for i in range(1,len(insulin)-1):
    print(insulin[0:i])
```

aaa


