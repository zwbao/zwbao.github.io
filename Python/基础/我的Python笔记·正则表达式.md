## 本文内容

- 使用正则表达式来表达共有序列
- 使用正则表达式在字符串中搜索子字符串

## 问题

在蛋白质序列中搜索磷酸化motif，并返回第一个出现的motif


```python
import re
seq = "VSVLTMFRYAGWLDRLYMLVGTQLAAIIHGVALPLMMLI"
pattern = re.compile('[ST]Q')
match = pattern.search(seq)
if match:
    print (match.start(),match.end())
    print (match.group())
else:
    print ("no match")
```

    21 23
    TQ
    

## 程序是如何工作的

- 程序第一行导入了模块`re`；
- `pattern = re.compile('[ST]Q')`设定一个正则表达对象来匹配`SQ`或`TQ`；
- 若序列中包含pattern子串，则返回Match对象，否则返回None。

### 编译正则表达式

- 在Python中使用正则表达式需要使用模块`re`；
- `compile()`编译字符串并且转换成正则表达对象。当然，也可以传递参数到`compile()`中。

### 模式匹配

- `search()`函数扫描字符串并寻找正则表达式第一次匹配的位置，需要注意的是`search()`方法返回的是`Match object`而不是直接返回字符串


```python
import re
seq = "ATAGAGAGGTAGAGTAAG"
pattern = re.compile('TA[AG]')
match = pattern.search(seq)
match
```




    <_sre.SRE_Match object; span=(1, 4), match='TAG'>



- `match.group()` 返回第一个匹配的字符串；
- `match.span()` 返回第一个匹配结果的起点和终点；
- `match.start()` 返回第一个匹配结果的起始位置；
- `match.end()` 返回匹配结果的结束位置。


```python
import re
seq = "ATAGAGAGGTAGAGTAAG"
pattern = re.compile('TA[AG]')
match = pattern.search(seq)
print(match.group(),match.span(),match.start(),match.end())
```

    TAG (1, 4) 1 4
    

#### 如果我想要找所有的匹配结果呢？

- `findall()`可返回包含所有匹配的一个字符串列表；
- `finditer()`也可以返回所有匹配对象并且返回迭代器形式，即可以用`for`循环进行遍历。


```python
all = pattern.findall(seq)
all
```




    ['TAG', 'TAG', 'TAA']




```python
iter = pattern.finditer(seq)
iter
```




    <callable_iterator at 0x234fecd0780>




```python
for i in iter:
    print(i.group())
    print(i.span())
    print(i.start())
    print(i.end())
```

    TAG
    (1, 4)
    1
    4
    TAG
    (9, 12)
    9
    12
    TAA
    (14, 17)
    14
    17
    

### 分组

有时候会将一个正则表达式分为若干个子组，来匹配不同的部分，比如想要知道`.`匹配了什么氨基酸，可以将她用圆括号括起来以创建一个组，然后使用`group()`方法得到相匹配的氨基酸类型。
- `group()`方法不填写参数或参数为0，则返回完全匹配的字符串；子组从1开始自左向右编号。



```python
import re
seq = "VSVLTMFRYSTALDRLYMLVGTQLAAIIHGVALPLMMLI"
pattern = re.compile('R(.)([ST])([^p])')
match = pattern.search(seq)
print(match.group())
print(match.group(1))
print(match.group(2))
print(match.group(3))
```

    RYST
    Y
    S
    T
    

- 也可以向`group()`方法中传递多个参数，得到各个子组的元组。


```python
print(match.group(3,2,1))
```

    ('T', 'S', 'Y')
    

- `groups()`方法返回一个包含所有与子组相关的元组。


```python
print(match.groups())
```

    ('Y', 'S', 'T')
    

### 修改字符串

re模块提供了三种修改字符串的方法：
- `split(s)`
- `sub(pattern, repl, string, count=0, flags=0)`
- `subn(r,s,[c])`

`split(s)`方法将分割符合正则表达式的字符串，产生一个列表。在下面的例子中，实现了一个字符串在`|`处分割，又因为`|`是元字符所以需要加反斜杠转义。


```python
import re
separator = re.compile('\|')
annotation = 'ATOM:CA|RES:ALA|CHAIN:B|NUMRES:166'
columes = separator.split(annotation)
print(columes)
```

    ['ATOM:CA', 'RES:ALA', 'CHAIN:B', 'NUMRES:166']
    

`sub(pattern, repl, string, count=0, flags=0)`方法：
- pattern : 正则中的模式字符串；
- repl : 替换的字符串，也可为一个函数；
- string : 要被查找替换的原始字符串；
- count : 模式匹配后替换的最大次数，默认 0 表示替换所有的匹配。

