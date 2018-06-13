# 我的Python笔记·模块化编程（一）

## 如何自定义和调用函数

### 定义一个函数

```python
def my_function(arg1, arg2, ...):
    '''documentation'''
    <instructions>
    return value1, value2, ...
```

`my_function`是函数名，函数名可以是除了保留字以外的任意名字。

- 查看保留字：
```python
import keyword
keyword.kwlist
```

`(arg1, arg2, ...)`为函数的参数，是可选的（有些函数可以不需要参数）。
`'''documentation'''`是用三引号引起来的可选注释。
`return`就是将结果返回到调用的地方，并把程序的控制权一起返回。
`return value1, value2, ...`是指函数返回的值。

例如，一个计算两数之和的函数可以这么写：

```python
def add(num1,num2):
    '''calculates the sum of two numbers'''
    result = num1 + num2
    return result
```

`add`函数取两个参数相加并返回结果。

### 调用函数

调用一个函数时，必须在函数后加圆括号。如果函数需要参数，则需要在括号中填入参数。比如，调用前面的`add`函数：

```python
result = add(12,8)
print(result)
```

### Python函数的十大注意事项

- `def`为定义函数的语句标识；
- 必须使用圆括号来定义和调用函数；
- 函数代码块以冒号开头，后面紧跟缩进；
- 最后一个缩进是函数定义结束的标志；
- 传递给函数的参数序列是一个元组，函数也会以元组的形式返回多个值；
- 可以在函数内定义变量；
- `return`会退出函数，也可以选择将值传回调用方；
- `return`可以无返回值，而函数也可以无返回语句，这两种情况下，默认的返回值为`None`；
- 可以在函数中插入用三引号引起来的文档字符串。这些字符串会在函数调用时被忽略，但可以利用函数对象的`__doc__`属性来检索；
- 函数内的变量为局部变量，而不再脚本或模块的全局空间内。当一个函数被调用时，首先会在函数命名空间内搜索函数内的变量名，如果在函数和体内未找到对象名称，然后就会在脚本或模块的全局空间内进行搜索。

### lambda函数（又称匿名函数）

lambda表达式，通常是在需要一个函数，但是又不想费神去命名一个函数的场合下使用，也就是指匿名函数。lambda所表示的匿名函数的内容应该是很简单的，如果复杂的话，不如就重新定义一个函数吧。

比如上面的`add`函数用lambda写就是这样：

```python
result = lambda x,y : x + y
print(result(12,8))
```

lambda函数不包括`return`，而是包括一个表达式，而且总是会返回该表达式的值。且可以有任意多个参数，并返回单个表达式的值。

当然，也不是任何情况下lambda函数都要比常规函数更清晰明了，Python之禅中有这么一句话：Explicit is better than implicit，就是说那种方式更清晰就用哪一种方式，不要盲目的都使用lambda表达式。

### `struct`模块

`struct`模块是Python的内置模块，可以在自定义格式的基础上将字符串转化为元组。`struct`方法`pack(format, v1, v2, ...)`可以根据字符串格式返回由`v1, v2, ...`值压缩成的单一字符串。例如：

```python
import struct
format = '2s3s'
a = struct.pack(format, '10','100')
a
```

表示两个字符的字符串后面跟一个三个字符的字符串。

需要注意的是上面的例子在Python 2.x中可以编译成功，因为`s`表示的是字符串：
![struct2.x](_v_images/_struct2x_1527427559_24083.png)

而在Python 3.x中则会报错：
```python
---------------------------------------------------------------------------
error                                     Traceback (most recent call last)
<ipython-input-6-5811fc7b09a8> in <module>()
      1 import struct
      2 format = '2s3s'
----> 3 a = struct.pack(format, '10','100')
      4 a

error: argument for 's' must be a bytes object
```

这是因为在Python 3.x中，格式化字符串的`s`参数，在Python中的类型是bytes类型。所以把字符串的地方转为字节类型即可解决：

```python
import struct
format = '2s3s'
a = struct.pack(format, b'10',b'100')
a
```

![struct](_v_images/_struct_1527427484_5431.png)

在表中可以看到，`c`,`s`,`p`的Python type都是bytes类型，所以在使用这些类型的时候，要将需要pack的字符串写成bytes型的。

方法`unpack(format, string)`按照根据`format`编码格式将`string`解压为元组，字符串中包含的字符与格式字符数必须相等。

```python
import struct
format = '2s3s'
line = '10100'
col = struct.unpack(format, line)
col
```

同样的，在Python 3.x需要写成：

```python
import struct
format = '2s3s'
line = '10100'
col = struct.unpack(format, bytes(line.encode('utf-8')))
col
```

这里除了要把字符串的地方转为字节类型,还要先转成`utf-8`的编码，否则报错`string argument without an encoding`。

方法`calcsize(fmt)`返回给定格式化字符串的总字符数。

```python
import struct
format = '2s3s'
struct.calcsize(format) # 5
```

## 问题

在PDB数据库中下载1tld.pdb文件。已知胰蛋白酶活性部位为 Asp 102、His 57和Ser 195。识别pdb文件中的这三个残基并提取其坐标保存到文件中。

PDB文件的格式：
![pdb](https://wx1.sinaimg.cn/mw690/c5d7b0ebgy1frqu7z4bu8j20r80sln0l.jpg)

```python
pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s10s2s3s'
```

### 答案

代码使用两个函数来完成这项任务：一个从PDB文件中解析单行，调用`struct`模块来将单行进行切分，得到所需的数据（残基名称，序号，坐标）输出为元组；另一个则用于处理整个文件，进行条件判断和写入。

正如上一部分所说的在Python3中`s`的Python type是bytes类型，所以当进行条件判断时需要进行`decode`解码。

![decode](https://images2017.cnblogs.com/blog/764761/201708/764761-20170828175555827-950234753.png)

- Python3答案

```python
import struct
pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s10s2s3s'

def parse_atom_line(line):
    '''return an ATOM line parsed to a tuple'''
    tmp = struct.unpack(pdb_format, bytes(line.encode('utf-8')))
    atom = tmp[3].strip()
    res_type = tmp[5].strip()
    res_num = tmp[8].strip()
    chain = tmp[7].strip()
    x = float(tmp[11].strip())
    y = float(tmp[12].strip())
    z = float(tmp[13].strip())
    return chain, res_type, res_num, atom, x, y, z

def main(pdb_file, residues, outfile):
    '''writes residues from a PDB file to an output file'''
    pdb = open(pdb_file)
    outfile = open(outfile, "w")
    for line in pdb:
        if line.startswith('ATOM'):
            res_data = parse_atom_line(line)
            res2type = bytes.decode(res_data[1])
            res2num = bytes.decode(res_data[2])
            for aa,num in residues:
                if res2type == aa and res2num == num:
                    outfile.write(line)
    outfile.close()

residues = [('ASP', '102'), ('HIS', '57'), ('SER', '195')]
main("1TLD.pdb", residues, "trypsin_triad.pdb")
```

- Python2答案

```python
import struct
pdb_format = '6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s10s2s3s'

def parse_atom_line(line):
    '''return an ATOM line parsed to a tuple'''
    tmp = struct.unpack(pdb_format, line)
    atom = tmp[3].strip()
    res_type = tmp[5].strip()
    res_num = tmp[8].strip()
    chain = tmp[7].strip()
    x = float(tmp[11].strip())
    y = float(tmp[12].strip())
    z = float(tmp[13].strip())
    return chain, res_type, res_num, atom, x, y, z

def main(pdb_file, residues, outfile):
    '''writes residues from a PDB file to an output file'''
    pdb = open(pdb_file)
    outfile = open(outfile, "w")
    for line in pdb:
        if line.startswith('ATOM'):
            res_data = parse_atom_line(line)
            for aa,num in residues:
                if res_data[1] == aa and res_data[2] == num:
                    outfile.write(line)
    outfile.close()

residues = [('ASP', '102'), ('HIS', '57'), ('SER', '195')]
main("1TLD.pdb", residues, "trypsin_triad.pdb")
```

aaa
