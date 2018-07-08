# NumPy入门（一）

## 简介

NumPy（Numerical Python）是 Python 中的一个线性代数库。对每一个数据科学或机器学习 Python 包而言，这都是一个非常重要的库，SciPy（Scientific Python）、Matplotlib（plotting library）、Scikit-learn 等都在一定程度上依赖 NumPy。

如果已经安装了`Anaconda`，那么你已经可以使用`NumPy`；也可以按照`NumPy`的官方网站（www.numpy.org）的教程进行安装。

## Start

### 导入NumPy

按照传统，通常用`np`作为别名来导入`NumPy`：

```python
import numpy as np
```

### 从Python列表创建数组：

可以用`np.array`从Python列表创建数组：


```python
# 整型数组：
np.array([1, 4, 2, 5, 3])
```




    array([1, 4, 2, 5, 3])


`NumPy`要求数组必须包含同一类型的数据，如果不匹配则会向上转换。比如整型会变成浮点型：


```python
np.array([3.14, 4, 2, 3])
```




    array([3.14, 4.  , 2.  , 3.  ])


可以用`dtype`关键字来设置数值类型：


```python
np.array([1, 2, 3, 4], dtype='float32')
```




    array([1., 2., 3., 4.], dtype=float32)


`NumPy`数组也可以被指定为多维：


```python
np.array([range(i , i+3) for i in [2, 4, 6]])
```




    array([[2, 3, 4],
           [4, 5, 6],
           [6, 7, 8]])


### 从头创建数组

`NumPy`有许多内置的方法从头创建数组，以下是几个示例：


```python
# 创建一个值都为0，长度为10的数组
np.zeros(10, dtype=int)
```




    array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])



```python
# 创建一个值都为1，3*5的浮点型数组
np.ones((3, 5),dtype=float)
```




    array([[1., 1., 1., 1., 1.],
           [1., 1., 1., 1., 1.],
           [1., 1., 1., 1., 1.]])



```python
# 创建一个值都为3.14，3*5的浮点型数组
np.full((3, 5), 3.14)
```




    array([[3.14, 3.14, 3.14, 3.14, 3.14],
           [3.14, 3.14, 3.14, 3.14, 3.14],
           [3.14, 3.14, 3.14, 3.14, 3.14]])



```python
# 创建一个值为线性序列，从0到20，步长为2的数组
np.arange(0, 20, 2)
```




    array([ 0,  2,  4,  6,  8, 10, 12, 14, 16, 18])



```python
# 创建一个有5个元素的数组，且这5个数均匀分配到0~1
np.linspace(0, 1, 5)
```




    array([0.  , 0.25, 0.5 , 0.75, 1.  ])



```python
# 创建一个3*3的，在0~1均匀分布的随机数组成的数组
np.random.random((3, 3))
```




    array([[0.36709154, 0.97609935, 0.53257059],
           [0.98030801, 0.75542908, 0.41482827],
           [0.91319197, 0.37346273, 0.55469639]])



```python
# 创建一个3*3的，均值为0，方差为1的正态分布的随机数数组
np.random.normal(0, 1, (3, 3))
```




    array([[-0.01152638,  0.88988951, -0.21171668],
           [ 0.08373947, -2.41123452,  0.50703231],
           [ 1.98705654, -1.61103069, -0.97971459]])



```python
# 创建一个3*3的、`[0, 10)`区间的随机整型数组
np.random.randint(0, 10, (3, 3))
```




    array([[6, 3, 8],
           [8, 8, 0],
           [1, 2, 1]])



```python
# 创建一个由3个整型数组成的未初始化的数组
# 数组的值是内存空间中的任意值
np.empty(3)
```




    array([9.80485013e-312, 0.00000000e+000, 9.76154329e-313])


### NumPy标准数据类型

`NumPy`数组包含同一类型的值，所以了解下这些数据类型是十分必要的：

| Data type	    | Description |
|---------------|-------------|
| ``bool_``     | Boolean (True or False) stored as a byte |
| ``int_``      | Default integer type (same as C ``long``; normally either ``int64`` or ``int32``)|
| ``intc``      | Identical to C ``int`` (normally ``int32`` or ``int64``)|
| ``intp``      | Integer used for indexing (same as C ``ssize_t``; normally either ``int32`` or ``int64``)|
| ``int8``      | Byte (-128 to 127)|
| ``int16``     | Integer (-32768 to 32767)|
| ``int32``     | Integer (-2147483648 to 2147483647)|
| ``int64``     | Integer (-9223372036854775808 to 9223372036854775807)|
| ``uint8``     | Unsigned integer (0 to 255)|
| ``uint16``    | Unsigned integer (0 to 65535)|
| ``uint32``    | Unsigned integer (0 to 4294967295)|
| ``uint64``    | Unsigned integer (0 to 18446744073709551615)|
| ``float_``    | Shorthand for ``float64``.|
| ``float16``   | Half precision float: sign bit, 5 bits exponent, 10 bits mantissa|
| ``float32``   | Single precision float: sign bit, 8 bits exponent, 23 bits mantissa|
| ``float64``   | Double precision float: sign bit, 11 bits exponent, 52 bits mantissa|
| ``complex_``  | Shorthand for ``complex128``.|
| ``complex64`` | Complex number, represented by two 32-bit floats|
| ``complex128``| Complex number, represented by two 64-bit floats|

### NumPy数组的属性

介绍一些有用的数组属性。首先用`NumPy`的随机数生成器设置一组种子值，确保每次程序运行时可以生成相同的随机数组，再定义三个随机数组：一个一维数组、一个二维数组和一个三维数组。


```python
np.random.seed(0) # 设置种子

x1 = np.random.randint(10, size=6)
x2 = np.random.randint(10, size=(3, 4))
x3 = np.random.randint(10, size=(3, 4, 5))
```

每个数组有`ndim`数组的维度、`shape`数组每个维度的大小、`size`数组的总大小属性和`dtype`数组的数据类型属性：


```python
print("x3 ndim: ", x3.ndim)
print("x3 shape:", x3.shape)
print("x3 size: ", x3.size)
print("dtype:", x3.dtype)
```

    x3 ndim:  3
    x3 shape: (3, 4, 5)
    x3 size:  60
    dtype: int32

其他属性包括数组元素的字节大小`itemsize`以及数组总字节大小`nbytes`：


```python
print("itemsize:", x3.itemsize, "bytes")
print("nbytes:", x3.nbytes, "bytes")
```

    itemsize: 4 bytes
    nbytes: 240 bytes

### 数组索引：获取三个元素

和Python列表一样，对于`NumPy`数组你也可以通过中括号来进行索引（从0开始计数）：


```python
x1
```




    array([5, 0, 3, 3, 7, 9])



```python
x1[0]
```




    5



```python
x1[4]
```




    7


同样可以用负值来进行末尾索引：


```python
x1[-1]
```




    9


在多维数组中，用逗号分隔的索引元组获取元素：


```python
x2
```




    array([[3, 5, 2, 4],
           [7, 6, 8, 8],
           [1, 6, 7, 7]])



```python
x2[0, 0]
```




    3



```python
x2[2, -1]
```




    7


也可以用以上索引方式修改元素值：


```python
x2[0, 0] = 12
x2
```




    array([[12,  5,  2,  4],
           [ 7,  6,  8,  8],
           [ 1,  6,  7,  7]])


注意，`NumPy`数组是固定类型的，如果把一个浮点值插入一个整型数组时，浮点数会截短变为整型。


```python
x1[0] = 3.14159
x1
```




    array([3, 0, 3, 3, 7, 9])


### 数组切片：获取子数组

为了获取数组`x`的切片，可以用以下方式：`x[start:stop:step]`。

#### 一维子数组


```python
x = np.arange(10)
x
```




    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])



```python
x[:5] # 前5个元素
```




    array([0, 1, 2, 3, 4])



```python
x[4:7] # 中间的子数组
```




    array([4, 5, 6])



```python
x[::2] # 每隔一个元素取
```




    array([0, 2, 4, 6, 8])



```python
x[1::2] # 每隔一个元素取，从索引1开始
```




    array([1, 3, 5, 7, 9])


当步长值为负时，`start`和`stop`参数是被交换的，所以这可以很方便地将数组逆序排列：


```python
x[::-1] # 逆序排列
```




    array([9, 8, 7, 6, 5, 4, 3, 2, 1, 0])



```python
x[5::-2] # 从索引5开始每隔一个元素逆序
```




    array([5, 3, 1])


#### 多维子数组


```python
x2
```




    array([[12,  5,  2,  4],
           [ 7,  6,  8,  8],
           [ 1,  6,  7,  7]])



```python
x2[:2, :3] # 两行三列
```




    array([[12,  5,  2],
           [ 7,  6,  8]])



```python
x2[:3, ::2] # 所有行，每隔一列
```




    array([[12,  2],
           [ 7,  8],
           [ 1,  7]])



```python
x2[::-1, ::-1] # 逆序
```




    array([[ 7,  7,  6,  1],
           [ 8,  8,  6,  7],
           [ 4,  2,  5, 12]])


#### 获取数组的行和列


```python
print(x2[:, 0]) # x2的第一列
```

    [12  7  1]


```python
print(x2[0,:]) # x2的第一行
```

    [12  5  2  4]

在获取行时也可以：


```python
print(x2[0])
```

    [12  5  2  4]

#### 非副本视图的子数组

数组切片返回的是数组数据的**视图**，而不是数组数据的**副本**。


```python
print(x2)
```

    [[12  5  2  4]
     [ 7  6  8  8]
     [ 1  6  7  7]]

从中抽取一个2*2的子数组：


```python
x2_sub = x2[:2, :2]
print(x2_sub)
```

    [[12  5]
     [ 7  6]]

如果修改这个子数组，原始数组也会被改变：


```python
x2_sub[0, 0] = 99
print(x2_sub)
```

    [[99  5]
     [ 7  6]]


```python
print(x2)
```

    [[99  5  2  4]
     [ 7  6  8  8]
     [ 1  6  7  7]]

这种特性意味着，当我们处理非常大的数据集时，可以直接获取这些数据的片段进行处理。

#### 创建数组的副本

有时候我们也想明确地复制数组里的数据或子数据，这时候可以使用`copy()`方法实现：


```python
x2_sub_copy = x2[:2, :2].copy()
print(x2_sub_copy)
```

    [[99  5]
     [ 7  6]]

如果修改这个数组，原始的数组不会被改变：


```python
x2_sub_copy[0, 0] = 42
print(x2_sub_copy)
```

    [[42  5]
     [ 7  6]]


```python
print(x2)
```

    [[99  5  2  4]
     [ 7  6  8  8]
     [ 1  6  7  7]]

### 数组的变形

使用`reshape()`方法进行变形，例如你希望将数字1~9放入一个3*3的矩阵中：


```python
grid = np.arange(1, 10).reshape((3, 3))
print(grid)
```

    [[1 2 3]
     [4 5 6]
     [7 8 9]]

注意，原始数组的大小必须和变形后数组的大小一致。你也可以通过`reshape()`将一个一维数组变为二维的行或列的矩阵：


```python
x = np.array([1, 2, 3])

# 通过变形获得的行向量
x.reshape((1, 3))
```




    array([[1, 2, 3]])


也可以通过`newaxis`来进行相同的操作：


```python
# 通过newaxis获得的行向量
x[np.newaxis, :]
```




    array([[1, 2, 3]])



```python
# 通过变形获得的列向量
x.reshape((3, 1))
```




    array([[1],
           [2],
           [3]])



```python
# 通过newaxis获得的列向量
x[:, np.newaxis]
```




    array([[1],
           [2],
           [3]])


### 数组的拼接和分裂

以上操作都是针对单一数组的，有时候我们也需要把多个数组成为一个，或将一个数组拆分成多个。

#### 数组的拼接


```python
x = np.array([1, 2, 3])
y = np.array([3, 2, 1])
np.concatenate([x, y])
```




    array([1, 2, 3, 3, 2, 1])


也可以一次性拼接两个以上数组：


```python
z = [99, 99, 99]
np.concatenate([x, y, z])
```




    array([ 1,  2,  3,  3,  2,  1, 99, 99, 99])


`np.concatenate`也可以用于二维数组的拼接：


```python
grid = np.arange(1,7).reshape((2, 3))
grid
```




    array([[1, 2, 3],
           [4, 5, 6]])



```python
# 沿着第一个轴拼接
np.concatenate([grid, grid])
```




    array([[1, 2, 3],
           [4, 5, 6],
           [1, 2, 3],
           [4, 5, 6]])



```python
# 沿着第二个轴拼接（从0开始索引）
np.concatenate([grid, grid], axis=1)
```




    array([[1, 2, 3, 1, 2, 3],
           [4, 5, 6, 4, 5, 6]])


沿着固定维度处理初级使用`np.vstack`（垂直）和`np.hstack`（水平）会更简洁：


```python
x = np.array([1, 2, 3])
grid = np.array([[9, 8, 7],
                 [6, 5, 4]])

# 垂直拼接数据
np.vstack([x, grid])
```




    array([[1, 2, 3],
           [9, 8, 7],
           [6, 5, 4]])



```python
# 水平拼接数据
y = np.array([[99],
              [99]])
np.hstack([grid, y])
```




    array([[ 9,  8,  7, 99],
           [ 6,  5,  4, 99]])


#### 数组的分裂

同样，分裂可以通过`np.split`、`np.hsplit`和`np.vsplit`来实现，参数为数组和分裂点的位置索引：


```python
x = [1, 2, 3, 99, 99, 3, 2, 1]
x1, x2, x3 = np.split(x, [3, 5])
print(x1, x2, x3)
```

    [1 2 3] [99 99] [3 2 1]

`np.hsplit`和`np.vsplit`的用法也类似：


```python
grid = np.arange(16).reshape((4, 4))
grid
```




    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])



```python
upper, lower = np.vsplit(grid, [2])
print(upper)
print(lower)
```

    [[0 1 2 3]
     [4 5 6 7]]
    [[ 8  9 10 11]
     [12 13 14 15]]


```python
left, right = np.hsplit(grid, [2])
print(left)
print(right)
```

```
    [[ 0  1]
     [ 4  5]
     [ 8  9]
     [12 13]]
    [[ 2  3]
     [ 6  7]
     [10 11]
     [14 15]]
```

> 笔记整理自《Python数据科学手册》，本书的英文版以及一些资料已在[Github开源](https://github.com/jakevdp/PythonDataScienceHandbook)。
