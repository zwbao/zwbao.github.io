
# Pandas数据处理（三）

## 处理缺失值

### Pandas的缺失值

1. `None`：Python对象类型的缺失值


```python
import numpy as np
import pandas as pd
```


```python
vals1 = np.array([1, None, 3, 4])
vals1
```




    array([1, None, 3, 4], dtype=object)



`None`是一个Python对象，这里的`dtype=object`表示NumPy认为该数组由Python对象生成，其类型为`object`。使用Python对象构成的数组意味着不能对该数组进行`sum`或`min`等操作，会出现报错：


```python
vals1.sum()
```


    ---------------------------------------------------------------------------

    TypeError                                 Traceback (most recent call last)

    <ipython-input-3-30a3fc8c6726> in <module>()
    ----> 1 vals1.sum()
    

    c:\programdata\anaconda3\lib\site-packages\numpy\core\_methods.py in _sum(a, axis, dtype, out, keepdims)
         30 
         31 def _sum(a, axis=None, dtype=None, out=None, keepdims=False):
    ---> 32     return umr_sum(a, axis, dtype, out, keepdims)
         33 
         34 def _prod(a, axis=None, dtype=None, out=None, keepdims=False):
    

    TypeError: unsupported operand type(s) for +: 'int' and 'NoneType'



```python
vals1.min()
```


    ---------------------------------------------------------------------------

    TypeError                                 Traceback (most recent call last)

    <ipython-input-5-662df191e581> in <module>()
    ----> 1 vals1.min()
    

    c:\programdata\anaconda3\lib\site-packages\numpy\core\_methods.py in _amin(a, axis, out, keepdims)
         27 
         28 def _amin(a, axis=None, out=None, keepdims=False):
    ---> 29     return umr_minimum(a, axis, None, out, keepdims)
         30 
         31 def _sum(a, axis=None, dtype=None, out=None, keepdims=False):
    

    TypeError: '<=' not supported between instances of 'int' and 'NoneType'


2. `NaN`：数值类型的缺失值

`NaN`全称为`Not a Number`（不是一个数字），`NaN`是一种特殊的浮点数，而不是整数、字符串以及其他类型。


```python
vals2 = np.array([1, np.nan, 3, 4])
vals2.dtype
```




    dtype('float64')



无论对`NaN`进行何种操作都会返回`NaN`：


```python
1 + np.nan
```




    nan




```python
0 * np.nan
```




    nan




```python
vals2.sum(), vals2.min(), vals2.max()
```

    c:\programdata\anaconda3\lib\site-packages\numpy\core\_methods.py:29: RuntimeWarning: invalid value encountered in reduce
      return umr_minimum(a, axis, None, out, keepdims)
    c:\programdata\anaconda3\lib\site-packages\numpy\core\_methods.py:26: RuntimeWarning: invalid value encountered in reduce
      return umr_maximum(a, axis, None, out, keepdims)
    




    (nan, nan, nan)



可以用一些特殊的函数忽略缺失值的影响：


```python
np.nansum(vals2), np.nanmin(vals2), np.nanmax(vals2)
```




    (8.0, 1.0, 4.0)



3. Pandas中`NaN`与`None`的差异

Pandas有时可以将他们看做是等价的，可以互相替换：


```python
pd.Series([1, np.nan, 2, None])
```




    0    1.0
    1    NaN
    2    2.0
    3    NaN
    dtype: float64



Pandas也会将没有值的数据自动转换为`NA`：


```python
x = pd.Series(range(2), dtype=int)
x
```




    0    0
    1    1
    dtype: int32




```python
x[0] = None
x
```




    0    NaN
    1    1.0
    dtype: float64



### 处理缺失值

`isnull()`：创建一个布尔类型的掩码标签缺失值。

`notnull()`：与`isnull()`操作相反。

`dropna()`：返回一个剔除缺失值的数据。

`fillna()`：范湖一个填充了缺失值的数据。

#### 发现缺失值

Pandas可以用`isnull()`和`notnull()`来发现缺失值：


```python
data = pd.Series([1, np.nan, 'hello', None])
```


```python
data.isnull()
```




    0    False
    1     True
    2    False
    3     True
    dtype: bool



布尔类型可以直接作为索引使用：


```python
data[data.notnull()]
```




    0        1
    2    hello
    dtype: object



#### 剔除缺失值

`dropna()`和`fillna()`用于剔除缺失值：


```python
data.dropna()
```




    0        1
    2    hello
    dtype: object



而在`DataFrame`上使用时需要设置一些参数，因为没法在`DataFrame`中剔除一个值，要么是剔除有缺失值的行，要么是列。


```python
df = pd.DataFrame([[1,      np.nan, 2],
                   [2,      3,      5],
                   [np.nan, 4,      6]])
df
```

![df](http://oo3g995ih.bkt.clouddn.com/blog/180813/Kb8C8D8JH1.png?imageslim)

```python
df.dropna()
```

![df.dropna()](http://oo3g995ih.bkt.clouddn.com/blog/180813/ciH4DfcGe4.png?imageslim)

可以按不同的坐标轴剔除缺失值：


```python
df.dropna(axis='columns') # 或 axis=1
```

![df.dropna(axis='columns')](http://oo3g995ih.bkt.clouddn.com/blog/180813/J3JmcFHFGm.png?imageslim)

默认设置是`how='any'`，只要有缺失值就会剔除整行或整列，还可以设置`how='all'`，只删除全部是缺失值的行或列：


```python
df[3] = np.nan
df
```

![df](http://oo3g995ih.bkt.clouddn.com/blog/180813/0KI6lBC79f.png?imageslim)


```python
df.dropna(axis='columns', how='all')
```

![df](http://oo3g995ih.bkt.clouddn.com/blog/180813/GDiGGm6bdK.png?imageslim)

还可通过`thresh`参数设置行或列中非缺失值的最小数量：


```python
df.dropna(axis='rows', thresh=3)
```

![df](http://oo3g995ih.bkt.clouddn.com/blog/180813/9GH2lJd3EH.png?imageslim)

3. 填充缺失值

虽然可以用`isnull()`来建立掩码填充缺失值，但也可以直接用`fillna()`方法，它将返回填充了缺失值后的数组副本。


```python
data = pd.Series([1, np.nan, 2, None, 3], index=list('abcde'))
data
```




    a    1.0
    b    NaN
    c    2.0
    d    NaN
    e    3.0
    dtype: float64




```python
data.fillna(0) # 用0来填充缺失值
```




    a    1.0
    b    0.0
    c    2.0
    d    0.0
    e    3.0
    dtype: float64




```python
data.fillna(method='ffill') # 用缺失值前面的有效值来填充（forward-fill）
```




    a    1.0
    b    1.0
    c    2.0
    d    2.0
    e    3.0
    dtype: float64




```python
data.fillna(method='bfill') # 用缺失值后面的有效值来填充（back-fill）
```




    a    1.0
    b    2.0
    c    2.0
    d    3.0
    e    3.0
    dtype: float64



`DataFrame`与`Series`类似，只是需要添加坐标轴参数`axis`：


```python
df
```

![df](http://oo3g995ih.bkt.clouddn.com/blog/180813/63jI2KLHG4.png?imageslim)

```python
df.fillna(method='ffill', axis=1)
```

![df.fillna](http://oo3g995ih.bkt.clouddn.com/blog/180813/lCaGB6K9KJ.png?imageslim)

## 层级索引

用更简单的方法分析美国各州在两个不同年份的数据。

### Pandas多级索引


```python
import pandas as pd
import numpy as np
```


```python
index = [('California', 2000), ('California', 2010),
         ('New York', 2000), ('New York', 2010),
         ('Texas', 2000), ('Texas', 2010)]
```


```python
index = pd.MultiIndex.from_tuples(index)
index
```




    MultiIndex(levels=[['California', 'New York', 'Texas'], [2000, 2010]],
               labels=[[0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]])



使用`MultiIndex`从元组创建多级索引，`MultiIndex`中的`levels`属性表示索引的等级。


```python
populations = [33871648, 37253956,
               18976457, 19378102,
               20851820, 25145561]
pop = pd.Series(populations, index=index)
pop
```




    California  2000    33871648
                2010    37253956
    New York    2000    18976457
                2010    19378102
    Texas       2000    20851820
                2010    25145561
    dtype: int64



可以使用切片的方法来获取2010年的全部数据：


```python
pop[:, 2010]
```




    California    37253956
    New York      19378102
    Texas         25145561
    dtype: int64



### 高维数据的多级索引

`unstack()`方法可以快速将一个多级索引的`Series`转化为普通索引的`DataFrame`：


```python
pop_df = pop.unstack()
pop_df
```

![pop_df](http://oo3g995ih.bkt.clouddn.com/blog/180813/LAeImCkIag.png?imageslim)

`stack()`方法实现相反的效果：


```python
pop_df.stack()
```




    California  2000    33871648
                2010    37253956
    New York    2000    18976457
                2010    19378102
    Texas       2000    20851820
                2010    25145561
    dtype: int64



有了多级索引，我们就可以用`Series`或`DataFrame`表示三维或更高维度的数据了。例如要添加一列显示每一年各州的人口统计指标（18岁以下的人口），使用`MultiIndex`就十分简单：


```python
pop_df = pd.DataFrame({'total': pop,
                       'under18': [9267089, 9284094,
                                   4687374, 4318033,
                                   5906301, 6879014]})
pop_df
```

![pop_df](http://oo3g995ih.bkt.clouddn.com/blog/180813/F9A9AKB5aa.png?imageslim)

计算18岁以下人口所占的比例：


```python
f_u18 = pop_df['under18'] / pop_df['total']
f_u18.unstack()
```

![f_u18.unstack](http://oo3g995ih.bkt.clouddn.com/blog/180813/9DJFgI6dJ4.png?imageslim)

### 多级索引的创建方法

添加多级索引最直接的方法就是将`index`参数设置为至少二维数组：


```python
df = pd.DataFrame(np.random.rand(4, 2),
                  index=[['a', 'a', 'b', 'b'], [1, 2, 1, 2]],
                  columns=['data1', 'data2'])
df
```

![df](http://oo3g995ih.bkt.clouddn.com/blog/180813/d57581dflJ.png?imageslim)

如果将元组作为字典的键传递给Pandas，Pandas会默认转换为`MultiIndex`：


```python
data = {('California', 2000): 33871648,
        ('California', 2010): 37253956,
        ('Texas', 2000): 20851820,
        ('Texas', 2010): 25145561,
        ('New York', 2000): 18976457,
        ('New York', 2010): 19378102}
pd.Series(data)
```




    California  2000    33871648
                2010    37253956
    New York    2000    18976457
                2010    19378102
    Texas       2000    20851820
                2010    25145561
    dtype: int64



#### 显式地创建多级索引

- 通过包含不同等级的若干简单数组的列表来构建多级索引：


```python
pd.MultiIndex.from_arrays([['a', 'a', 'b', 'b'], [1, 2, 1, 2]])
```




    MultiIndex(levels=[['a', 'b'], [1, 2]],
               labels=[[0, 0, 1, 1], [0, 1, 0, 1]])



- 通过包含多个索引值的元组的列表构建多级索引：


```python
pd.MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1), ('b', 2)])
```




    MultiIndex(levels=[['a', 'b'], [1, 2]],
               labels=[[0, 0, 1, 1], [0, 1, 0, 1]])



- 用两个笛卡尔积构建多级索引：


```python
pd.MultiIndex.from_product([['a', 'b'], [1, 2]])
```




    MultiIndex(levels=[['a', 'b'], [1, 2]],
               labels=[[0, 0, 1, 1], [0, 1, 0, 1]])



- 直接提供`levels`和`labels`创建多级索引：


```python
pd.MultiIndex(levels=[['a', 'b'], [1, 2]],
              labels=[[0, 0, 1, 1], [0, 1, 0, 1]])
```




    MultiIndex(levels=[['a', 'b'], [1, 2]],
               labels=[[0, 0, 1, 1], [0, 1, 0, 1]])



在创建`Series`或`DataFrame`时可以将这些对象作为`index`参数，或通过`reindex`方法来对已有的`Series`或`DataFrame`重新索引。

#### 多级索引的等级名称

使用`names`参数来设置等级名称：


```python
pop.index.names = ['state', 'year']
pop
```




    state       year
    California  2000    33871648
                2010    37253956
    New York    2000    18976457
                2010    19378102
    Texas       2000    20851820
                2010    25145561
    dtype: int64



#### 多级列索引


```python
# 多级行列索引
index = pd.MultiIndex.from_product([[2013, 2014], [1, 2]],
                                   names=['year', 'visit'])
columns = pd.MultiIndex.from_product([['Bob', 'Guido', 'Sue'], ['HR', 'Temp']],
                                     names=['subject', 'type'])

# 模拟数据
data = np.round(np.random.randn(4, 6), 1)
data[:, ::2] *= 10
data += 37

# 创建DataFrame
health_data = pd.DataFrame(data, index=index, columns=columns)
health_data
```

![health_data](http://oo3g995ih.bkt.clouddn.com/blog/180813/F208gg0jLm.png?imageslim)

比如要查看一个人的全部检查信息：


```python
health_data['Guido']
```

![Guido](http://oo3g995ih.bkt.clouddn.com/blog/180813/hla91flbEm.png?imageslim)

对于多个维度的数据，使用多级行列索引进行查询会非常方便。

### 多级索引的取值与切片

#### `Series`多级索引


```python
pop
```




    state       year
    California  2000    33871648
                2010    37253956
    New York    2000    18976457
                2010    19378102
    Texas       2000    20851820
                2010    25145561
    dtype: int64



可以通过多个级别的索引值来获得单个元素：


```python
pop['California', 2000]
```




    33871648



也可以只取索引的某一层级，获得的是一个新的`Series`：


```python
pop['California']
```




    year
    2000    33871648
    2010    37253956
    dtype: int64




```python
pop[:, 2010]
```




    state
    California    37253956
    New York      19378102
    Texas         25145561
    dtype: int64



还可以进行局部切片：


```python
pop.loc['California':'New York']
```




    state       year
    California  2000    33871648
                2010    37253956
    New York    2000    18976457
                2010    19378102
    dtype: int64



通过布尔掩码选择数据：


```python
pop[pop > 22000000]
```




    state       year
    California  2000    33871648
                2010    37253956
    Texas       2010    25145561
    dtype: int64



用花哨的索引选择数据：


```python
pop[['California', 'Texas']]
```




    state       year
    California  2000    33871648
                2010    37253956
    Texas       2000    20851820
                2010    25145561
    dtype: int64



#### `DataFrame`多级索引

由于`DataFrame`的基本索引是列索引，所以所有的操作都在列的水平上进行：


```python
health_data
```

![health_data](http://oo3g995ih.bkt.clouddn.com/blog/180813/DKHlmIKK3e.png?imageslim)

```python
health_data['Guido', 'HR']
```




    year  visit
    2013  1        33.0
          2        44.0
    2014  1        35.0
          2        52.0
    Name: (Guido, HR), dtype: float64



之前介绍的`loc`、`iloc`、`ix`索引器都可使用：


```python
health_data.iloc[:2, :2]
```

![health_data.iloc](http://oo3g995ih.bkt.clouddn.com/blog/180813/9a1K9E928B.png?imageslim)

```python
health_data[['Guido','Sue']]
```

![Guido&Sue](http://oo3g995ih.bkt.clouddn.com/blog/180813/cLeA70KF97.png?imageslim)

但这种索引方法不是很方便，在元组中使用切片还会导致语法错误：


```python
health_data.loc[(:, 1), (:, 'HR')]
```


      File "<ipython-input-125-fb34fa30ac09>", line 1
        health_data.loc[(:, 1), (:, 'HR')]
                         ^
    SyntaxError: invalid syntax
    


更好的办法就是使用`IndexSlice`对象：


```python
idx = pd.IndexSlice
health_data.loc[idx[:, 1], idx[:, 'HR']]
```

![idx](http://oo3g995ih.bkt.clouddn.com/blog/180813/aDdgiH31Aj.png?imageslim)

```python
# 获取Bob的所有数据
health_data.loc[idx[:], idx['Bob',:]]
```

![Bob](http://oo3g995ih.bkt.clouddn.com/blog/180813/K91cgililG.png?imageslim)


```python
# 获取2013年Bob的数据
health_data.loc[idx[2013, :], idx['Bob',:]]
```

![2013-Bob](http://oo3g995ih.bkt.clouddn.com/blog/180813/9HIG7Gb2hl.png?imageslim)

### 多级索引行列转换

#### 有序的索引和无序的索引

如果多级索引不是有序的索引，那么大多数切片操作都会失败：


```python
index = pd.MultiIndex.from_product([['a', 'c', 'b'], [1, 2]])
data = pd.Series(np.random.rand(6), index=index)
data.index.names = ['char', 'int']
data
```




    char  int
    a     1      0.887551
          2      0.753588
    c     1      0.130569
          2      0.514826
    b     1      0.217901
          2      0.661305
    dtype: float64



如果想使用局部切片，就会出现错误：


```python
try:
    data['a':'b']
except KeyError as e:
    print(type(e))
    print(e)
```

    <class 'pandas.errors.UnsortedIndexError'>
    'Key length (1) was greater than MultiIndex lexsort depth (0)'
    

因为局部切片和其他类似的操作都要求索引是有序的，所以Pandas提供了许多方法来进行排序，如`sort_index()`和`sortlevel()`：


```python
data = data.sort_index()
data
```




    char  int
    a     1      0.887551
          2      0.753588
    b     1      0.217901
          2      0.661305
    c     1      0.130569
          2      0.514826
    dtype: float64




```python
data['a':'b']
```




    char  int
    a     1      0.887551
          2      0.753588
    b     1      0.217901
          2      0.661305
    dtype: float64



#### 索引`stack`和`unstack`

可以通过`level`参数来设置转换的索引层级：


```python
pop
```




    state       year
    California  2000    33871648
                2010    37253956
    New York    2000    18976457
                2010    19378102
    Texas       2000    20851820
                2010    25145561
    dtype: int64




```python
pop.unstack(level=0)
```

![level=0](http://oo3g995ih.bkt.clouddn.com/blog/180813/gJgA20kEh8.png?imageslim)

```python
pop.unstack(level=1)
```

![level=1](http://oo3g995ih.bkt.clouddn.com/blog/180813/EeAbdcJCka.png?imageslim)

`unstack`是`stack`的逆操作，同时使用这两种方法，数据不会发生任何变化：


```python
pop.unstack().stack()
```




    state       year
    California  2000    33871648
                2010    37253956
    New York    2000    18976457
                2010    19378102
    Texas       2000    20851820
                2010    25145561
    dtype: int64



#### 索引的设置与重置

对上面的`Series`使用`reset_index`方法会生成一个包含之前索引的`DataFrame`:


```python
pop_flat = pop.reset_index()
pop_flat
```

![pop_flat](http://oo3g995ih.bkt.clouddn.com/blog/180813/I2iFc7ec8J.png?imageslim)

也可以使用数据的`name`属性为列设置名称：


```python
pop_flat = pop.reset_index(name='population')
pop_flat
```

![pop_flat](http://oo3g995ih.bkt.clouddn.com/blog/180813/2DlEkA78ce.png?imageslim)

可以通过`set_index`方法来自动生成多级索引：


```python
pop_flat.set_index(['state', 'year'])
```

![set_index](http://oo3g995ih.bkt.clouddn.com/blog/180813/8EG7a17fEb.png?imageslim)

### 多级索引的数据累计方法

对于多级索引，可以通过`level`来实现对数据子集的统计操作：(以体检数据为例)


```python
health_data
```

![health_data](http://oo3g995ih.bkt.clouddn.com/blog/180813/fd131935DJ.png?imageslim)

获得每一项指标的每年的平均值，将`level`设置为`year`：


```python
data_mean = health_data.mean(level='year')
# 相当于 data_mean = health_data.mean(axis=0, level='year')
data_mean
```

![data_mean](http://oo3g995ih.bkt.clouddn.com/blog/180813/0gGci1F1eJ.png?imageslim)

如果再设置`axis`就可以对列索引进行操作：


```python
data_mean.mean(axis=1, level='type')
```

![data_mean.mean](http://oo3g995ih.bkt.clouddn.com/blog/180813/iKl1F66hlI.png?imageslim)

> 笔记整理自《Python数据科学手册》，本书的英文版以及一些资料已在[Github开源](https://github.com/jakevdp/PythonDataScienceHandbook)（https://github.com/jakevdp/PythonDataScienceHandbook ）。
