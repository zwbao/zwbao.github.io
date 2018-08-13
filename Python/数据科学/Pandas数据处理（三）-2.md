# Pandas数据处理（三）

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
