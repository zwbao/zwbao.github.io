## 安装并使用Pandas

`Pandas`是一个基于`NumPy`的库，其中纳入了大量库和一些标准的数据模型，提供了高效地操作大型数据集所需的工具。

在安装`Pandas`之前需要安装`NumPy`，详细的安装方法请参考`Pandas`官方文档（http://pandas.pydata.org/ ）。如果你已经安装了Anaconda，`Pandas`已包含在其中。

和导入`NumPy` 使用别名 `np` 一样，`Pandas` 也通常使用别名 `pd` 导入：


```python
import pandas as pd
```

## Pandas的对象

首先导入`NumPy`和`Pandas`：


```python
import numpy as np
import pandas as pd
```

### Pandas的Series对象

`Pandas`的`Series`对象是一个带索引的一维数组：


```python
data = pd.Series([0.25, 0.5, 0.75, 1.0])
data
```




    0    0.25
    1    0.50
    2    0.75
    3    1.00
    dtype: float64



我们可以通过`values`属性和`index`属性来分别获取数据和索引：


```python
data.values
```




    array([0.25, 0.5 , 0.75, 1.  ])




```python
data.index
```




    RangeIndex(start=0, stop=4, step=1)



和`NumPy`一样数据也可以通过中括号索引标签获取：


```python
data[1]
```




    0.5




```python
data[1:3]
```




    1    0.50
    2    0.75
    dtype: float64



#### Series是通用的NumPy数组

NumPy数组和Series对象之间的本质差异在于索引：NumPy数组通过隐式索引获取数值，而Series对象除了隐式索引，还可以用显式索引来获取数值。显式索引可以让索引不仅仅是整数，可以是任意想要的类型：


```python
data = pd.Series([0.25, 0.5, 0.75, 1.0],
                 index=['a', 'b', 'c', 'd'])
data
```




    a    0.25
    b    0.50
    c    0.75
    d    1.00
    dtype: float64



用显式索引获取数值：


```python
data['b']
```




    0.5



使用不连续的索引：


```python
data = pd.Series([0.25, 0.5, 0.75, 1.0],
                 index=[2, 5, 3, 7])
data
```




    2    0.25
    5    0.50
    3    0.75
    7    1.00
    dtype: float64




```python
data[5]
```




    0.5



#### Series是特殊的字典

直接用Python中的字典创建一个Series对象：


```python
population_dict = {'California': 38332521,
                   'Texas': 26448193,
                   'New York': 19651127,
                   'Florida': 19552860,
                   'Illinois': 12882135}
population = pd.Series(population_dict)
population
```




    California    38332521
    Florida       19552860
    Illinois      12882135
    New York      19651127
    Texas         26448193
    dtype: int64



典型的自定数值获取方式依然有效：


```python
population['California']
```




    38332521



和字典不同的是，Series对象还支持数组形式的操作，比如切片：


```python
population['California':'Illinois']
```




    California    38332521
    Florida       19552860
    Illinois      12882135
    dtype: int64



#### 创建Series对象

`Series`对象一般是这样的形式：

```python
>>> pd.Series(data, index=index)
```

`index`是可选参数，`data`支持多种数据类型。例如，`data`可以是列表或NumPy数组，`index`默认是整数序列：


```python
pd.Series([2, 4, 6])
```




    0    2
    1    4
    2    6
    dtype: int64



`data`也可以是一个标量，Series对象会自动重复填充到每个索引上：


```python
pd.Series(5, index=[100, 200, 300])
```




    100    5
    200    5
    300    5
    dtype: int64



`data`还可以是字典，`index`默认是排序的键：


```python
pd.Series({2:'a', 1:'b', 3:'c'})
```




    1    b
    2    a
    3    c
    dtype: object



每一种形式都可以通过显式索引筛选需要的结果：


```python
pd.Series({2:'a', 1:'b', 3:'c'}, index=[3, 2])
```




    3    c
    2    a
    dtype: object



### Pandas的DataFrame对象

#### DataFrame是通用的NumPy数组

可以把`Series`看做有更为灵活的索引的一维数组，那么`DataFrame`就是既有灵活行索引也有灵活列索引的二维数组，用上面美国五个州面积的数据来演示：


```python
area_dict = {'California': 423967, 'Texas': 695662, 'New York': 141297,
             'Florida': 170312, 'Illinois': 149995}
area = pd.Series(area_dict)
area
```




    California    423967
    Florida       170312
    Illinois      149995
    New York      141297
    Texas         695662
    dtype: int64



现在我们可以用一个字典创建包含这些信息的二维数组：


```python
states = pd.DataFrame({'population': population,
                       'area': area})
states
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>area</th>
      <th>population</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>California</th>
      <td>423967</td>
      <td>38332521</td>
    </tr>
    <tr>
      <th>Florida</th>
      <td>170312</td>
      <td>19552860</td>
    </tr>
    <tr>
      <th>Illinois</th>
      <td>149995</td>
      <td>12882135</td>
    </tr>
    <tr>
      <th>New York</th>
      <td>141297</td>
      <td>19651127</td>
    </tr>
    <tr>
      <th>Texas</th>
      <td>695662</td>
      <td>26448193</td>
    </tr>
  </tbody>
</table>
</div>



和`Series`对象一样，`DataFrame`也有`index`属性可以获得索引标签:


```python
states.index
```




    Index(['California', 'Florida', 'Illinois', 'New York', 'Texas'], dtype='object')



`DataFrame`还有`columns`属性，来存放列索引：


```python
states.columns
```




    Index(['area', 'population'], dtype='object')



#### DataFrame是特殊的字典

需要注意的是，DataFrame是一列映射一个`Series`数据。例如，通过`'area'`列属性可以返回包含面积数据的`Series`对象：


```python
states['area']
```




    California    423967
    Florida       170312
    Illinois      149995
    New York      141297
    Texas         695662
    Name: area, dtype: int64



#### 创建DataFrame对象

- 通过单个`Series`对象创建。`DataFrame`是一组`Series`对象的集合，单个`Series`对象创建一个单列的`DataFrame`：


```python
pd.DataFrame(population, columns=['population'])
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>population</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>California</th>
      <td>38332521</td>
    </tr>
    <tr>
      <th>Florida</th>
      <td>19552860</td>
    </tr>
    <tr>
      <th>Illinois</th>
      <td>12882135</td>
    </tr>
    <tr>
      <th>New York</th>
      <td>19651127</td>
    </tr>
    <tr>
      <th>Texas</th>
      <td>26448193</td>
    </tr>
  </tbody>
</table>
</div>



- 通过字典列表创建：


```python
data = [{'a': i, 'b': 2 * i}
        for i in range(3)]
pd.DataFrame(data)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>a</th>
      <th>b</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>



要是字典中有些键不存在，Pandas会用`NaN`来填充：


```python
pd.DataFrame([{'a': 1, 'b': 2}, {'b': 3, 'c': 4}])
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>a</th>
      <th>b</th>
      <th>c</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1.0</td>
      <td>2</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NaN</td>
      <td>3</td>
      <td>4.0</td>
    </tr>
  </tbody>
</table>
</div>



- 通过由`Series`对象组成的字典创建：


```python
pd.DataFrame({'population': population,
              'area': area})
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>area</th>
      <th>population</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>California</th>
      <td>423967</td>
      <td>38332521</td>
    </tr>
    <tr>
      <th>Florida</th>
      <td>170312</td>
      <td>19552860</td>
    </tr>
    <tr>
      <th>Illinois</th>
      <td>149995</td>
      <td>12882135</td>
    </tr>
    <tr>
      <th>New York</th>
      <td>141297</td>
      <td>19651127</td>
    </tr>
    <tr>
      <th>Texas</th>
      <td>695662</td>
      <td>26448193</td>
    </tr>
  </tbody>
</table>
</div>



- 通过NumPy二维数组创建：


```python
pd.DataFrame(np.random.rand(3, 2),
             columns=['foo', 'bar'],
             index=['a', 'b', 'c'])
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>foo</th>
      <th>bar</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>a</th>
      <td>0.357701</td>
      <td>0.280100</td>
    </tr>
    <tr>
      <th>b</th>
      <td>0.817338</td>
      <td>0.278349</td>
    </tr>
    <tr>
      <th>c</th>
      <td>0.062861</td>
      <td>0.917770</td>
    </tr>
  </tbody>
</table>
</div>



如果不指定行列索引值，那么行列默认都是整数索引值：


```python
pd.DataFrame(np.random.rand(3, 2))
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.150498</td>
      <td>0.096370</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.299034</td>
      <td>0.158888</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0.384735</td>
      <td>0.402477</td>
    </tr>
  </tbody>
</table>
</div>



- 通过NumPy结构化数组创建：


```python
A = np.zeros(3, dtype=[('A', 'i8'), ('B', 'f8')])
A
```




    array([(0, 0.), (0, 0.), (0, 0.)], dtype=[('A', '<i8'), ('B', '<f8')])




```python
pd.DataFrame(A)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>A</th>
      <th>B</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>



### Pandas的Index对象

`Series`对象和`DataFrame`对象都使用了灵活的显式索引。Pandas的Index对象可以看做是一个不可变数组或有序集合。创建一个简单的整数列表的`Index`对象：


```python
ind = pd.Index([2, 3, 5, 7, 11])
ind
```




    Int64Index([2, 3, 5, 7, 11], dtype='int64')



#### 将Index看做不可变数组


```python
ind[1] # 标准的取值操作
```




    3




```python
ind[::2] # 切片
```




    Int64Index([2, 5, 11], dtype='int64')



`Index`对象还有一些和NumPy数组相似的属性：


```python
print(ind.size, ind.shape, ind.ndim, ind.dtype)
```

    5 (5,) 1 int64
    

`Index`对象和NumPy数组的不同之处在于，`Index`对象的索引是不可变的：


```python
ind[1] = 0
```


    ---------------------------------------------------------------------------

    TypeError                                 Traceback (most recent call last)

    <ipython-input-44-906a9fa1424c> in <module>()
    ----> 1 ind[1] = 0
    

    C:\ProgramData\Anaconda3\lib\site-packages\pandas\core\indexes\base.py in __setitem__(self, key, value)
       1722 
       1723     def __setitem__(self, key, value):
    -> 1724         raise TypeError("Index does not support mutable operations")
       1725 
       1726     def __getitem__(self, key):
    

    TypeError: Index does not support mutable operations


这种特性会在整合不同数据时更加安全。

#### 将Index看做有序集合

`Index`对象可以使用数据结构中的许多习惯用法，比如并集、交集、差集等等：


```python
indA = pd.Index([1, 3, 5, 7, 9])
indB = pd.Index([2, 3, 5, 7, 11])
```


```python
indA & indB # 交集
```




    Int64Index([3, 5, 7], dtype='int64')




```python
indA | indB # 并集
```




    Int64Index([1, 2, 3, 5, 7, 9, 11], dtype='int64')




```python
indA ^ indB # 异或
```




    Int64Index([1, 2, 9, 11], dtype='int64')


> 笔记整理自《Python数据科学手册》，本书的英文版以及一些资料已在[Github开源](https://github.com/jakevdp/PythonDataScienceHandbook)（https://github.com/jakevdp/PythonDataScienceHandbook ）。
