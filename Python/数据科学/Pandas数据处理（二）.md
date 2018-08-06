
# Pandas数据处理（二）

## 数值取值与选择

### Series数据选择方法

`Series`对象和一维的NumPy数组和标准的Python字典在许多方面都一样。

#### 将Series看做字典

`Series`对象和字典一样有键值对的映射：


```python
import pandas as pd
data = pd.Series([0.25, 0.5, 0.75, 1.0], index=['a', 'b', 'c', 'd'])
data
```




    a    0.25
    b    0.50
    c    0.75
    d    1.00
    dtype: float64




```python
data['b']
```




    0.5



还可以用Python字典的表达方式来检测键和值：


```python
'a' in data
```




    True




```python
data.keys()
```




    Index(['a', 'b', 'c', 'd'], dtype='object')




```python
list(data.items())
```




    [('a', 0.25), ('b', 0.5), ('c', 0.75), ('d', 1.0)]



和通过增加新的键来扩展字典一样，也可以通过这种方法扩展`Series`对象：


```python
data['e'] = 1.25
data
```




    a    0.25
    b    0.50
    c    0.75
    d    1.00
    e    1.25
    dtype: float64



#### 将Series看做一维数组

`Series`对象也有着类似数组的操作，包括索引、掩码、花哨的索引等操作：


```python
# 根据显式索引切片
data['a':'c']
```




    a    0.25
    b    0.50
    c    0.75
    dtype: float64




```python
# 根据隐式索引切片
data[0:2]
```




    a    0.25
    b    0.50
    dtype: float64




```python
# 掩码
data[(data > 0.3) & (data < 0.8)]
```




    b    0.50
    c    0.75
    dtype: float64




```python
# 花哨的索引
data[['a', 'e']]
```




    a    0.25
    e    1.25
    dtype: float64



#### 索引器：loc、iloc和ix

有时候，整数索引会导致混乱，比如，当你的`Series`对象是显式整数索引时，`data[1]`这样的取值操作会使用显式索引，而`data[1:3]`这样的切片操作却会使用隐式索引:


```python
data = pd.Series(['a', 'b', 'c'], index=[1, 3, 5])
data
```




    1    a
    3    b
    5    c
    dtype: object




```python
# 取值操作是显式索引
data[1]
```




    'a'




```python
# 切片操作是隐式索引
data[1:3]
```




    3    b
    5    c
    dtype: object



为了避免这些混乱，Pandas提供了索引器来取值：

- `loc`属性表示取值和切片都为显式：


```python
data.loc[1]
```




    'a'




```python
data.loc[1:3]
```




    1    a
    3    b
    dtype: object



- `iloc`属性表示取值和切片都为隐式：


```python
data.iloc[1]
```




    'b'




```python
data.iloc[1:3]
```




    3    b
    5    c
    dtype: object



使用`loc`和`iloc`可以使代码更易维护。

### DataFrame数据选择方法

#### 将DataFrame看做字典

可以将`DataFrame`当做由若干`Series`对象构成的字典。用之前的美国五州面积和人口数据来演示：


```python
area = pd.Series({'California': 423967, 'Texas': 695662,
                  'New York': 141297, 'Florida': 170312,
                  'Illinois': 149995})
pop = pd.Series({'California': 38332521, 'Texas': 26448193,
                 'New York': 19651127, 'Florida': 19552860,
                 'Illinois': 12882135})
data = pd.DataFrame({'area':area, 'pop':pop})
data
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
      <th>pop</th>
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



两个`Series`对象构成列，可以通过字典形式来取值：


```python
data['area']
```




    California    423967
    Florida       170312
    Illinois      149995
    New York      141297
    Texas         695662
    Name: area, dtype: int64



也可以使用属性的形式：


```python
data.area
```




    California    423967
    Florida       170312
    Illinois      149995
    New York      141297
    Texas         695662
    Name: area, dtype: int64



和`Series`对象一样，加一列可以这样做：


```python
data['density'] = data['pop'] / data['area']
data
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
      <th>pop</th>
      <th>density</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>California</th>
      <td>423967</td>
      <td>38332521</td>
      <td>90.413926</td>
    </tr>
    <tr>
      <th>Florida</th>
      <td>170312</td>
      <td>19552860</td>
      <td>114.806121</td>
    </tr>
    <tr>
      <th>Illinois</th>
      <td>149995</td>
      <td>12882135</td>
      <td>85.883763</td>
    </tr>
    <tr>
      <th>New York</th>
      <td>141297</td>
      <td>19651127</td>
      <td>139.076746</td>
    </tr>
    <tr>
      <th>Texas</th>
      <td>695662</td>
      <td>26448193</td>
      <td>38.018740</td>
    </tr>
  </tbody>
</table>
</div>



#### 将DataFrame看做二维数组

可以把`DataFrame`看做一个增强版的二维数组，用`values`属性查看数据：


```python
data.values
```




    array([[4.23967000e+05, 3.83325210e+07, 9.04139261e+01],
           [1.70312000e+05, 1.95528600e+07, 1.14806121e+02],
           [1.49995000e+05, 1.28821350e+07, 8.58837628e+01],
           [1.41297000e+05, 1.96511270e+07, 1.39076746e+02],
           [6.95662000e+05, 2.64481930e+07, 3.80187404e+01]])



可以把许多数组操作方式用在`DataFrame`上，比如进行转置：


```python
data.T
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
      <th>California</th>
      <th>Florida</th>
      <th>Illinois</th>
      <th>New York</th>
      <th>Texas</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>area</th>
      <td>4.239670e+05</td>
      <td>1.703120e+05</td>
      <td>1.499950e+05</td>
      <td>1.412970e+05</td>
      <td>6.956620e+05</td>
    </tr>
    <tr>
      <th>pop</th>
      <td>3.833252e+07</td>
      <td>1.955286e+07</td>
      <td>1.288214e+07</td>
      <td>1.965113e+07</td>
      <td>2.644819e+07</td>
    </tr>
    <tr>
      <th>density</th>
      <td>9.041393e+01</td>
      <td>1.148061e+02</td>
      <td>8.588376e+01</td>
      <td>1.390767e+02</td>
      <td>3.801874e+01</td>
    </tr>
  </tbody>
</table>
</div>



通过索引器切片：


```python
data.iloc[:3, :2]
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
      <th>pop</th>
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
  </tbody>
</table>
</div>




```python
data.loc[:'Illinois', :'pop']
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
      <th>pop</th>
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
  </tbody>
</table>
</div>



也可以在`loc`索引器中结合掩码和花哨的索引：


```python
data.loc[data.density > 100, ['pop', 'density']]
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
      <th>pop</th>
      <th>density</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Florida</th>
      <td>19552860</td>
      <td>114.806121</td>
    </tr>
    <tr>
      <th>New York</th>
      <td>19651127</td>
      <td>139.076746</td>
    </tr>
  </tbody>
</table>
</div>



#### 其他取值方法


```python
data['Florida':'Illinois']
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
      <th>pop</th>
      <th>density</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Florida</th>
      <td>170312</td>
      <td>19552860</td>
      <td>114.806121</td>
    </tr>
    <tr>
      <th>Illinois</th>
      <td>149995</td>
      <td>12882135</td>
      <td>85.883763</td>
    </tr>
  </tbody>
</table>
</div>




```python
data[1:3]
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
      <th>pop</th>
      <th>density</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Florida</th>
      <td>170312</td>
      <td>19552860</td>
      <td>114.806121</td>
    </tr>
    <tr>
      <th>Illinois</th>
      <td>149995</td>
      <td>12882135</td>
      <td>85.883763</td>
    </tr>
  </tbody>
</table>
</div>




```python
data[data.density > 100]
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
      <th>pop</th>
      <th>density</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Florida</th>
      <td>170312</td>
      <td>19552860</td>
      <td>114.806121</td>
    </tr>
    <tr>
      <th>New York</th>
      <td>141297</td>
      <td>19651127</td>
      <td>139.076746</td>
    </tr>
  </tbody>
</table>
</div>



## Pandas数值运算方法

Pandas除了继承NumPy的功能之外，对于一元运算（如函数与三角函数）还会在输出结果中保留索引和列标签，对于二元运算（如加法和乘法），Pandas会自动对齐索引。

### 通用函数：保留索引


```python
import pandas as pd
import numpy as np
```


```python
rng = np.random.RandomState(42)
ser = pd.Series(rng.randint(0, 10, 4))
ser
```




    0    6
    1    3
    2    7
    3    4
    dtype: int32




```python
df = pd.DataFrame(rng.randint(0, 10, (3, 4)),
                  columns=['A', 'B', 'C', 'D'])
df
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
      <th>C</th>
      <th>D</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>6</td>
      <td>9</td>
      <td>2</td>
      <td>6</td>
    </tr>
    <tr>
      <th>1</th>
      <td>7</td>
      <td>4</td>
      <td>3</td>
      <td>7</td>
    </tr>
    <tr>
      <th>2</th>
      <td>7</td>
      <td>2</td>
      <td>5</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>



对这两个对象都可以使用NumPy的通用函数，生成的结果是一个保留索引的Pandas对象：


```python
np.exp(ser) # e^x
```




    0     403.428793
    1      20.085537
    2    1096.633158
    3      54.598150
    dtype: float64




```python
np.sin(df * np.pi / 4)
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
      <th>C</th>
      <th>D</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>-1.000000</td>
      <td>7.071068e-01</td>
      <td>1.000000</td>
      <td>-1.000000e+00</td>
    </tr>
    <tr>
      <th>1</th>
      <td>-0.707107</td>
      <td>1.224647e-16</td>
      <td>0.707107</td>
      <td>-7.071068e-01</td>
    </tr>
    <tr>
      <th>2</th>
      <td>-0.707107</td>
      <td>1.000000e+00</td>
      <td>-0.707107</td>
      <td>1.224647e-16</td>
    </tr>
  </tbody>
</table>
</div>



### 通用函数：索引对齐

当进行二元计算时，Pandas会在计算过程中自动对齐两个对象的索引。

#### Series索引对齐

比如要整合两组数据，一个是美国面积最大的三个州的面积数据，另一个是美国人口最多的三个州的数据：


```python
area = pd.Series({'Alaska': 1723337, 'Texas': 695662,
                  'California': 423967}, name='area')
population = pd.Series({'California': 38332521, 'Texas': 26448193,
                        'New York': 19651127}, name='population')
```

看看人口除以面积是什么结果：


```python
population / area
```




    Alaska              NaN
    California    90.413926
    New York            NaN
    Texas         38.018740
    dtype: float64



可见结果是两个数组索引的并集，缺失的位置用`NaN`填充。也可以自定义缺失的数据：


```python
population.div(area, fill_value=1)
```




    Alaska        5.802696e-07
    California    9.041393e+01
    New York      1.965113e+07
    Texas         3.801874e+01
    dtype: float64



#### DataFrame索引对齐


```python
A = pd.DataFrame(rng.randint(0, 20, (2, 2)),
                columns=list('AB'))
A
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
      <td>6</td>
      <td>8</td>
    </tr>
    <tr>
      <th>1</th>
      <td>6</td>
      <td>17</td>
    </tr>
  </tbody>
</table>
</div>




```python
B = pd.DataFrame(rng.randint(0, 10, (3, 3)),
                 columns=list('BAC'))
B
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
      <th>B</th>
      <th>A</th>
      <th>C</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>6</td>
      <td>7</td>
      <td>2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0</td>
      <td>3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>7</td>
      <td>3</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
A + B
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
      <th>C</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>14.0</td>
      <td>11.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>14.0</td>
      <td>26.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>



你会发现，两个对象的索引是不同顺序的，但Pandas会自动按顺序排列进行计算。

### 通用函数：DataFrame与Series的运算

`DataFrame`和`Series`的运算与NumPy中二维数组和一维数组的运算规则是一样的。

例如，根据NumPy的广播规则（参见NumPy入门（二）），二维数组减去一行数据会按行计算。


```python
A = rng.randint(10, size=(3, 4))
A
```




    array([[5, 5, 9, 3],
           [5, 1, 9, 1],
           [9, 3, 7, 6]])




```python
A - A[0] # 减去第一行
```




    array([[ 0,  0,  0,  0],
           [ 0, -4,  0, -2],
           [ 4, -2, -2,  3]])



在Pandas里默认也是按行运算的：


```python
df = pd.DataFrame(A, columns=list('QRST'))
df - df.iloc[0]
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
      <th>Q</th>
      <th>R</th>
      <th>S</th>
      <th>T</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0</td>
      <td>-4</td>
      <td>0</td>
      <td>-2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>4</td>
      <td>-2</td>
      <td>-2</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>



如果想按列运算，可以通过`axis`参数设置：


```python
df.subtract(df['R'], axis=0)
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
      <th>Q</th>
      <th>R</th>
      <th>S</th>
      <th>T</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>-2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4</td>
      <td>0</td>
      <td>8</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>6</td>
      <td>0</td>
      <td>4</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>



``DataFrame``/``Series`` 的运算，结果的索引都会自动对齐：


```python
halfrow = df.iloc[0, ::2]
halfrow
```




    Q    5
    S    9
    Name: 0, dtype: int32




```python
df - halfrow
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
      <th>Q</th>
      <th>R</th>
      <th>S</th>
      <th>T</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.0</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>4.0</td>
      <td>NaN</td>
      <td>-2.0</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
</div>


> 笔记整理自《Python数据科学手册》，本书的英文版以及一些资料已在[Github开源](https://github.com/jakevdp/PythonDataScienceHandbook)（https://github.com/jakevdp/PythonDataScienceHandbook ）。