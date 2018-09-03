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

> 笔记整理自《Python数据科学手册》，本书的英文版以及一些资料已在[Github开源](https://github.com/jakevdp/PythonDataScienceHandbook)（https://github.com/jakevdp/PythonDataScienceHandbook ）。
