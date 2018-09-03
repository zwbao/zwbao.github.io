## 合并数据集：Concat 与 Append 操作

首先导入 Pandas 和 NumPy：


```python
import pandas as pd
import numpy as np
```

为了快速地创建DataFrame，可创建函数`make_df`：


```python
def make_df(cols, ind):
    data = {c: [str(c) + str(i) for i in ind]
           for c in cols}
    return pd.DataFrame(data, ind)

# 示例
make_df('ABC', range(3))
```

![make_df 示例](http://oo3g995ih.bkt.clouddn.com/blog/180819/DC5lk3e4e5.png?imageslim)

为了并排地显示数据，定义一个类：


```python
class display(object):
    """Display HTML representation of multiple objects"""
    template = """<div style="float: left; padding: 10px;">
    <p style='font-family:"Courier New", Courier, monospace'>{0}</p>{1}
    </div>"""
    def __init__(self, *args):
        self.args = args
        
    def _repr_html_(self):
        return '\n'.join(self.template.format(a, eval(a)._repr_html_())
                         for a in self.args)
    
    def __repr__(self):
        return '\n\n'.join(a + '\n' + repr(eval(a))
                           for a in self.args)
```

### 回顾：NumPy数组的合并

合并`Series`和`DataFrame`与合并NumPy数组类似，可使用NumPy入门（一）中提及的`np.concatenate`：


```python
x = np.array([1, 2, 3])
y = np.array([3, 2, 1])
np.concatenate([x, y])
```




    array([1, 2, 3, 3, 2, 1])



第一个参数是需要合并的数组列表或元组，第二个可选参数是`axis`可设置合并的坐标轴方向：


```python
x = [[1, 2],
     [3, 4]]
np.concatenate([x, x], axis = 1)
```




    array([[1, 2, 1, 2],
           [3, 4, 3, 4]])



### 通过 pd.concat 实现简易合并

Pandas中的`pd.concat`与`np.concatenate`类似，但可选参数更多，功能更为强大。

`pd.concat`可以简单地合并一维的对象：


```python
ser1 = pd.Series(['A', 'B', 'C'], index=[1, 2, 3])
ser2 = pd.Series(['D', 'E', 'F'], index=[4, 5, 6])
pd.concat([ser1, ser2])
```




    1    A
    2    B
    3    C
    4    D
    5    E
    6    F
    dtype: object



也可以合并高维数据：


```python
df1 = make_df('AB', [1, 2])
df2 = make_df('AB', [3, 4])
display('df1', 'df2', 'pd.concat([df1, df2])')
```

![合并高维数据](http://oo3g995ih.bkt.clouddn.com/blog/180819/2LbadeA8jh.png?imageslim)

默认情况下，DataFrame的合并设置是`axis=0`，即按行合并。我们也可以进行设置：


```python
df3 = make_df('AB', [0, 1])
df4 = make_df('CD', [0, 1])
display('df3', 'df4', "pd.concat([df3, df4], axis=1)") # 使用axis='col'报错
```

![合并设置](http://oo3g995ih.bkt.clouddn.com/blog/180819/KKF0fggcik.png?imageslim)

#### 索引重复

`np.concatenate`和`pd.concat`的主要差异就在于`pd.concat`在合并时会**保留索引**，即使索引是重复的：


```python
x = make_df('AB', [0, 1])
y = make_df('AB', [2, 3])
y.index = x.index # 复制索引
display('x', 'y', 'pd.concat([x, y])')
```

![索引重复](http://oo3g995ih.bkt.clouddn.com/blog/180819/4BHg4mJLAc.png?imageslim)

虽然`pd.concat`允许保留相同的索引，但有时候我们并不想要这样做，`pd.concat`也提供了一些方法来解决这个问题：

- 捕捉索引重复。可使用`verify_integrity`参数来检查合并的结果中是否出现了重复的索引：


```python
try:
    pd.concat([x, y], verify_integrity=True)
except ValueError as e:
    print("ValueError:", e)
```

    ValueError: Indexes have overlapping values: [0, 1]
    

- 忽略索引。有时候索引无关紧要，可使用`ignore_index`参数，在合并是忽略索引：


```python
display('x', 'y', 'pd.concat([x, y], ignore_index=True)')
```

![忽略索引](http://oo3g995ih.bkt.clouddn.com/blog/180819/f6EJ5L4Ec4.png?imageslim)

- 增加多级索引。通过`keys`参数来设置多级索引标签：


```python
display('x', 'y', "pd.concat([x, y], keys=['x', 'y'])")
```

![设置多级索引](http://oo3g995ih.bkt.clouddn.com/blog/180819/2mcb85mF01.png?imageslim)

#### 类似 join 的合并

前面的例子中，合并的数据都有着相同的列名。而在实际工作中，需要合并的数据集往往有着不同的列名，`pd.concat`也提供了一些选项来解决这些问题：


```python
df5 = make_df('ABC', [1, 2])
df6 = make_df('BCD', [3, 4])
display('df5', 'df6', 'pd.concat([df5, df6])')
```

![不同列名的数据集](http://oo3g995ih.bkt.clouddn.com/blog/180819/iF38l1aDcD.png?imageslim)

默认情况下，某个位置缺失的数据会用`NaN`表示。如果不想这样可以使用``join``和``join_axes``参数来设置合并方式。默认的合并方式是并集合并(``join='outer'``)，我们也可以使用`join='inner'`来进行交集合并：


```python
display('df5', 'df6',
        "pd.concat([df5, df6], join='inner')")
```

![交集合并](http://oo3g995ih.bkt.clouddn.com/blog/180819/4k82c2558g.png?imageslim)

另一种合并方式是使用`join_axes`参数直接确定结果使用的列名：


```python
display('df5', 'df6',
        "pd.concat([df5, df6], join_axes=[df5.columns])")
```

![直接确定列名](http://oo3g995ih.bkt.clouddn.com/blog/180819/9FJJ3aGcED.png?imageslim)

#### append方法

可以直接使用`append()`方法来对两个数组直接合并：


```python
display('df1', 'df2', 'df1.append(df2)')
```

![append方法](http://oo3g995ih.bkt.clouddn.com/blog/180819/8dA7km06hf.png?imageslim)

需要注意的是与Python列表中的``append()``或``extend()``方法不同，Pandas的``append()``方法不会直接更新原有的值，而是创建一个新的对象，所以如果需要进行多个`append`操作，这样的效率就会比较低，建议先创建一个DataFrame列表，再用`concat()`进行合并。

## 合并数据集：合并与连接

### 数据连接的类型

`pd.merge()`函数可实现可实现三种数据连接：一对一、多对一和多对多。

#### 一对一连接

一对一连接和之间介绍的按列合并十分类似，如下示例：


```python
# 同一所公司员工不同信息的DataFrame

df1 = pd.DataFrame({'employee': ['Bob', 'Jake', 'Lisa', 'Sue'],
                    'group': ['Accounting', 'Engineering', 'Engineering', 'HR']})
df2 = pd.DataFrame({'employee': ['Lisa', 'Bob', 'Jake', 'Sue'],
                    'hire_date': [2004, 2008, 2012, 2014]})
display('df1', 'df2')
```

![示例数据](http://oo3g995ih.bkt.clouddn.com/blog/180819/6jic31B3L2.png?imageslim)

为了将两个DataFrame合并成一个，可以用`pd.merge()`函数实现：


```python
df3 = pd.merge(df1, df2)
df3
```

![pd.merge()函数](http://oo3g995ih.bkt.clouddn.com/blog/180819/Lj8BEAkIj3.png?imageslim)

`pd.merge()`方法会发现两个数据集都包含了`employee`列，会自动以这列作为键进行合并。需要注意的是，`pd.merge()`会默认丢弃原来的行索引，不过也可以自定义（之后会提及）。

#### 多对一连接

多对一连接是指，在需要连接的两列中有一列的值有重复，而多对一连接则会保留重复值：


```python
df4 = pd.DataFrame({'group': ['Accounting', 'Engineering', 'HR'],
                    'supervisor': ['Carly', 'Guido', 'Steve']})
display('df3', 'df4', 'pd.merge(df3, df4)')
```

![示例数据](http://oo3g995ih.bkt.clouddn.com/blog/180819/kIaIeHHBB8.png?imageslim)

![多对一连接](http://oo3g995ih.bkt.clouddn.com/blog/180819/9bfhHKicm6.png?imageslim)

#### 多对多连接

多对多连接就是输入的两列都包含重复值：


```python
df5 = pd.DataFrame({'group': ['Accounting', 'Accounting',
                              'Engineering', 'Engineering', 'HR', 'HR'],
                    'skills': ['math', 'spreadsheets', 'coding', 'linux',
                               'spreadsheets', 'organization']})
display('df1', 'df5', "pd.merge(df1, df5)")
```

![示例数据](http://oo3g995ih.bkt.clouddn.com/blog/180819/bHKkHLeL72.png?imageslim)

![多对多连接](http://oo3g995ih.bkt.clouddn.com/blog/180819/LDlam3hhdm.png?imageslim)

### 设置数据合并的键

`pd.merge()`默认会将两个数据集的一个或多个共同的列作为键进行合并，但通常要合并的列不是同名的，所以`pd.merge()`也提供了一些参数来处理这个问题。

#### 参数 on 的用法

这个参数只有在两个数据集有共同的列名才可使用，可将`on`参数设置为列名字符串或包含多个列名的列表：


```python
display('df1', 'df2', "pd.merge(df1, df2, on='employee')")
```

![on 的用法](http://oo3g995ih.bkt.clouddn.com/blog/180819/jFbcdFcka5.png?imageslim)

#### left_on 与 right_on 参数

有时候也需要合并两个不同列名的数据集，例如前面的数据集中的那一列不叫`employee`而是`name`，这时就可以用`left_on`与`right_on`参数来指定列名：


```python
df3 = pd.DataFrame({'name': ['Bob', 'Jake', 'Lisa', 'Sue'],
                    'salary': [70000, 80000, 120000, 90000]})
display('df1', 'df3', 'pd.merge(df1, df3, left_on="employee", right_on="name")')
```

![left_on 与 right_on 参数](http://oo3g995ih.bkt.clouddn.com/blog/180819/fLd3eGfhGC.png?imageslim)

获取的结果中会产生一个多余的列，可用`drop()`方法将这列去掉：


```python
pd.merge(df1, df3, left_on="employee", right_on="name").drop('name', axis=1)
```

![去掉重复的列](http://oo3g995ih.bkt.clouddn.com/blog/180819/kc8dLHmC59.png?imageslim)

#### left_index 与 right_index 参数

除了合并列以外，你可能还需要合并索引：


```python
df1a = df1.set_index('employee')
df2a = df2.set_index('employee')
display('df1a', 'df2a')
```

![示例数据](http://oo3g995ih.bkt.clouddn.com/blog/180819/kKb31A63FK.png?imageslim)

可以通过`left_index`和`right_index`参数将所有设置为键进行合并：


```python
display('df1a', 'df2a',
        "pd.merge(df1a, df2a, left_index=True, right_index=True)")
```

![合并索引](http://oo3g995ih.bkt.clouddn.com/blog/180819/f235C1LciL.png?imageslim)

还有一种更方便的方法就是用`join()`方法：


```python
display('df1a', 'df2a', 'df1a.join(df2a)')
```

![join 方法](http://oo3g995ih.bkt.clouddn.com/blog/180819/455KKLd4Gi.png?imageslim)

如果想将索引与列混合使用，则可以结合``left_index``和``right_on``或``left_on``和``right_index``：


```python
display('df1a', 'df3', "pd.merge(df1a, df3, left_index=True, right_on='name')")
```

![根据左边的索引和右边的列进行合并](http://oo3g995ih.bkt.clouddn.com/blog/180819/amkABcdBlC.png?imageslim)

### 设置数据连接的集合操作规则

看看下面的例子：


```python
df6 = pd.DataFrame({'name': ['Peter', 'Paul', 'Mary'],
                    'food': ['fish', 'beans', 'bread']},
                   columns=['name', 'food'])
df7 = pd.DataFrame({'name': ['Mary', 'Joseph'],
                    'drink': ['wine', 'beer']},
                   columns=['name', 'drink'])
display('df6', 'df7', 'pd.merge(df6, df7)')
```

![示例数据](http://oo3g995ih.bkt.clouddn.com/blog/180819/ag1CALl8km.png?imageslim)

我们合并的两个数据集在`name`列中只有一个共同的值，默认情况下，结果中只包含两个输入的**交集**，这种连接方式称为内连接（inner join）。我们可以用`how`参数设置连接方式，默认值为`inner`：


```python
pd.merge(df6, df7, how='inner')
```

![取交集](http://oo3g995ih.bkt.clouddn.com/blog/180819/e3hE9EChE6.png?imageslim)

`how`参数支持的连接方式还有`outer`、`left`以及`right`。外连接（outer join）返回两个数据集的并集，所有的缺失值都用`NaN`填充：


```python
display('df6', 'df7', "pd.merge(df6, df7, how='outer')")
```

![取并集](http://oo3g995ih.bkt.clouddn.com/blog/180819/E7f45I4kld.png?imageslim)

左连接和右连接返回的结果则只包含左边的列和右边的列：


```python
display('df6', 'df7', "pd.merge(df6, df7, how='left')")
```

![只包含左边的列](http://oo3g995ih.bkt.clouddn.com/blog/180819/0lkIcALLFK.png?imageslim)

```python
display('df6', 'df7', "pd.merge(df6, df7, how='right')")
```

![只包含右边的列](http://oo3g995ih.bkt.clouddn.com/blog/180819/EK696E90L0.png?imageslim)

### 重复列名：suffixes 参数

有时候两个输入的数据集包含着同名的列：


```python
df8 = pd.DataFrame({'name': ['Bob', 'Jake', 'Lisa', 'Sue'],
                    'rank': [1, 2, 3, 4]})
df9 = pd.DataFrame({'name': ['Bob', 'Jake', 'Lisa', 'Sue'],
                    'rank': [3, 1, 4, 2]})
display('df8', 'df9', 'pd.merge(df8, df9, on="name")')
```

![重复列名](http://oo3g995ih.bkt.clouddn.com/blog/180819/IJBei5fbJi.png?imageslim)

由于输出结果中有两个重复的列名，因此`pd.merge()`函数会自动添加后缀`_x`或`_y`，我们也可以通过`suffixes`参数来自定义后缀：


```python
display('df8', 'df9', 'pd.merge(df8, df9, on="name", suffixes=["_L", "_R"])')
```

![自定义后缀](http://oo3g995ih.bkt.clouddn.com/blog/180819/hG0hBAK5jk.png?imageslim)

更多关于Pandas的合并操作参考 Pandas 官方文档：[Merge, Join, and Concatenate](http://pandas.pydata.org/pandas-docs/stable/merging.html)

## 实例：美国各州的统计数据

- 数据下载：http://github.com/jakevdp/data-USstates/

首先用 Pandas 的`read_csv()`函数来读入数据集：


```python
pop = pd.read_csv('state-population.csv')
areas = pd.read_csv('state-areas.csv')
abbrevs = pd.read_csv('state-abbrevs.csv')

display('pop.head()', 'areas.head()', 'abbrevs.head()')
```

![读入数据集](http://oo3g995ih.bkt.clouddn.com/blog/180819/5elegfIcKD.png?imageslim)

现在我们想要获取美国各州的人口密度排名，首先需要将`pop`中的`state/region`列与`abbrevs`中的`abbreviation`列进行合并，来获取各州名称对应的全称，除此之外还需要用`how='outer'`来确保信息没有丢失：


```python
merged = pd.merge(pop, abbrevs, how='outer',
                  left_on='state/region', right_on='abbreviation')
merged = merged.drop('abbreviation', axis=1) #去掉 abbreviation 列
merged.head()
```

![merged](http://oo3g995ih.bkt.clouddn.com/blog/180819/0e2C5kGei5.png?imageslim)

接下来检查一下数据是否有缺失：


```python
merged.isnull().any()
```




    state/region    False
    ages            False
    year            False
    population       True
    state            True
    dtype: bool



可以发现部分`population`是缺失值，看看这些数据：


```python
merged[merged['population'].isnull()].head()
```

![查看缺失值](http://oo3g995ih.bkt.clouddn.com/blog/180819/4GgeCaK3bF.png?imageslim)

好像所有的人口缺失值都出现在2000年之前的波多黎各。

我们也发现一些州的数据也有缺失，看看哪些州有缺失：

- 用法参见 [Pandas | 数值取值与选择及数值运算方法](https://zwbao.github.io/#/Python/数据科学/Pandas数据处理（二）)


```python
merged.loc[merged['state'].isnull(), 'state/region'].unique()
```




    array(['PR', 'USA'], dtype=object)



现在我们来快速地填充对应的全称：


```python
merged.loc[merged['state/region'] == 'PR', 'state'] = 'Puerto Rico'
merged.loc[merged['state/region'] == 'USA', 'state'] = 'United States'
merged.isnull().any()
```




    state/region    False
    ages            False
    year            False
    population       True
    state           False
    dtype: bool



可以看到`state`列已经没有缺失值了。接下来可以把面积的数据也合并进来：


```python
final = pd.merge(merged, areas, on='state', how='left')
final.head()
```

![合并面积数据](http://oo3g995ih.bkt.clouddn.com/blog/180819/EG50k14cdf.png?imageslim)

再看看哪些列有缺失值：


```python
final.isnull().any()
```




    state/region     False
    ages             False
    year             False
    population        True
    state            False
    area (sq. mi)     True
    dtype: bool



可以看到`area`列还有缺失值：


```python
final.loc[final['area (sq. mi)'].isnull(), 'state'].unique()
```




    array(['United States'], dtype=object)



相同的操作还可以用掩码来完成：

- 用法参见 [NumPy | 将布尔数组作为掩码](https://zwbao.github.io/#/Python/数据科学/NumPy入门（三）?id=%e6%af%94%e8%be%83%e3%80%81%e6%8e%a9%e7%a0%81%e5%92%8c%e5%b8%83%e5%b0%94%e9%80%bb%e8%be%91)


```python
final['state'][final['area (sq. mi)'].isnull()].unique()
```




    array(['United States'], dtype=object)



我们发现面积的数据不包括全美国的面积数据，他也不在我们计算范围之内，所以需要去掉这个缺失值：


```python
final.dropna(inplace=True) # 使用inplace参数直接修改final的值
final.head()
```

![去掉缺失值](http://oo3g995ih.bkt.clouddn.com/blog/180819/I17hLcKbEK.png?imageslim)

现在所有的数据都已经准备好了。先选择2000年各州的人口及总人口数据，再用`query()`函数进行快速计算：


```python
data2010 = final.query("year == 2010 & ages == 'total'")
data2010.head()
```

![用 query 筛选数据](http://oo3g995ih.bkt.clouddn.com/blog/180819/LcchfCBada.png?imageslim)

计算人口密度并按序排列。首先需要对索引进行重置，然后再计算结果：


```python
data2010.set_index('state', inplace=True)
density = data2010['population'] / data2010['area (sq. mi)']
```


```python
density.sort_values(ascending=False, inplace=True) # ascending默认为True进行升序排列
density.head()
```




    state
    District of Columbia    8898.897059
    Puerto Rico             1058.665149
    New Jersey              1009.253268
    Rhode Island             681.339159
    Connecticut              645.600649
    dtype: float64



我们可以发现人口密度最高的地区是华盛顿特区的哥伦比亚地区，在各州的人口密度中，新泽西州是最高的。还可以看看人口密度最低的几个州的数据：


```python
density.tail()
```




    state
    South Dakota    10.583512
    North Dakota     9.537565
    Montana          6.736171
    Wyoming          5.768079
    Alaska           1.087509
    dtype: float64



可以看到，人口密度最低的州是阿拉斯加，刚刚超过1万人/平方公里。

> 笔记整理自《Python数据科学手册》，本书的英文版以及一些资料已在[Github开源](https://github.com/jakevdp/PythonDataScienceHandbook)（https://github.com/jakevdp/PythonDataScienceHandbook ）。
