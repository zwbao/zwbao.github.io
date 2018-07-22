# NumPy入门（三）


## 比较、掩码和布尔逻辑

### 示例：统计下雨天数

用`Pandas`加载2014年西雅图市的日降水统计数据：


```python
import numpy as np
import pandas as pd

rainfall = pd.read_csv('data/Seattle2014.csv')['PRCP'].values
inches = rainfall / 254 # 1/10mm -> inches
inches.shape
```




    (365,)



可见，这个数组包含了365个值，用`Matplotlib`生成下雨天数的直方图：


```python
%matplotlib inline
import matplotlib.pyplot as plt
import seaborn; seaborn.set() # 设置绘图风格
```


```python
plt.hist(inches, 40);
```


![直方图](http://oo3g995ih.bkt.clouddn.com/blog/180722/3fC1a9cF63.png?imageslim)


这幅直方图并不能很好地传递一些信息，比如一年中有多少天在下雨，下雨天的平均降水量是多少，有多少天的降水量超过了半英寸。

#### 深入数据

为了回答上面的一些问题，传统的方法是对数据循环计数，效率很低；但因为我们已学习了`NumPy`，就可以用通用函数来替代循环，以快速地实现数组的逐元素运算。

### 和通用元素类似的比较操作

`NumPy`除了实现了对数组逐元素的运算操作外，还实现了比如小于和大于的逐元素比较的通用函数：


```python
x = np.array([1, 2, 3, 4, 5])
x < 3
```




    array([ True,  True, False, False, False])




```python
x > 3
```




    array([False, False, False,  True,  True])




```python
x <= 3
```




    array([ True,  True,  True, False, False])




```python
x >= 3
```




    array([False, False,  True,  True,  True])




```python
x != 3
```




    array([ True,  True, False,  True,  True])




```python
x == 3
```




    array([False, False,  True, False, False])



对两个数组的逐元素比较：


```python
(2 * x) == (x ** 2)
```




    array([False,  True, False, False, False])



和算术运算符一样，比较运算符在`NumPy`中也是借助通用函数来实现的：

| Operator	    | Equivalent ufunc    || Operator	   | Equivalent ufunc    |
|---------------|---------------------||---------------|---------------------|
|``==``         |``np.equal``         ||``!=``         |``np.not_equal``     |
|``<``          |``np.less``          ||``<=``         |``np.less_equal``    |
|``>``          |``np.greater``       ||``>=``         |``np.greater_equal`` |

比较运算符也可用于任意形状、大小的数组：


```python
rng = np.random.RandomState(0)
x = rng.randint(10, size=(3, 4))
x
```




    array([[5, 0, 3, 3],
           [7, 9, 3, 5],
           [2, 4, 7, 6]])




```python
x < 6
```




    array([[ True,  True,  True,  True],
           [False, False,  True,  True],
           [ True,  True, False, False]])



### 操作布尔数组

#### 统计记录的个数

如需要统计布尔数组中`True`记录的个数，可以使用``np.count_nonzero`` 函数：


```python
np.count_nonzero(x < 6)
```




    8



可以看到有8条记录是小于6的。另一种实现方式是`np.sum`，``False`` 会被解释为 ``0``，``True`` 会被解释为 ``1``：


```python
np.sum(x < 6)
```




    8



用``sum()``的好处是，这个求和可以按行或列进行：


```python
# 每行有多少小于6的？
np.sum(x < 6, axis=1)
```




    array([4, 2, 2])



这是矩阵每一行小于6的个数。

快速检查任意或所有值是否为`True`， 可以用``np.any`` 或 ``np.all` ：


```python
np.any(x > 8)
```




    True




```python
np.any(x < 0)
```




    False




```python
np.all(x < 10)
```




    True




```python
np.all(x == 6)
```




    False



也可以沿特定坐标：


```python
np.all(x < 8, axis=1)
```




    array([ True, False,  True])



#### 布尔运算符

统计降水量在0.5英寸和1英寸之间的天数：


```python
np.sum((inches > 0.5) & (inches < 1))
```




    29




```python
# 用“非”的方式来算
np.sum(~((inches <= 0.5) | (inches >= 1)))
```




    29



布尔运算符及其对应的通用函数：

| Operator	    | Equivalent ufunc    || Operator	    | Equivalent ufunc    |
|---------------|---------------------||---------------|---------------------|
|``&``          |``np.bitwise_and``   ||&#124;         |``np.bitwise_or``    |
|``^``          |``np.bitwise_xor``   ||``~``          |``np.bitwise_not``   |

利用这些工具就可以回答上面关于天气的一些问题：


```python
print("Number days without rain:      ", np.sum(inches == 0))
print("Number days with rain:         ", np.sum(inches != 0))
print("Days with more than 0.5 inches:", np.sum(inches > 0.5))
print("Rainy days with < 0.2 inches  :", np.sum((inches > 0) &
                                                (inches < 0.2)))
```

    Number days without rain:       215
    Number days with rain:          150
    Days with more than 0.5 inches: 37
    Rainy days with < 0.2 inches  : 75
    

### 将布尔数组作为掩码


```python
x
```




    array([[5, 0, 3, 3],
           [7, 9, 3, 5],
           [2, 4, 7, 6]])



用前面介绍的方法可以用得到一个布尔数组：


```python
x < 5
```




    array([[False,  True,  True,  True],
           [False, False,  True, False],
           [ True,  True, False, False]])



我们可以通过掩码操作将这些值从数组中选出：


```python
x[x < 5]
```




    array([0, 3, 3, 3, 2, 4])



试着对西雅图降水数据进行一些统计：


```python
# 为所有下雨天创建掩码
rainy = (inches > 0)

# 构建一个包括整个夏天日期的掩码（6月21日为第172天）
summer = ( np.arange(365) - 172 < 90) & ( np.arange(365) -172 > 0)

print("Median precip on rainy days in 2014 (inches):   ",
      np.median(inches[rainy]))
print("Median precip on summer days in 2014 (inches):  ",
      np.median(inches[summer]))
print("Maximum precip on summer days in 2014 (inches): ",
      np.max(inches[summer]))
print("Median precip on non-summer rainy days (inches):",
      np.median(inches[rainy & ~summer]))
```

    Median precip on rainy days in 2014 (inches):    0.19488188976377951
    Median precip on summer days in 2014 (inches):   0.0
    Maximum precip on summer days in 2014 (inches):  0.8503937007874016
    Median precip on non-summer rainy days (inches): 0.20078740157480315
    

## 花哨的索引

简而言之就是构建一个索引数组来获得多个数组元素，例如：


```python
rand = np.random.RandomState(42)
x = rand.randint(100, size=10)
print(x)
```

    [51 92 14 71 60 20 82 86 74 74]
    

当我们需要其中三个元素，可以通过以下方法实现：


```python
[x[3], x[7], x[2]]
```




    [71, 86, 14]



另一种方法就是通过传递索引来获得：


```python
ind = [3, 7, 4]
x[ind]
```




    array([71, 86, 60])



索引的结果形状与索引数组的形状一致：


```python
ind = np.array([[3, 7],
                [4, 5]])
x[ind]
```




    array([[71, 86],
           [60, 20]])



也对多个维度适用：


```python
X = np.arange(12).reshape((3, 4))
X
```




    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11]])



第一个索引指的是行，第二个索引指的是列：


```python
row = np.array([0, 1, 2])
col = np.array([2, 1, 3])
X[row, col]
```




    array([ 2,  5, 11])



索引值的配对遵循广播的规则。需要特别注意的是，索引返回的值是广播后的索引数组的形状，而不是被索引的数组形状。

### 组合索引

可以将花哨的索引和简单的索引组合使用：


```python
X[2, [2, 0, 1]]
```




    array([10,  8,  9])



也可以和切片组合使用：


```python
X[1:, [2, 0, 1]]
```




    array([[ 6,  4,  5],
           [10,  8,  9]])



更可以和掩码组合使用：


```python
mask = np.array([1, 0, 1, 0], dtype=bool)
X[row[:, np.newaxis], mask]
```




    array([[ 0,  2],
           [ 4,  6],
           [ 8, 10]])



### 用花哨的索引修改值



```python
x = np.arange(10)
i = np.array([2, 1, 8, 4])
x[i] = 99
print(x)
```

    [ 0 99 99  3 99  5  6  7 99  9]
    


```python
x[i] -= 10
print(x)
```

    [ 0 89 89  3 89  5  6  7 89  9]
    

不过，操作中重复的索引也会导致一些出乎意料的结果产生：


```python
x = np.zeros(10)
x[[0, 0]] = [4, 6]
print(x)
```

    [6. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    

实际上是先赋值``x[0] = 4``，然后赋值 ``x[0] = 6``。

## 数组的排序

### NumPy中的快速排序： ``np.sort`` 和 ``np.argsort``

如果想在不修改原始输入数组的基础上返回一个排好序的数组，可以使用 ``np.sort``：


```python
x = np.array([2, 1, 4, 3, 5])
np.sort(x)
```




    array([1, 2, 3, 4, 5])



如果希望用排好序的数组替代原始数组，可用数组的`sort`方法：


```python
x.sort()
print(x)
```

    [1 2 3 4 5]
    

另一个相关函数是``np.argsort``， 返回原始数组排好序的索引值：


```python
x = np.array([2, 1, 4, 3, 5])
i = np.argsort(x)
print(i)
```

    [1 0 3 2 4]
    

这些索引可以用于创建有序的数组（花哨的索引）：


```python
x[i]
```




    array([1, 2, 3, 4, 5])



#### 沿着行或列排序

可通过`axis`参数，沿着行或列进行排序：


```python
rand = np.random.RandomState(42)
X = rand.randint(0, 10, (4, 6))
print(X)
```

    [[6 3 7 4 6 9]
     [2 6 7 4 3 7]
     [7 2 5 4 1 7]
     [5 1 4 0 9 5]]
    


```python
# 对X的每一列排序
np.sort(X, axis=0)
```




    array([[2, 1, 4, 0, 1, 5],
           [5, 2, 5, 4, 3, 7],
           [6, 3, 7, 4, 6, 7],
           [7, 6, 7, 4, 9, 9]])




```python
# 对X的每一行排序
np.sort(X, axis=1)
```




    array([[3, 4, 6, 6, 7, 9],
           [2, 3, 4, 6, 7, 7],
           [1, 2, 4, 5, 7, 7],
           [0, 1, 4, 5, 5, 9]])



**需要注意的是，这会把行或列当成独立的数组，任何行或列的值之间的关系会丢失！**

### 部分排序：分隔

有时候我们并不希望对整个数组进行排序，仅仅想找到数组中第*K*小的值，可以使用``np.partition`` 。该函数的输入时数组和数字*K*，输出是一个新数组：


```python
x = np.array([7, 2, 3, 1, 6, 5, 4])
np.partition(x, 3)
```




    array([2, 1, 3, 4, 6, 5, 7])



以排序后的第3个数，即3进行分区，分区后的结果即是：小于3的元素2,1位于3的前面，大于等于3的元素都位于3的后面（任意顺序）。

与排序类似，也可以按任意轴进行分隔：


```python
np.partition(X, 2, axis=1)
```




    array([[3, 4, 6, 7, 6, 9],
           [2, 3, 4, 7, 6, 7],
           [1, 2, 4, 5, 7, 7],
           [0, 1, 4, 5, 9, 5]])




> 笔记整理自《Python数据科学手册》，本书的英文版以及一些资料已在[Github开源](https://github.com/jakevdp/PythonDataScienceHandbook)（https://github.com/jakevdp/PythonDataScienceHandbook ）。
