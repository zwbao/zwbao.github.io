# NumPy入门（二）

## Numpy数组的计算：通用函数



```python
import numpy as np
```

## NumPy的通用函数

### 数组的运算

NumPy继承了Python原生的算术运算符，除了加减乘除以外还包括指数和取模等操作：


```python
x = np.arange(4)
print("x     =", x)
print("x + 5 =", x + 5)
print("x - 5 =", x - 5)
print("x * 2 =", x * 2)
print("x / 2 =", x / 2)
print("x // 2 =", x // 2)  # floor division
print("-x     = ", -x)
print("x ** 2 = ", x ** 2)
print("x % 2  = ", x % 2)
```

    x     = [0 1 2 3]
    x + 5 = [5 6 7 8]
    x - 5 = [-5 -4 -3 -2]
    x * 2 = [0 2 4 6]
    x / 2 = [0.  0.5 1.  1.5]
    x // 2 = [0 0 1 1]
    -x     =  [ 0 -1 -2 -3]
    x ** 2 =  [0 1 4 9]
    x % 2  =  [0 1 0 1]
    

当然也可以对这些运算符任意组合（注意优先级）：


```python
-(0.5*x + 1) ** 2
```




    array([-1.  , -2.25, -4.  , -6.25])



所有运算符都是NumPy内置函数的简单封装器，比如 `np.add(x, 2)`就相当于`x + 2`，以下是所有NumPy实现的算术运算符：

| Operator	    | Equivalent ufunc    | Description                           |
|---------------|---------------------|---------------------------------------|
|``+``          |``np.add``           |Addition (e.g., ``1 + 1 = 2``)         |
|``-``          |``np.subtract``      |Subtraction (e.g., ``3 - 2 = 1``)      |
|``-``          |``np.negative``      |Unary negation (e.g., ``-2``)          |
|``*``          |``np.multiply``      |Multiplication (e.g., ``2 * 3 = 6``)   |
|``/``          |``np.divide``        |Division (e.g., ``3 / 2 = 1.5``)       |
|``//``         |``np.floor_divide``  |Floor division (e.g., ``3 // 2 = 1``)  |
|``**``         |``np.power``         |Exponentiation (e.g., ``2 ** 3 = 8``)  |
|``%``          |``np.mod``           |Modulus/remainder (e.g., ``9 % 4 = 1``)|

### 绝对值


```python
x = np.array([-2, -1, 0, 1, 2])
np.absolute(x)
```




    array([2, 1, 0, 1, 2])



该函数也可以使用`np.abs`：


```python
np.abs(x)
```




    array([2, 1, 0, 1, 2])



### 三角函数

首先定义一个角度数组：


```python
theta = np.linspace(0, np.pi, 3)
```

接下来就可以进行一些三角函数的计算啦~


```python
print("theta      = ", theta)
print("sin(theta) = ", np.sin(theta))
print("cos(theta) = ", np.cos(theta))
print("tan(theta) = ", np.tan(theta))
```

    theta      =  [0.         1.57079633 3.14159265]
    sin(theta) =  [0.0000000e+00 1.0000000e+00 1.2246468e-16]
    cos(theta) =  [ 1.000000e+00  6.123234e-17 -1.000000e+00]
    tan(theta) =  [ 0.00000000e+00  1.63312394e+16 -1.22464680e-16]
    

### 指数和对数


```python
x = [1, 2, 3]
print("x     =", x)
print("e^x   =", np.exp(x))
print("2^x   =", np.exp2(x))
print("3^x   =", np.power(3, x))
```

    x     = [1, 2, 3]
    e^x   = [ 2.71828183  7.3890561  20.08553692]
    2^x   = [2. 4. 8.]
    3^x   = [ 3  9 27]
    


```python
x = [1, 2, 4, 10]
print("x        =", x)
print("ln(x)    =", np.log(x))
print("log2(x)  =", np.log2(x))
print("log10(x) =", np.log10(x))
```

    x        = [1, 2, 4, 10]
    ln(x)    = [0.         0.69314718 1.38629436 2.30258509]
    log2(x)  = [0.         1.         2.         3.32192809]
    log10(x) = [0.         0.30103    0.60205999 1.        ]
    

## 高级的通用函数特性

### 指定输出

所有的通用函数都可以用`out`参数来指定计算结果的存放位置：


```python
x = np.arange(5)
y = np.empty(5)
np.multiply(x, 10, out=y)
print(y)
```

    [ 0. 10. 20. 30. 40.]
    

这个特性也可以用作数组视图，比如将计算结果写入指定数组的每个元素的位置：


```python
y = np.zeros(10)
np.power(2, x, y[::2])
print(y)
```

    [ 1.  0.  2.  0.  4.  0.  8.  0. 16.  0.]
    

对于计算量较大的数组，使用`out`可以有效节约内存。

### 聚合

`reduce`方法可以对给定的元素和操作重复执行，直到得到单个结果。例如。对`add`函数调用`reduce`方法会返回数组中所有函数的和：


```python
x = np.arange(1, 6)
np.add.reduce(x)
```




    15



同样，对于`multiply`函数调用`reduce`会返回数组中所有元素的乘积：


```python
np.multiply.reduce(x)
```




    120



如果需要储存每次计算的中间结果，可以使用`accumulate`：


```python
np.add.accumulate(x)
```




    array([ 1,  3,  6, 10, 15], dtype=int32)




```python
np.multiply.accumulate(x)
```




    array([  1,   2,   6,  24, 120], dtype=int32)



注意：在一些特殊的情况中，NumPy也提供了专用的函数来实现以上`reduce`的功能（例如：``np.sum``, ``np.prod``, ``np.cumsum``, ``np.cumprod``）。

### 外积

任何通用函数都可以用`outer`方法获得两个不同输入数组所有元素对函数的运算结果，比如结合`multiply`可以一行代码生成一个乘法表：


```python
x = np.arange(1, 6)
np.multiply.outer(x, x)
```




    array([[ 1,  2,  3,  4,  5],
           [ 2,  4,  6,  8, 10],
           [ 3,  6,  9, 12, 15],
           [ 4,  8, 12, 16, 20],
           [ 5, 10, 15, 20, 25]])



## 聚合：最小值、最大值和其他值

### 数组值求和

NumPy中的`sum`函数和Python本身的`sum`函数十分类似，结果也一样，但由于NumPy的`sum`函数在编译码中执行操作，所以NumPy会更快一些：


```python
big_array = np.random.rand(1000000)
%timeit sum(big_array)
%timeit np.sum(big_array)
```

    133 ms ± 1.01 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
    799 µs ± 5.32 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
    

需要注意的是，NumPy中的`sum`函数和Python本身的`sum`函数并不完全相同，`np.sum`函数是知道数组维度的。

### 最小值和最大值

Python也有内置的`min`函数和`max`函数，NumPy中也有对应的语法并且执行得更快：


```python
%timeit min(big_array)
%timeit np.min(big_array)
```

    51.3 ms ± 346 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)
    392 µs ± 5.52 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
    

还有一种更简单地调用这些函数的方法，就是数组直接调用这些方法：


```python
print(big_array.min(), big_array.max(), big_array.sum())
```

    4.86147444966889e-07 0.9999992378264679 499824.3720491337
    

#### 多维度聚合

当你的数据储存在二维数组中时，默认情况下`np.sum`会返回整个数组的聚合结果，调用`axis`参数可以指定沿哪个轴进行聚合。例如：通过`axis=0`参数找到每一列的最小值：


```python
M = np.random.random((3, 4))
print(M)
M.min(axis=0)
```

    [[0.62475337 0.3414042  0.92886579 0.69921013]
     [0.63710828 0.35998738 0.46732804 0.25822278]
     [0.69088793 0.78272269 0.81428264 0.31513276]]
    




    array([0.62475337, 0.3414042 , 0.46732804, 0.25822278])



函数会返回四个值，对应四列数字的计算值。

也可以找每一行的最大值：


```python
M.max(axis=1)
```




    array([0.92886579, 0.63710828, 0.81428264])



#### 其他聚合函数

大多数聚合在计算时都会对NaN值进行忽略。

以下是NumPy中可用的聚合函数：

|Function Name      |   NaN-safe Version  | Description                                   |
|-------------------|---------------------|-----------------------------------------------|
| ``np.sum``        | ``np.nansum``       | Compute sum of elements                       |
| ``np.prod``       | ``np.nanprod``      | Compute product of elements                   |
| ``np.mean``       | ``np.nanmean``      | Compute mean of elements                      |
| ``np.std``        | ``np.nanstd``       | Compute standard deviation                    |
| ``np.var``        | ``np.nanvar``       | Compute variance                              |
| ``np.min``        | ``np.nanmin``       | Find minimum value                            |
| ``np.max``        | ``np.nanmax``       | Find maximum value                            |
| ``np.argmin``     | ``np.nanargmin``    | Find index of minimum value                   |
| ``np.argmax``     | ``np.nanargmax``    | Find index of maximum value                   |
| ``np.median``     | ``np.nanmedian``    | Compute median of elements                    |
| ``np.percentile`` | ``np.nanpercentile``| Compute rank-based statistics of elements     |
| ``np.any``        | N/A                 | Evaluate whether any elements are true        |
| ``np.all``        | N/A                 | Evaluate whether all elements are true        |

### 示例：美国总统的身高是多少

我们使用`Pandas`包来读取文件并抽取身高信息（在之后的笔记中会详细介绍），`president_heights.csv`文件可在github（https://github.com/jakevdp/PythonDataScienceHandbook） 下载到。


```python
import pandas as pd
```


```python
data = pd.read_csv('./president_heights.csv')
heights = np.array(data['height(cm)'])
print(heights)
```

    [189 170 189 163 183 171 185 168 173 183 173 173 175 178 183 193 178 173
     174 183 183 168 170 178 182 180 183 178 182 188 175 179 183 193 182 183
     177 185 188 188 182 185]
    

得到身高的数组后，就可以利用NumPy中的函数来计算一些统计值了：


```python
print("Mean height:       ", heights.mean())
print("Standard deviation:", heights.std())
print("Minimum height:    ", heights.min())
print("Maximum height:    ", heights.max())
```

    Mean height:        179.73809523809524
    Standard deviation: 6.931843442745892
    Minimum height:     163
    Maximum height:     193
    

还可以计算分位数：


```python
print("25th percentile:   ", np.percentile(heights, 25))
print("Median:            ", np.median(heights))
print("75th percentile:   ", np.percentile(heights, 75))
```

    25th percentile:    174.25
    Median:             182.0
    75th percentile:    183.0
    

使用`Matplotlib`包进行快速画图（也会在之后的笔记中介绍）：


```python
# 内嵌画图
%matplotlib inline
import matplotlib.pyplot as plt
import seaborn; seaborn.set() # 设置绘图风格
```


```python
plt.hist(heights)
plt.title('Height Distribution of US Presidents')
plt.xlabel('height (cm)')
plt.ylabel('number');
```


![mark](http://oo3g995ih.bkt.clouddn.com/blog/180714/944gb9glAb.png?imageslim)


## 数组的计算：广播

### 广播的介绍

广播可以对不同大小的数组的对应元素进行逐个计算。比如，可以把一个标量和一个数组相加：


```python
a = np.array([0, 1, 2])
a + 5
```




    array([5, 6, 7])



接着，观察下一个一维数组和一个二维数组相加：


```python
M = np.ones((3, 3))
M
```




    array([[1., 1., 1.],
           [1., 1., 1.],
           [1., 1., 1.]])




```python
M + a
```




    array([[1., 2., 3.],
           [1., 2., 3.],
           [1., 2., 3.]])



为了方便理解，上面的例子可以可视化为下图：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180714/Kb3h93if37.png?imageslim)

### 广播的规则

- 如果两个数组的维度数不相同，那么小维度数组将会在其左边补1；
- 如果两个数组的形状在任何一个维度上都不匹配，那么数组的形状会沿着维度为1的维度扩展以匹配另外一个数组的形状；
- 如果两个数组的形状在任何一个维度上都不匹配并且没有一个维度等于1，那么会发生异常。

下面来看几个示例：

#### 示例1

将一个二维数组与一个一维数组相加：


```python
M = np.ones((2, 3))
a = np.arange(3)
```

这两个数组的形状如下：

- ``M.shape = (2, 3)``
- ``a.shape = (3,)``

根据第一条规则，数组 ``a`` 的维度更小，所以在其左边补上1：

- ``M.shape -> (2, 3)``
- ``a.shape -> (1, 3)``

根据第二条规则，第一个维度不匹配，所以需要扩展这个维度以匹配数组：

- ``M.shape -> (2, 3)``
- ``a.shape -> (2, 3)``

现在这两个数组的形状匹配了，他们的形状都为 ``(2, 3)``：


```python
M + a
```




    array([[1., 2., 3.],
           [1., 2., 3.]])



#### 示例2


```python
a = np.arange(3).reshape((3, 1))
b = np.arange(3)
```

可以看到这两个数组均需要进行广播，他们的形状为：

- ``a.shape = (3, 1)``
- ``b.shape = (3,)``

根据第一条规则，需要在数组 ``b``左边补上1：

- ``a.shape -> (3, 1)``
- ``b.shape -> (1, 3)``

根据第二条规则，需要扩展他们的维度：

- ``a.shape -> (3, 3)``
- ``b.shape -> (3, 3)``

可以看到以下结果：


```python
a + b
```




    array([[0, 1, 2],
           [1, 2, 3],
           [2, 3, 4]])



#### 示例3


```python
M = np.ones((3, 2))
a = np.arange(3)
```




    array([0, 1, 2])



这两个数组的形状如下：

- ``M.shape = (3, 2)``
- ``a.shape = (3,)``

同样，根据第一条规则在数组 ``a`` 左边补上1：

- ``M.shape -> (3, 2)``
- ``a.shape -> (1, 3)``

根据第二条规则，扩展数组 ``a`` 的第一个维度来匹配 ``M``:

- ``M.shape -> (3, 2)``
- ``a.shape -> (3, 3)``

现在，根据第三条规则，最终形状还是不匹配，所以会看到如下结果：


```python
M + a
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-60-8cac1d547906> in <module>()
    ----> 1 M + a
    

    ValueError: operands could not be broadcast together with shapes (3,2) (3,) 


### 广播的实际应用

#### 数组的归一化

假设，现在你有10个观测值，每个观测值包含3个数值：


```python
X = np.random.random((10, 3))
```

我们可以计算每个特征的均值，用`mean`函数沿着第一个维度（列）进行聚合：


```python
Xmean = X.mean(0)
Xmean
```




    array([0.65764505, 0.60113647, 0.68073558])



现在通过从X数组的元素中减去这个均值实现归一化：


```python
X_centered = X - Xmean
```

为了进一步核对处理是否正确，可以看看归一化的数组的均值是否接近0：


```python
X_centered.mean(0)
```




    array([-2.22044605e-17, -7.77156117e-17, -3.33066907e-17])



#### 画一个二维函数

广播另一个十分有用的地方就是，它可以基于二维函数显示图像：


```python
# x和y表示0~5区间50个步长的序列
x = np.linspace(0, 5, 50)
y = np.linspace(0, 5, 50)[:, np.newaxis] # 转换为列向量，可参见上一次的笔记

z = np.sin(x) ** 10 + np.cos(10 + y * x) * np.cos(x)
```

使用Matplotlib来画出这个二维数组：


```python
%matplotlib inline
import matplotlib.pyplot as plt
```


```python
plt.imshow(z, origin='lower', extent=[0, 5, 0, 5],
           cmap='viridis')
plt.colorbar();
```

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180714/a2iALE4CGk.png?imageslim)



> 笔记整理自《Python数据科学手册》，本书的英文版以及一些资料已在[Github开源](https://github.com/jakevdp/PythonDataScienceHandbook)（https://github.com/jakevdp/PythonDataScienceHandbook）。
