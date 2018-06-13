# 我的Python笔记·模块化编程（二）

## 问题

1856年，孟德尔开始了长达8年的豌豆实验，从此奠定了遗传学的基础。示例程序将使用类`class`来模拟豌豆杂交实验。显性等位基因（黄色）和隐性基因（绿色）分别由`G`和`g`来表示。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180604/56B9bB8Ikh.png?imageslim)

```python
class Pea:

    def __init__(self, genotype):
        self.genotype = genotype
    def get_phenotype(self):
        if "G" in self.genotype:
            return "yellow"
        else:
            return "green"

    def create_offspring(self, other):
        offspring = []
        new_genotype = ""
        for haplol in self.genotype:
            for haplol2 in other.genotype:
                new_genotype = haplol + haplol2
                offspring.append(Pea(new_genotype))
        return offspring

    def __repr__(self):
        return self.get_phenotype() + '[%s]' % self.genotype

yellow = Pea("GG")
green = Pea("gg")
f1 = yellow.create_offspring(green)
f2 = f1[0].create_offspring(f1[1])

print(f1)
print(f2)
```

- 程序输出为：

```python
[yellow[Gg], yellow[Gg], yellow[Gg], yellow[Gg]]
[yellow[GG], yellow[Gg], yellow[gG], green[gg]]
```

## 程序是如何工作的

`Pea`类定义了每颗豌豆的属性和行为。每颗豌豆都有自己的基因型（`GG`，`Gg`，`gg`），和表型（`yellow` or `green`）。最后，每颗豌豆都可以与第二颗豌豆生成后代。

### 类

众所周知，Python是一种面向对象的语言，一个对象的特征称之为属性（attribute），它所具有的行为称之为方法（method）。所以在python中，我们可以把具有相同属性和方法的对象归为一个类（class）。

- 对象=属性+方法

#### 创建类

**类**定义了其所代表事物的特征和行为，但不包含特定的数据，当在这一类中赋予特定的具体数据时则称之为**实例**。使用`class`来创建一个新类，`class`之后为类的名称并以冒号结尾。下面一整个缩进的代码块都属于该类。

```python
class Pea:
```

#### 构造函数`__init__`

构造函数`__init__`（注意`init`前后各有两条下划线）具有初始化的作用，也就是当该类被实例化的时候就会执行该函数。我们可以把要先初始化的属性放到这个函数里面。

```python
class Pea:
    def __init__(self, genotype):
        self.genotype = genotype
```

比如，通过赋值`self.genotype = genotype`使每个豌豆都有其特定的基因型（也可以称之为属性），`__init__`方法的第一个参数永远是`self`，表示创建的类实例本身，因此，在`__init__`方法内部，就可以把各种属性绑定到`self`，因为`self`就指向创建的实例本身。

这里`self`就是指类本身，`self.genotype`就是`Pea`类的属性变量，是`Pea`类所有。而`genotype`是外部传来的参数，不是`Pea`类所自带的。所以`self.genotype = genotype`的意思就是把外部传来的参数`genotype`的值赋值给`Pea`类的属性变量`self.genotype`。

#### 创建实例

```python
yellow = Pea("GG")
print(yellow.genotype)
```

创建实例时与调用函数类似，将构造函数所需要的参数`genotype`作为参数传递。`self`只是告诉类，自己由哪个实例调用。本例中yellow是调用Pea类的实例（变量），GG则为yellow实例中genotype的具体值。

#### 类以属性的形式包含数据

实例内的数据储存在属性中。可以通过`.`来读取（如`yellow.genotype`），属性也可以像变量一样动态地改变，例如重新赋值更新`yellow.genotype = 'Gg'`

#### 类包含的方法

类中的函数称之为**方法**，在`Pea`类中包含了计算表型的方法：

```python
    def get_phenotype(self):
        if "G" in self.genotype:
            return "yellow"
        else:
            return "green"
```

`get_phenotype()`方法使用了`genotype`属性中的数据，`self`则指实例本身。

#### `__repr__`方法打印类和实例

打印从类中创建的对象时，Python通常会显示这样的信息：

```python
<Pea object a FFFFx234234ou>
```

这里给出的信息很不充分，这时候就可以在类中添加一种称为`__repr__()`的特殊方法来打印对象：

```python
def __repr__(self):
    return self.get_phenotype() + '[%s]' % self.genotype
```

上面的方法会返回一个包含实例基因型和表型的字符串。除`self`参数以外`__repr__()`方法不需要其他参数。打印`Pea`实例时，`__repr__()`会被自动调用。编写`__repr__()`方法时，只需要包括你认为必要的信息即可。

#### 创建子类

一个类的属性和方法可以由其他类继承。被继承的类称为父类，继承的类称为自类。比如想要豌豆包含注释，则可以定义一个从`Pea`类继承的`CommentedPea`类：

```python
class CommentedPea(Pea):

    def __init__(self, genotype, comment):
        Pea.__init__(self, genotype)
        self.comment = comment

    def __repr__(self):
        return '%s [%s] (%s)' % (self.get_phenotype(),
        self.genotype, self.comment)

yellow1 = CommentedPea('GG', 'homozygote')
yellow2 = CommentedPea('Gg', 'heterozygote')
print(yellow1)
```

- 程序输出为：

```python
yellow [GG] (homozygote)
```

除了新加入的，`CommentedPea`类与`Pea`类的方法和属性一样，注意`__init__()`方法用于调用父方法。
