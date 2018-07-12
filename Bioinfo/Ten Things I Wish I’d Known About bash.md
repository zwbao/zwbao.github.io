> 原文链接：[Ten Things I Wish I’d Known About bash](https://zwischenzugs.com/2018/01/06/ten-things-i-wish-id-known-about-bash/)

\# 为我的注释

## \`\` vs  $()

这两种写法做了同一件事，比较这两行：

```
$ echo `ls`
$ echo $(ls)
```

为什么会存在着两者形式？这困扰了我好久。

这两种形式都将其中包含的命令替换为命令的执行结果。
\# 即命令替换，详见：[我的Linux笔记·bash的特性](https://vip.biotrainee.com/d/252 )

他们的主要区别在于哪一个嵌套更简单。

看看下面哪个更容易阅读（和写）？

```
$ echo `echo \`echo \\\`echo inside\\\`\``
$ echo $(echo $(echo $(echo inside)))
```

如果你还有兴趣深入，请看[这里](http://mywiki.wooledge.org/BashFAQ/082)或[这里](https://stackoverflow.com/questions/4708549/what-is-the-difference-between-command-and-command-in-shell-programming)。

## globbing vs regexps

这是另一个要是你从来没有想过或者探究过就会混淆的话题。

globs 和 regexps 看起来很像，但是他们不尽相同。

思考一下这个命令：

```
$ rename -n 's/(.*)/new$1/' *
```

两个星号意思并不一样。

第一个被shell忽略（因为它在引号中），并被理解为零个或更多个字符，所以它被解释为正则表达式。

第二个由shell解释（因为它不在引号中），并被当前工作文件夹中所有文件的列表替换。它被解释为glob。

所以通过`man bash`，你可以找出为什么这两个命令产生了不同的结果：

```
$ ls *
$ ls .*    # 以点开头的文件
```

第二个看起来更像一个正则表达式。但它不是！

\# globbing 的内容详见：[我的Linux笔记·bash的特性](https://vip.biotrainee.com/d/252 )
\# 正则表达式的内容详见：[我的Linux笔记·grep及正则表达式](https://vip.biotrainee.com/d/256 )

## 退出状态码

不是每个人都知道，每当你在bash中运行一个shell命令时，都会有一个“退出状态码”返回给bash。

通常，如果一个命令“成功”，你会得到一个`0`的退出状态码。如果它不成功，你会得到一个非零的代码。 `1`是一个“一般错误”，其他数字可以给你更多的信息（例如，哪个信号杀死了它）。

但这些规则并不总是成立：

```
$ grep not_there /dev/null
$ echo $？
```

`$?`是一个特殊的bash变量，它在运行后被设置为每个命令的退出状态码。

`grep`使用退出状态码来表示它是否匹配到了东西。我必须每次回头看看：是否匹配到了或者返回了不是`0`的值？

好好理解它，会使你在学习bash的过程中豁然开朗。

\# 关于退出状态码的内容详见：[我的Linux笔记·bash脚本编程](https://vip.biotrainee.com/d/257 )

## if statements, [ and [[

这是另一个“找不同”，就像上面的反引号。

看看接下来的命令将输出什么？

```
if grep not_there /dev/null
then
    echo hi
else
    echo lo
fi
```

作为其使用了退出码的副作用，grep的返回码使得像这样的代码更直观地工作。

下面是几个选项：

- a) `hihi`
- b) `lolo`
- c) 别的东西

```
if [ $(grep not_there /dev/null) = '' ]
then
    echo -n hi
else
    echo -n lo
fi
if [[ $(grep not_there /dev/null) = '' ]]
then
    echo -n hi
else
    echo -n lo
fi
```

`[`和`[[`之间的区别是另一个我从来没有真正理解的东西。 `[`是条件判断的原始形式，之后`[[`被引入，它更灵活和直观。在上面的第一个`if`块中，因为`$(grep not_there /dev/null)`的结果为空，导致了下面这种比较：

```
[ = '' ]
```

这是没有意义的。而双括号的形式可以解决这个问题。

这就是为什么你偶尔会在bash脚本中看到这样的写法：

```
if [ x$(grep not_there /dev/null) = 'x' ]
```

所以如果命令没有返回值，它仍然可以运行，而不会报错。虽然没有必要，但这就是它存在的原因。

## set

bash 有许多可配置的选项，我一直使用其中的两个：

```
set -e
```

任何命令返回一个非零的退出状态码（见上），就退出脚本。
\# `set -e`告诉bash如果任何语句的执行结果不是true则应该退出。这样的好处是防止错误像滚雪球般变大导致一个致命的错误，而这些错误本应该在之前就被处理掉。

这将输出在运行时所运行的命令：

```
set -x
```

所以一个脚本可能是这样开头的：

```
#!/bin/bash
set -e
set -x
grep not_there /dev/null
echo $?
```

那么这个脚本会输出什么？

## <()

这是我最喜欢的。但它从未被充分利用，也许是因为它可能最初看起来令人困惑，但我一直在使用它。

它和`$()`类似，里面的命令的输出被重新使用。

但是在这种情况下，该输出被视为一个文件。该文件又可以用作将文件作为参数的命令的参数。

是不是有些困惑？看看下面的例子：

你有没有做过这样的事情？

```
$ grep somestring file1 > /tmp/a
$ grep somestring file2 > /tmp/b
$ diff /tmp/a /tmp/b
```

这写的没错，但你也可以这样写：

```
diff <(grep somestring file1) <(grep somestring file2)
```

是不是看起来更简洁些？

## 引用

引用是bash中的一个棘手话题，就像在许多软件环境中一样。

首先，看看引号的作用：

```
A='123'
echo "$A"
echo '$A'
```

很简单 ：双引号引用变量，而单引号就是字面的意思。

那么试试看，下面的结果是什么？

```
mkdir -p tmp
cd tmp
touch a
echo "*"
echo '*'
```

是不是很惊讶？我也是。

\# 这两行确实就是输出星号，不过若是输入`echo *`就会列出当前文件夹下的所有文件（不包括隐藏文件）。

## 排名前三的捷径

`man bash`中列出了很多快捷方式，找到完整的列表并不难。这份清单都是我经常使用的，并按照使用频率排序。

我不推荐一下子记住所有东西，但我建议你先选择一个，并试着记住并使用它，直到你会无意识地输入它，然后下一个。我也会跳过一些最明显的（如`!!` - 重复上一个命令和`~` - 你的家目录）。

- `!$`

我每天都会使用好多次。它会重复最后一个命令的最后一个参数。使你不必再一次输入同样的参数，它可以节省很多的精力：

```
grep somestring /long/path/to/some/file/or/other.txt
vi !$
```

- `!:1-$`

加点魔法就可以更进一步。它把所有的参数传给前面的命令。所以：

```
grep isthere /long/path/to/some/file/or/other.txt
egrep !:1-$
fgrep !:1-$
```

翻译一下：`!`意思是“看前一个命令”；`:`是一个分隔符；`1`表示“取第一个词”， `-`表示“直到”，`$`表示“最后一个词”。

注意：您可以使用`!*`来实现同样的功能。了解了上述内容，你也可以控制特定连续的参数子集，如：`!:2-3`。

## 启动顺序

bash 运行启动脚本的顺序会导致很多头疼的事情发生。我一直把这张图放在身边（引自[这个](https://blog.flowblok.id.au/2013-02/shell-startup-scripts.html)伟大的页面）：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180207/BmbD401Bg2.png?imageslim)

它展示了 bash 在顶部运行哪些脚本，是基于 bash 运行的环境（用不同的颜色表示）所决定的。

所以，如果你在一个本地（非远程）的，非登录式shell（如，当你从命令行运行bash本身），就好像是在“绿色”通道，这是读取文件的顺序：

```
/etc/bash.bashrc
~/.bashrc
[bash runs, then terminates]
~/.bash_logout
```

这可以节省很多调试的时间。

## getopts（cheapci）

如果你逐渐深入 bash，你最终可能会写一些非常复杂的工具。如果这样，那么使用`getopts`有着巨大的好处。

为了好玩，我曾经写过一个名为`cheapci`的[脚本](https://github.com/ianmiell/cheapci/blob/master/cheapci)来做和 Jenkins 差不多的工作。
\# Jenkins是一个开源软件项目，是基于Java开发的一种持续集成工具，用于监控持续重复的工作，旨在提供一个开放易用的软件平台，使软件的持续集成变成可能。

这个[代码](https://github.com/ianmiell/cheapci/blob/master/cheapci#L70-L96)实现了读取[两个必需的和14个非必需的参数](https://github.com/ianmiell/cheapci/blob/master/cheapci#L33-L51)。要好好地学习这一点，而不是写一堆代码，因为随着项目的发展，你的代码会很迅速地变得非常混乱。

\# 关于`getopts`参见：[使用getopts处理shell中的输入参数](http://www.linuxidc.com/Linux/2016-06/132102.htm)
