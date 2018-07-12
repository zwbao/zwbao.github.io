# 用Atom优雅地写Python

Atom是由 GitHub 的程序员们打造的称为“属于21世纪”的代码编辑器。它开源免费跨平台（支持 Windows、Mac、Linux 三大桌面平台），并且整合 GIT 并提供类似 SublimeText 的包管理功能，人们可以非常方便地安装和管理各种插件，并将 Atom 打造成真正适合自己的开发工具。

作为一个现代的代码编辑器，Atom 支持各种编程语言的代码高亮(HTML / CSS / Javascript / PHP / Python / C / C++ / Objective C / Java / JSON / Perl / CoffeeScript / Go / Sass / YAML / Markdown 等等)、 与大多数其他编辑器相比，Atom的语言支持已经算是覆盖非常全面了。另外，它的代码补全功能（也叫Snippets） 也非常好用，你只需输入几个字符即可展开成各种常用代码，可以极大提高编程效率。

Atom 编辑器可以和 GIT 完美结合，所有对代码、文本的修改都能体现在编辑器的界面上。比如在文件内新写的代码会在左边标记为绿色，删除的标记为红色，修改的标记为黄色。在左边的目录导航也能方便的看到文件改动：有改动的文件其文件名和所在文件夹名都会被标记为高亮显示。编辑器底部会显示当前所在分支和对文件的修改行数统计，对于 GIT 用户来说非常方便。

接下来，我将带大家一步步把Atom打造成Python的IDE~

## 安装Atom

下载安装Atom：https://atom.io/

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180505/FfEHIGihBB.png?imageslim)

## 在Atom中愉快地运行Python

首先要介绍的Package是**Hydrogen**，一个能在 Atom 中模拟 Jupyter 编辑方式的插件，有了它，我们就可以摆脱浏览器，回到IDE里愉快的借助Jupyter写代码了！

![hydrogen](http://oxnc5ug9u.bkt.clouddn.com/pic/171026/BidI31mDIl.gif)

安装`Hydrogen`后还需要**安装 Jupyter**，小编建议新手直接安装下载安装相应版本的**Anaconda**（https://www.anaconda.com/download/），Anaconda中已包含了Jupyter，它会自动帮你安装配置。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180505/EFFaJhKLF8.png?imageslim)

当然你也可以直接通过`pip`来安装Jupyter：

- 如果你安装的是Python3（推荐）：
```
python3 -m pip install --upgrade pip
python3 -m pip install jupyter
```

- 如果你安装的是Python2：
```
python -m pip install --upgrade pip
python -m pip install jupyter
```

安装完 Jupyter后，就可以选中要运行的代码，使用快捷键`ctrl `+ `Enter`来执行命令了。

![hydrogen](https://pic3.zhimg.com/v2-659bcc57132003e3b2395427d8305703_b.gif)

## Kite
kite是一款Python开发必备神器，集代码**自动补全**，**帮助文档**，**示例代码模板**等功能于一身。小编为大家从YouTube搬运回来了Kite的官方介绍视频，大家可以自己感受下这款工具的强大。

<div align=center>
<iframe width="560" height="315" src="https://www.youtube.com/embed/bF50YPyUKTQ" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
</div>

## 一些常用的插件

### Minimap

让你了解当前屏幕所处相对位置。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180505/d9Cei2mJiE.png?imageslim)

### regex-railroad-diagram

正则表达式图形化插件，只要把鼠标放上去直接自动图形显示。

![regex](https://pic1.zhimg.com/80/f5394efec12c3e3173816ed77fa381df_hd.jpg)

### file-icons

为不同类型的文件添加一个漂亮的小图标，颜控必备。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180505/Dd6kL5KFHi.png?imageslim)

> 得益于atom强大的维护团队，其插件系统发展迅速，还好多好用的插件等有着我们去发掘，通过简单的安装就能使Atom怀十八般武艺，无所不能~