# Linux常用命令
## ls（list，列出，列表）：
ls命令用列出指定路径下的所有文件。
### 一些概念：
- 目录也是一种文件，即路径映射文件。
- 路径（为实现层次化文件管理的一种机制）：从指定起始点到目的地所经过的位置。
- 文件系统（file system）

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180122/KhJe6ldLmH.png?imageslim)

从根往下找的即为绝对路径；相对于当前所处位置的路径即相对路径。

### 选项：
```shell
-l 
显示结果：
    第一列有十位：
        第一位为文件类型：
        -：普通文件（f，file）
        d：目录文件
        b：块设备文件（block）
        c：字符设备文件（character）
        p：命令管道（pipe）
        s：套接字文件（socket）
        l：符号链接文件（symbolic link file）
    后九位为文件权限（mode）：每3位一组，每一组：rwx（读，写，执行）没有则“-”
    之后列分别为：
        文件硬链接次数
        文件的属性（owner）
        文件的属组（group）
        文件大小（size）：默认单位是字节
        时间戳（timestamp）：最近一次被修改的时间。
            可分为：
                访问时间（access）；
                   修改时间（modify）：改变内容；
                   改变时间（change）：改变文件的属性数据（即元数据，metadata）
        文件名
-h：做单位的转换
-a：显示以 . 开头的隐藏文件，包括 .（当前目录） ，..（父目录）
-A：显示以 . 开头的隐藏文件，不包括 .（当前目录） ，..（父目录）
-d：显示目录自身属性
-i：显示 index node（iNode）
-r：逆序显示文件（本来按照的是字母顺序）
-R：递归（recursive）显示，即子目录内容也显示
```
### pwd（printing working dictionary）：
pwd命令以绝对路径的方式显示用户当前工作目录。命令将当前目录的全路径名称（从根目录）写入标准输出。全部目录使用`/`分隔。第一个`/`表示根目录，最后一个目录是当前目录。执行pwd命令可立刻得知您目前所在的工作目录的绝对路径名称。
### cd：
cd命令用来切换工作目录至dirname。 其中dirName表示法可为绝对路径或相对路径。若目录名称省略，则变换至使用者的home directory(也就是刚login时所在的目录)。另外，~也表示为home directory的意思，.则是表示目前所在的目录，..则表示目前目录位置的上一层目录。

```bash
cd    进入用户主目录；
cd ~  进入用户主目录；
cd -  返回进入此目录之前所在的目录。
```

### type:
type命令用来显示指定命令的类型，判断给出的指令是内部指令还是外部指令。
#### 命令类型：
- alias：别名。
- keyword：关键字，Shell保留字。
- function：函数，Shell函数。
- builtin：内建命令，Shell内建命令。
- file：文件，磁盘文件，外部命令。
- unfound：没有找到。
### date:
date命令是显示或设置系统时间与日期。
### 如何获得命令的使用帮助
#### help：
help命令用于显示shell内部命令的帮助信息。help命令只能显示shell内部的命令帮助信息。而对于外部命令的帮助信息只能使用man或者info命令查看。
#### man：
man命令是Linux下的帮助指令，通过man指令可以查看Linux中的指令帮助、配置文件帮助和编程帮助等信息。
- man内容：
    - \[ ]：表示可选
    - |：表示多选
    - <>：表示必选
    - …：表示可以使用多次
    - { }：表示分组
    - NAME：命令名称及功能简要说明
    - SYNOPSIS：用法说明，包括可用的选项
    - DESCRIPTION：命令功能的详尽说明，可能包括每一个选项的意义
    - OPTIONS：说明每一个选项的意义
    - FILES：此命令相关的配置文件
    - BUGS
    - EXAMPLES：使用示例
    - SEE ALSO：另外参照
- 翻屏：
    - 向后翻一屏：SPACE
    - 向前翻一屏：b
    - 向后翻一行：ENTER
    - 向前翻一行：k
- 查找：不分大小写
    - /KEYWORD：向后找
    - ？KEYWORD：向后找
    - n：下一个
    - N：前一个
- q：退出
- 分章节：man 2（显示第几章的内容） read
#### whatis：
whatis命令是用于查询一个命令执行什么功能，并将查询结果打印到终端上。
```
[root@localhost ~]# whatis ls
ls                   (1)  - list directory contents
ls                   (1p)  - list directory contents
[root@localhost ~]# whatis cp
cp                   (1)  - copy files and directories
cp                   (1p)  - copy files
[root@localhost ~]# whatis chown
chown                (1)  - change file owner and group
chown                (1p)  - change the file ownership
chown                (2)  - change ownership of a file
chown                (3p)  - change owner and group of a file
```
- 中间列数字表示的含义：
1. 用户命令（所有用户都可使用的命令）（/bin, /usr/bin, /usr/local/bin）
2. 系统调用
3. 库调用
4. 特殊文件（特殊文件）
5. 文件格式（配置文件的语法）
6. 游戏
7. 杂项（miscellaneous）
8. 管理命令（/sbin, /usr/sbin, /usr/local/sbin）：只有管理员有权限，能修改系统级别的配置
#### 命令的帮助文档：/usr/share/doc
#### Google
## which：
which命令用于查找并显示给定命令的绝对路径，环境变量PATH中保存了查找命令时需要遍历的目录。which指令会在环境变量$PATH设置的目录里查找符合条件的文件。也就是说，使用which命令，就可以看到某个系统命令是否存在，以及执行的到底是哪一个位置的命令。
## cal：
cal命令用于显示当前日历，或者指定日期的日历。
## echo：
echo命令用于在shell中打印shell变量的值，或者直接输出指定的字符串。linux的echo命令，在shell编程中极为常用, 在终端下打印变量value的时候也是常常用到的，因此有必要了解下echo的用法echo命令的功能是在显示器上显示一段文字，一般起到一个提示的作用。
```
# 使用-e选项时，若字符串中出现以下字符，则特别加以处理，而不会将它当成一般文字输出：
\a 发出警告声；
\b 删除前一个字符；
\c 最后不加上换行符号；
\f 换行但光标仍旧停留在原来的位置；
\n 换行且光标移至行首；
\r 光标移至行首，但不换行；
\t 插入tab；
\v 与\f相同；
\\ 插入\字符；
\nnn 插入nnn（八进制）所代表的ASCII字符；
```
## printf：
printf命令格式化并输出结果到标准输出。
