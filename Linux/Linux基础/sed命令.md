# sed命令
## sed

**sed**是一种行编辑器，逐行读取，它是文本处理中非常中的工具，能够完美的配合正则表达式使用，功能不同凡响。处理时，把当前处理的行存储在临时缓冲区中，称为“模式空间”（pattern space），接着用sed命令处理缓冲区中的内容，处理完成后，把缓冲区的内容送往屏幕。接着处理下一行，这样不断重复，直到文件末尾。文件内容并没有 改变，除非你使用重定向存储输出。Sed主要用来自动编辑一个或多个文件；简化对文件的反复操作；编写转换程序等。

```shell
sed [options] 'AddressCommand' file ...
```

### options：

```shell
-n：静默模式，不再默认显示模式空间中的内容 # 使用 p Command 会显示两次符合的行，用这个选项后只显示一次，因为 p 命令是显示符合条件的行，而原来在模式空间的也会一并输出，故显示两次
-i：直接修改原文件
-e SCRIPT -e SCRIPT：可以同时执行多个脚本
-f /PATH/TO/SED_SCRIPT
	sed -f /path/to/scripts  file
	以选项中指定的script文件来处理输入的文本文件
-r：表示使用扩展正则表达式
```

### Address：

1. StartLine,EndLine

```shell
#从第1行行到100行，用逗号隔开
比如1,100
$：最后一行
#倒数第二行
$-1
```
2. /RegExp/

```shell
#以root开头的行
/^root/
```
3. /pattern1/,/pattern2/

```shell
第一次被pattern1匹配到的行开始，至第一次被pattern2匹配到的行结束，这中间的所有行
```
4. LineNumber

  指定的行
5. StartLine, +N

```shell
从指定行开始，向后的N行，共 N+1 行；
```

### Command：

```shell
d：删除符合条件的行；
p：显示符合条件的行；
a \string：在指定的行后面追加新行，内容为 string，要以 \ 开头
	\n：可以用于换行
i \string：在指定的行前面添加新行，内容为string
r FILE：将指定的文件的内容添加至符合条件的行处
w FILE：将 Address 指定范围内的行另存至指定的文件中;
s/pattern/string/修饰符：查找并替换，默认只替换每行中第一次被模式匹配到的字符串
	加修饰符：
		g：全局替换
		i：忽略字符大小写
s/// 也可以写为 s###, s@@@ （这时候就匹配 / 就不用转义啦）
	后向引用：
	\(\), \1, \2
	#例：like-->liker
	     love-->lover
	sed 's#\(l..e\)#\1r#' sed.txt
	&：引用模式匹配整个串
	sed 's/l..e/&r/g' sed.txt
```

sed练习：

```shell
# 删除/etc/grub.conf文件中行首的空白符；
sed -r 's@^[[:spapce:]]+@@g' /etc/grub.conf
# 删除/etc/inittab文件中的空白行；
sed '/^$/d' /etc/inittab
# 删除/etc/inittab文件中开头的#号;
sed 's@^#@@g' /etc/inittab
# 删除某文件中开头的#号及后面的空白字符，但要求#号后面必须有空白字符;
sed -r 's@^#[[:space:]]+@@g' /etc/inittab
# 删除某文件中以空白字符后面跟#类的行中的开头的空白字符及#
sed -r 's@^[[:space:]]+#@@g' /etc/inittab
```

练习：
传递三个参数给脚本，第一个为整数，第二个为算术运算符，第三个为整数，将计算结果显示出来，要求保留两位精度。形如：

```shell
# 保留两位精度：
echo "scale=2;111/22" | bc
bc <<< "scale=2;111/22"
```
