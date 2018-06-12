# bash脚本编程
## 一些基础知识：
### 高级语言可分为：

- 静态语言：编译型语言

  强类型(变量)语言
  需要事先转换成可执行格式
  C、C++、JAVA、C#

- 动态语言：解释型语言（ on the fly）

  通常为弱类型语言
  边解释边执行
  PHP、SHELL、python、perl

- 面向过程：Shell, C

- 面向对象：JAVA, Python, perl, C++

### 变量类型：事先确定数据的存储格式和长度

- 之所以要区分变量类型因为：

    - 当10作为字符时：每个字符8位，共16位
    - 当10作为数字时：1010，4位，因为计算机的最小存储单元是字节，故最终占用8位

- 逻辑：1+1>2
- 逻辑运算：与、或、非、异或
1：真
0：假

- 与：
1 & 0  = 0
0 & 1 = 0
0 & 0 = 0
1 & 1 = 1

- 或

- 非：
! 真 = 假
! 假 = 真

- 短路逻辑运算：
    - 与：只要有一个为假，结果一定为假
    - 或：只要有一个为真，结果一定为真

- shell：是弱类型编程语言
    - 强类型：变量在使用前，必须事先声明，甚至还需要初始化（数值为0，字符串为空`NULL`）；
    - 弱类型：变量用时声明，甚至不区分类型（默认为字符串）；
- 变量赋值：VAR_NAME=VALUE

## bash的变量：

bash变量的类型：

```shell
本地变量(局部变量)：
本地变量：
set VARNAME=VALUE：作用域为整个bash进程，set 可略；
引用变量：${VARNAME}，不影响变量名时可以略去括号。
局部变量：
local VARNAME=VALUE：作用域为当前代码段；

环境变量：作用域为**当前shell进程及其子进程**，脚本在执行时会启动一个子shell进程（命令行中启动的脚本会继承当前shell环境变量；系统自动执行的脚本(非命令行启动)就需要自我定义需要各环境变量）
	export VARNAME=VALUE
	VARNAME=VALUE
	export VARNAME
	“导出”
位置变量：
	$1, $2, ... #表示命令的第 1 个参数，第 2 个参数...
	shift：前一个位置变量shift后可重复使用，shift(shift 1) 命令每执行一次，变量的个数($#)减一（之前的$1变量被销毁,之后的$2就变成了$1），而变量值提前一位。同理，shift n后，前n位参数都会被销毁。
	#输入1 2 后显示 1 2
	echo $1
	shift
	echo $1
	shift
特殊变量：
	$?：上一条命令的退出状态码
		程序执行，可能有两类返回值：
			1.程序执行结果
			2.程序状态返回代码（0-255）
				0：正确执行
				1-255：错误执行，1，2，127系统预留，有特殊意义
	$#：参数的个数
	$*：参数列表
	$@：参数列表
```
输出重定向的常用技巧：
```
/dev/null：软件设备，bit bucket，数据黑洞，一旦输出到这里则不会输出
【例】检验是否有以 student 为名的用户
id student &> /dev/null （因为此输出结果不重要）
echo $?
返回 0，即存在
```

撤消变量：注意没有`$` ，撤销的是变量本身而不是变量值
unset VARNAME

查看当shell中变量：不带任何选项和参数，包括环境变量和本地变量

```shell
set
```

查看当前shell中的环境变量：

```shell
printenv
env
export
```

在已有变量上添加变量：（如添加环境变量）

```shell
#在环境变量后面添加变量
export PATH=$PATH:/usr/local/apache/bin
#在环境变量前面添加变量
export PATH=/usr/local/apache/bin:$PATH
```

## 脚本：命令的堆砌，按实际需要，结合命令流程控制机制实现的源程序

shebang：魔数

```shell
#!/bin/bash
#注释行，不执行
```

### 条件判断：

​	如果用户不存在
​		添加用户，给密码并显示添加成功；
​	否则
​		显示如果已经没在，没有添加；

bash中如何实现条件判断？看`echo $?` 的值
条件测试类型：

- 整数测试
- 字符测试
- 文件测试

条件测试的表达式：
	[ expression ]：表达式两段必须要有空格
	[[ expression ]]
	test expression

整数比较:
```
-eq：测试两个整数是否相等；比如 $A -eq $B；
-ne：测试两个整数是否不等；不等，为真；相等，为假；
-gt：测试一个数是否大于另一个数；大于，为真；否则，为假；
-lt：测试一个数是否小于另一个数；小于，为真；否则，为假；
-ge：大于或等于
-le：小于或等于
```

命令的间逻辑关系：
```shell
逻辑与： &&
	第一个条件为假时，第二条件不用再判断，最终结果已经有；
	第一个条件为真时，第二条件必须得判断；
逻辑或： ||
```

```shell
#如果用户user6不存在，就添加用户user6
! id user6 && useradd user6
id user6 || useradd user6
#如果用户存在，就显示用户已存在；否则，就添加此用户；
id user1 && echo "user1 exists." || useradd user1
#如果用户不存在，就添加；否则，显示其已经存在；
! id user1 && useradd user1 || echo "user1 exists."
#如果用户不存在，添加并且给密码；否则，显示其已经存在；
! id user1 && useradd user1 && echo "user1" | passwd --stdin user1	|| echo "user1 exists."
#如果/etc/inittab文件的行数大于100，就显示好大的文件；
[ `wc -l /etc/inittab | cut -d' ' -f1` -gt 100 ] && echo "Large file."
```

取变量名称：
1. 只能包含字母、数字和下划线，并且不能数字开头；
2. 不应该跟系统中已有的环境变量重名；
3. 最好做到见名知义；



练习，写一个脚本，完成以下要求：

1. 添加3个用户user1, user2, user3；但要先判断用户是否存在，不存在而后再添加；最后显示当前系统上共有多少个用户；

   ```shell
   #!/bin/bash
   ! id user1 &> /dev/null && useradd user1 && echo "user1" |passed --stdin user1 &> /dev/null || echo "user1 exists."
   ! id user2 &> /dev/null && useradd user2 && echo "user2" |passed --stdin user2 &> /dev/null || echo "user2 exists."
   ! id user3 &> /dev/null && useradd user3 && echo "user3" |passed --stdin user3 &> /dev/null || echo "user3 exists."
   USERS=`wc -l /etc/passwd | cut -d：-f1`
   echo "$USERS users."
   ```

2. 添加完成后，显示一共添加了几个用户；当然，不能包括因为事先存在而没有添加的；

练习，写一个脚本，完成以下要求：
给定一个用户：
1. 如果其UID为0，就显示此为管理员；

2. 否则，就显示其为普通用户；

   ```shell
   #!/bin/bash
   NAME=user1
   USERID=`id -u $NAME`
   [ $USERID -eq 0] && echo "Admin" ||echo "Common user"
   ```

条件判断，控制结构：

- 单分支if的结构

  ```shell
  if 判断条件; then
  	statement1
  	statement2
  	...
  fi
  ```

- 双分支的if结构

  ```shell
  if 判断条件; then
  	statement1
  	statement2
  	...
  else
  	statement3
  	statement4
  fi

  #上个例题
  #!/bin/bash
  NAME=user1
  if [ `id -u $NAME` eq 0 ]; then
  	echo "Admin"
  else
  	echo "Common user"

  #上上个例题
  #!/bin/bash
  NAME=user1
  if id $NAME &> /dev/null; then
  	echo "$NAME exists."
  else
  	useradd $NAME
  	echo $NAME | passed --stdin $NAME &> /dev/null
  	echo "Add $NAME finished."
  fi
  ```

- 多分支的if语句

  ```shell
  if 判断条件1; then
    statement1
    ...
  elif 判断条件2; then
    statement2
    ...
  elif 判断条件3; then
    statement3
    ...
  else
    statement4
    ...
  fi
  ```

练习：写一个脚本
判断当前系统上是否有用户的默认shell为bash；如果有，就显示有多少个这类用户；否则，就显示没有这类用户；

```shell
#!/bin/bash
grep "bash$" /etc/passwd &> /dev/null
RETVAL=$?
if [ $RETVAL -eq 0 ]; then
	USERS=`grep "bash$" /etc/passwd |wc -l`
	echo "The shells of $USERS is bash."
else
	echo "No such user."
```

练习：写一个脚本
判断当前系统上是否有用户的默认shell为bash； 如果有，就显示其中一个的用户名；否则，就显示没有这类用户；

```shell
#!/bin/bash
grep "bash$" /etc/passwd &> /dev/null
RETVAL=$?
if [ $RETVAL -eq 0 ]; then
	AUSERS=`grep "bash$" /etc/passwd | head -1 | cut -d：-f1`
	echo "$AUSER is one of such users."
else
	echo "No such user."
```

练习：写一个脚本
给定一个文件，比如/etc/inittab；判断这个文件中是否有空白行；如果有，则显示其空白行数；否则，显示没有空白行。

```shell
#!/bin/bash

FILE=/etc/inittab
if [ ! -e $FILE ]; then
  echo "No $FILE."
  exit 8
fi

if grep "^$" $FILE &> /dev/null; then
  echo "Total blank lines：`grep "^$" $FILE | wc -l`."
else
  echo "No blank line."
fi
```

练习：写一个脚本
给定一个用户，判断其UID与GID是否一样；如果一样，就显示此用户为“good guy”；否则，就显示此用户为“bad guy”。

```shell
#!/bin/bash
USERNAME=user1
USERID=`id -u $USERNAME`
GROUPID=`id -g $USERNAME`
if [ $USERID -eq $GROUPID ]; then
  echo "Good guy."
else
  echo "Bad guy."
fi

#进一步要求：不使用id命令获得其id号；
#!/bin/bash
USERNAME=user1
if ! grep "^$USERNAME\>" /etc/passwd &> /dev/null; then
  echo "No such user：$USERNAME."
  exit 1
fi

USERID=`grep "^$USERNAME\>" /etc/passwd | cut -d：-f3`
GROUPID=`grep "^$USERNAME\>" /etc/passwd | cut -d：-f4`
if [ $USERID -eq $GROUPID ]; then
  echo "Good guy."
else
  echo "Bad guy."
fi
```


练习：写一个脚本
给定一个用户，获取其密码警告期限；而后判断用户最近一次修改密码时间距今天是否已经小于警告期限；如果小于，则显示“Warning”；否则，就显示“OK”。

提示：计算方法`bc`，最长使用期限减去已经使用的天数即为剩余使用期限；

```shell
# !/bin/bash
W=`grep "student" /etc/shadow | cut -d：-f6`
S=`date +%s`
T=`expr $S/86400`
L=`grep "^student" /etc/shadow | cut -d：-f5`
N=`grep "^student" /etc/shadow | cut -d：-f3`
SY=$[$L-$[$T-$N]]

if [ $SY -lt $W ]; then
  echo 'Warning'
else
  echo 'OK'
fi
```

练习：写一个脚本
判定命令历史中历史命令的总条目是否大于1000；如果大于，则显示“Some command will gone.”；否则显示“OK”。

```shell
#!/bin/bash
# $HISTORY=`history | wc -l` 不可以这么写，错误！！因为 history 命令看起来有一千多个但是 wc 只显示1000
# $history =`history | tail -1 | cut -d' ' -f1` 也不可以，错误！！tail 结果第一列是空格
$history =`history | tail -1 | cut -d' ' -f2`
if [ $HISTORY -gt 1000 ]; then
	echo "Some command will gone."
else
	echo "OK"
fi
```

### shell中如何进行算术运算：

bash默认会圆整：丢弃小数点后的所有内容

```shell
# 1. let 算术运算表达式
let C=$A+$B
# 2. $[算术运算表达式]
C=$[$A+$B]
# 3. $((算术运算表达式))
C=$(($A+$B))
# 4. expr 算术运算表达式，表达式中各操作数及运算符之间要有空格，而且要使用命令引用
C=`expr $A + $B`
```

### exit

**exit命令**同于退出shell，并返回给定值。在shell脚本中可以终止当前脚本执行。执行exit可使shell以指定的状态值退出。若不设置状态值参数，则shell以预设值退出。状态值0代表执行成功，其他值代表执行失败。如果脚本没有明确定义退出状态码，那么，最后执行的一条命令的退出码即为脚本的退出状态码（所以需要自己指定一个状态码）。

```shell
#给定一个用户，判断其UID与GID是否一样；如果一样，就显示此用户为“good guy”；否则，就显示此用户为“bad guy”。
#!/bin/bash
USERNAME=user1
if ! grep "^$USERNAME\>" /etc/passwd &> /dev/null; then
  echo "No such user：$USERNAME."
  exit 1 # 执行 echo $? 显示 1
fi

USERID=`grep "^$USERNAME\>" /etc/passwd | cut -d：-f3`
GROUPID=`grep "^$USERNAME\>" /etc/passwd | cut -d：-f4`
if [ $USERID -eq $GROUPID ]; then
  echo "Good guy."
else
  echo "Bad guy."
fi
```

### 测试方法：（不要少了expression两边的空格！！）

1. [ expression ]
2. [[ expression ]]
3. test expression

bash中常用的条件测试有三种：

- 整数测试：

  ```shell
  -eq：测试两个整数是否相等；比如 $A -eq $B
  -ne：测试两个整数是否不等；不等，为真；相等，为假；
  -gt：测试一个数是否大于另一个数；大于，为真；否则，为假；
  -lt：测试一个数是否小于另一个数；小于，为真；否则，为假；
  -ge：大于或等于
  -le：小于或等于
  ```

- 文件测试：

  ```shell
  -e FILE：测试文件是否存在；
  -f FILE：测试文件是否为普通文件；
  -d FILE：测试指定路径是否为目录；
  -r FILE：测试当前用户对指定文件是否有读取权限；
  -w FILE：测试当前用户对指定文件是否有写入权限；
  -x FILE：测试当前用户对指定文件是否有执行权限；
  ```

  练习（完善之前的练习）：写一个脚本
  给定一个文件，比如/etc/inittab
  判断这个文件中是否有空白行；
  如果有，则显示其空白行数；否则，显示没有空白行。

  ```shell
  #!/bin/bash
  FILE=/etc/inittab
  if [ ! -e $FILE ]; then
    echo "No $FILE."
    exit 8
  fi

  if grep "^$" $FILE &> /dev/null; then
    echo "Total blank lines：`grep "^$" $FILE | wc -l`."
  else
    echo "No blank line."
  fi
  ```


- 字符测试：

  ```shell
  == / =：测试是否相等，相等为真，不等为假；#等号两端一定要加空格，否则会识别为赋值
  !=：测试是否不等，不等为真，等为假
  >
  <
  -n string：测试指定字符串是否为空，空则真，不空则假
  -z string：测试指定字符串是否不空，不空为真，空则为假
  ```

- 组合测试条件

  ```shell
  -a: 与关系
  -o: 或关系
  !： 非关系
  ```

练习：写一个脚本
传递一个参数(单字符就行)给脚本，如参数为q、Q、quit或Quit，就退出脚本；否则，就显示用户的参数；

```shell
#!/bin/bash
#
if [ $1 = 'q' ]; then
	echo "Quiting..."
	exit 1
elif [ $1 = 'Q' ]; then
	echo "Quiting..."
	exit 2
elif [ $1 = 'quit' ]; then
	echo "Quiting..."
	exit 3
elif [ $1 = 'Quit' ]; then
	echo "Quiting..."
	exit 4
else
	echo $1
fi
```
### 测试脚本是否有语法错误：

```shell
bash -n 脚本
bash -x 脚本：单步执行
```

bash变量的类型：
```shell
本地变量(局部变量)：作用域为整个bash进程；
局部变量：作用域为当前代码段；
环境变量：作用域为当前shell进程及其子进程；
位置变量：
	$1, $2, ... #表示命令的第 1 个参数，第 2 个参数...
	shift：前一个位置变量shift后可重复使用，shift(shift 1) 命令每执行一次，变量的个数($#)减一（之前的$1变量被销毁,之后的$2就变成了$1），而变量值提前一位。同理，shift n后，前n位参数都会被销毁。
	#输入1 2 后显示 1 2
	echo $1
	shift
	echo $1
	shift
特殊变量：
	$?：上一条命令的退出状态码
	$#：参数的个数
	$*：参数列表
	$@：参数列表
```

练习：写一脚本
能接受一个参数(文件路径)
判定：此参数如果是一个存在的文件，就显示“OK.”；否则就显示"No such file."

```shell
#!/bin/bash
#
if [ $# -lt 1 ]; then
	echo "Usage：./filetest3.sh ARG1 [ARG2 ...]"
	exit 7
fi

if [ -e $1 ]; then
	echo "OK."
else
	echo "No such file."
fi
```

练习：写一个脚本
给脚本传递两个参数(整数)；显示此两者之和，之乘积；

```shell
# !/bin/bash
#
if [ $# -lt 2 ]; then
  echo "Usage：cacl.sh ARG1 ARG2"
  exit 8
fi
echo "The sum is：$[$1+$2]."
echo "The prod is：$[$1*$2]."
```

## for循环：

```shell
for 变量 in 列表; do
  循环体
done

遍历完成之后，退出；

如何生成列表：
{1..100}
# `seq [起始数 [步进长度]] 结束数` []表示可略，要用命令替换
seq 1 2 10
1 3 5 7 9
# 声明 SUM 是整型
declare -i SUM=0
	   -x # 声明是环境变量
	   
#计算从 1 加到 100
#!/bin/bash
#
declare -i SUM=0
for I in {1..100}; do
	let SUM=$[$SUM+$I]
done
echo "The sum is $SUM."
```

写一个脚本：

依次向/etc/passwd中的每个用户问好

```shell
LINES=`wc -l /etc/passwd | cut -d' ' -f1`
for I in `seq 1 $LILNES`; do echo "Hello, `head -n $I /etc/passwd | tail -1 |cut -d: -f1`"; done
```
## while循环：适用于循环次数未知的场景，要有退出条件

语法：

```shell
while CONDITION; do
  statement
  ...
done
```

练习：计算100以内所有正整数的和

```shell
#!/bin/bash
declare -i I=1
declare -i SUM=0

while [ $I -le 100 ]; do
  let SUM+=$I
  let I++
done

echo $SUM
```

练习：转换用户输入的字符为大写，除了quit:

```shell
#!/bin/bash
#
read -p "Input something: " STRING

while [ $STRING != 'quit' ]; do
  echo $STRING | tr 'a-z' 'A-Z'
  read -p "Input something: " STRING
done
```

练习：每隔5秒查看hadoop用户是否登录，如果登录，显示其登录；否则，显示当前时间，并说明hadoop尚未登录：

```shell
#!/bin/bash
#
who | grep "hadoop" &> /dev/null
RETVAL=$?

while [ $RETVAL -ne 0 ]; do
  echo "`date`, hadoop is not log." 
  sleep 5
  who | grep "hadoop" &> /dev/null
  RETVAL=$?
done

echo "hadoop is logged in."
```
## basename：

**basename命令**用于打印目录或者文件的基本名称。basename和[dirname](http://man.linuxde.net/dirname)命令通常用于shell脚本中的命令替换来指定和指定的输入文件名称有所差异的输出文件名称。

```shell
$0: 执行脚本时的脚本路径及名称
basename $0
```

## 生成随机数

RANDOM: 0-32768

```shell
echo $RANDOM
```

随机数生成器：熵池
/dev/random：取空了就会停下来，阻塞用户进程，等产生更多随机数，但它更安全
/dev/urandom：取空了，会通过软件的方式模拟产生随机数，不会阻塞用户进程

写一个脚本，利用RANDOM生成10个随机数，并找出其中的最大值，和最小值；

```shell
#!/bin/bash
#
declare -i MAX=0
declare -i MIN=0

for I in {1..10}; do
  MYRAND=$RANDOM
  [ $I -eq 1 ] && MIN=$MYRAND
  if [ $I -le 9 ]; then  # 前面九个随机数后加逗号，最后一个没有逗号；echo -n 表示不换行输出
	echo -n "$MYRAND,"
  else
	echo "$MYRAND"
  fi
  [ $MYRAND -gt $MAX ] && MAX=$MYRAND
  [ $MYRAND -lt $MIN ] && MIN=$MYRAND
done
echo $MAX, $MIN
```

## case语句：
选择结构，比多分支`if`更加清晰易懂。

```shell
case SWITCH in 
value1)
  statement
  ...
  ;;
value2)
  statement
  ...
  ;;
*)
  statement
  ...
  ;;
esac
```

写个脚本：

只接受参数start,stop,restart,status其中之一

```shell
#!/bin/bash
#
case $1 in
'start')
	echo "start server..." ;;
'stop')
	echo "stop server..." ;;
'restart')
	echo "restarting server..." ;;
'status')
	echo "Runing..." ;;
*)
	echo "`basename $0`(start|stop|restart|status)"
esac
```
## read:

**read命令**从键盘读取变量的值，通常用在shell脚本中与用户进行交互的场合。该命令可以一次读取多个变量的值，变量和输入的值都需要使用空格隔开。在read命令后面，如果没有指定变量名，读取的数据将被自动赋值给特定的变量REPLY

```shell
-p：给出提示；
-t：指定读取值时等待的时间（秒）。
```