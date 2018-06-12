# Linux操作系统基础
## 一些概念：
![mark](http://oo3g995ih.bkt.clouddn.com/blog/180122/13cEBaa0lA.png?imageslim)

### Shell：
即人机交互接口。将用户的行为翻译为计算机（内核）能理解的命令。用户与系统交互必须要通过shell。任何与shell相关的程序，只要shell关闭了也就关闭了。
#### Shell分为：
1. 图形用户界面（GUI，graphic user interface）：Linux 的GUI称为X-window（x表示一种图形显示协议），包含Gnome（使用C开发），KDE（使用C++开发），XFace（简介，轻量级桌面系统）三种。而在Windows系统中只有一种。
2. 命令行接口（界面=接口）（CLI，command line interface）：包括bash（最广泛，开源），zsh，ksh，tcsh等。
- 命令提示符（prompt）：
```shell
#    root
$    普通用户
```
### 命令
命令通过shell传递给内核，并由内核判断该程序是否有执行权限，以及是否能执行，从什么时候开始执行。（任何一个程序想要执行，必须要有执行入口）
#### 命令格式：
1. 命令（command）:
Linux命令类型：
- 1.内部命令（shell 内置）：基本的实现一些管理功能的命令。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180122/iAkbe4lbeJ.png?imageslim)

- 2.外部命令：文件系统中的某个路径下有一个与命令名称相应的可执行文件。
2. 选项（options，用以修改命令的执行方式的选项）：分为短选项（-character，多个选项可以组合）和长选项（-word，不能组合，必须分开）
3. 参数（arguments，命令的作用对象）：多个参数之间用空格隔开。
### Linux的基本原则
1.由目的单一的小程序组成；组合小程序完成复杂任务
2.一切皆文件
3.尽量避免获取用户接口
4.配置文件保存为纯文本格式
### 登录：
#### 虚拟终端（terminal）：
Ctrl+Alt+F1-F6 切换
#### 切换用户：
```shell
su（switch user）
su 用户名 （半切换）
exit（退出）
su -l 用户名 （完全切换）
passwd （修改当前用户密码，连输两次即可）
```
##### 设置密码要符合密码复杂性规则：
1. 使用4种类型（小写字母，大写字母，数字，符号）字符中的至少三种
2. 足够长，大于7位
3. 使用随机字符串
4. 定期更换
5. 循环周期足够大（不要使用最近使用的密码）