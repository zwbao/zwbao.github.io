# Linux基础

Linux基础-基本概念和技巧
1.多用户操作系统：同一系统上，同时多用户使用，通过权限隔离
2.home目录：用户默认目录，存放了用户环境配置相关的配置文件
3.文件：Linux一切皆文件，不同类型决定不同用途，与文件后缀无关，只是标识
4.文件属性：ls -l，类型，权限，链接数，用户，用户组，日期，大小，名称
5.文件基本权限：可读，可写，可执行，用户，用户组，其他人
6.绝对路径与相对路径：.. . / - ~
7.自动补全：TAB
8.软连接：利用ln 建立，相当于快捷方式，数据同时处理时可以用它以节约磁盘空间
9.历史命令：上下键，Ctrl+R，history，！+历史序号即可执行相应命令

Linux学习思想：去可视化

常用命令
1.文件目录操作：ls(list), cd(change dictionary), pwd, mkdir, rm, mv(最后一个参数即是移动的位置), cp（不能直接复制目录，加-r 可以）, touch（新建一个空文件，用于更改最后修改时间）, ln
2.系统信息：df（-h 看硬盘参数）, du（-h估计占用空间，即文件大小，-sh可看文件夹大小，ll无法做到）, free（-h查看内存使用情况）, ipconfig（查看网络配置，看服务器IP）, netstat（查看网络进程，系统上多少个进程与网络相关）, top（预览系统运行情况，使用频率高）
3.用户权限：chown（修改权限）, chgrp（改变用户组）, chmod（改变权限，用的多，chmod 750 yy 给自己所有权限，同组可读可执行，其他人无权限）
4.文本操作：grep（筛选，匹配文本，man 看文档！）, sed, awk, paste（以列的行使，把两个文件组合 ）, cat, diff, wc, vi, head, tail, less（按页显示，回车下一行，F翻页，-N 每列加序号，less -S(一行行显示，可以用左右键看)N（显示行号）按G会跳到最后一行）, more（会显示百分比 ）, zcat（直接查看被压缩的文件，zless也有相同功能）, 重定向，管道(uniq 可归一去除重复，uniq -c 输出不重复的数量)

5.进程/作业管理：ps（aux列出所有当前在运行的进程，更多地是与grep结合来筛选进程，如：ps aux | grep “python”）, top, kill（+PID 来结束进程）, jobs, nohup（防止连接断开后程序挂掉，如：nohup tophat &（&符号是指将程序放到后台跑））, bg（放到后台，比如bg %1，一般不会挂掉）, fg（把后台的作业放到前台来，比如 fg %1把刚刚放到后台的作业放到前台来，用jobs查看后台编号）, screen（打开一个窗口，不会停掉，比如在实验室打开screen，关掉后在宿舍使用 screen -a 打开窗口，继续干活）, ctrl+Z（放到后台，暂停）, ctrl+C 
6.其他：find（查找文件，如:find . -name xxx/”samtoo*”可使用通配符；还可以根据类型找：find . type d（d代表目录，f代表普通文件）；其他参数等）, which（查找命令位置，which python）, tar, gzip（不能对目录直接打包，用tar 打包：tar czf xxx.tar.gz xxx/ （z参数对应gz格式，j参数对应bz2格式）；解压：tar xf xxx.tar.gz（xf参数会自动识别用什么解压））, man（看帮助文档）, scp（在两个服务器之间传数据：scp xxx.tar.gz username@ip）, ssh（登录远程系统：ssh username@ip）, history, who（查看当前多少人登录，还可使用last命令查看近一段时间内登录的人，IP，时间）, wget（下载文件）, dos2unix（Windows的文件传输到Linux一般要经过此处理，有些编码改变）

技巧：从别处复制命令时，最好开头加一个#行命令。

Linux文本操作实例
![mark](http://oo3g995ih.bkt.clouddn.com/blog/180122/8hGbf8hIk1.png?imageslim)