# 我的Git笔记（三）

至此，我们已经了解了关于`git`的一些基本操作，现在我们将学习如何用`git`进行合作，其实质就是和你的伙伴们共享`commit`。那么如何来共享彼此的`commit`呢？这时候就需要用到远程库，比如著名的`GitHub`。

## 工作流程

让我们先看看用`git`进行合作的工作流程：假如你正在进行着一个项目，同时你也已经在你的本地仓库部署了这个项目，且进行了若干次`commit`，现在你想要把这些`commit`共享给你的同事：

1. 首先你会在服务器上创建一个远程库；
2. 接着，把你的`commit`提交`push`至这个库；

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180625/15Hd9D12id.png?imageslim)

3. 这样你的同事就可以通过`clone`这个库来得到你的项目；
4. 之后，你的同事又对这个项目做了一些改变，并把它们`commit`到她的本地库，然后也同时提交`push`到了远程库；

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180625/fKic5ic5dG.png?imageslim)

5. 你又可以拉取`pull`你的同事提交至远程库的版本。

在多数情况下，当你和你的同事处理同一文件的不同部分时，`Git`会自动找出合并这些更改的最佳方式。但要是你们修改了同一文件的同一个部分时，`Git`会出现错误，它无法协调你们的修改。这时候，我们只能手动解决这些冲突的部分。通常情况下，我们通过事先规划或逐步纳入更改来避免这种混乱局面。 同样，与合作伙伴进行良好的沟通和规划也可以有效防止Git合并冲突。此外，频繁地`push`或`pull`可以使所有成员都保持在最新版本的文件。

## 我们开始吧

### 在GitHub上建立一个远程库

首先请自行注册一个GitHub账号，由于你的本地Git仓库和GitHub仓库之间的传输是通过SSH加密的，所以，需要一点设置：

第1步：创建SSH Key。在用户主目录（Windows下为`C:\Users\Administrator`）下，看看有没有.ssh目录，如果有，再看看这个目录下有没有`id_rsa`和`id_rsa.pub`这两个文件，如果已经有了，可直接跳到下一步。如果没有，打开Shell（Windows下打开Git Bash），创建SSH Key：

```
$ ssh-keygen -t rsa -C "youremail@example.com"
```

你需要把邮件地址换成你自己的邮件地址，然后一路回车，使用默认值即可，由于这个Key也不是用于军事目的，所以也无需设置密码。

如果一切顺利的话，可以在用户主目录里找到`.ssh`目录，里面有`id_rsa`和`id_rsa.pub`两个文件，这两个就是SSH Key的秘钥对，`id_rsa`是私钥，不能泄露出去，`id_rsa.pub`是公钥，可以放心地告诉任何人。

第2步：登陆GitHub，打开“Settings”，“SSH Keys”页面：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180626/924iHg3I3g.png?imageslim)

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180626/2FEGJa82II.png?imageslim)

然后，点“New SSH Key”，填上任意Title，在Key文本框里粘贴`id_rsa.pub`文件的内容：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180626/hj4AJKGg52.png?imageslim)

点“Add Key”，你就应该看到已经添加的Key：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180626/EG2J8C66hd.png?imageslim)

GitHub需要识别出你推送的提交确实是你推送的，而不是别人冒充的，而Git支持SSH协议，所以，GitHub只要知道了你的公钥，就可以确认只有你自己才能推送。

当然，GitHub允许你添加多个Key。假定你有若干电脑，你一会儿在公司提交，一会儿在家里提交，只要把每台电脑的Key都添加到GitHub，就可以在每台电脑上往GitHub推送了。

现在，你可以在GitHub首页点击“New repository”新建一个名为`zmayssnps`的仓库：

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180626/Dh8IKg74l3.png?imageslim)

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180626/CHcbCFKGm1.png?imageslim)

#### 关于GitHub

- GitHub的公共仓库是免费的，任何人都可以看到，所以不要把敏感信息放进去。私人仓库需要额外付费，另一个办法是自己动手，搭一个Git服务器，因为是你自己的Git服务器，所以别人也是看不见的。如果你是学生的话可以申请一下[GitHub的学生大礼包](https://education.github.com/pack)包含了免费的私人仓库。
- 如果你使用GitHub进行合作，那么所有的成员都需要一个GitHub账号。
- 默认情况下，对于你创建的仓库只有你有`push`的权限，所以你需要在GitHub设置里把你的同事也添加进去。

### 使用`Git remote`和远程仓库进行联系

目前，我们这个仓库还是空的，我们可以根据GitHub的提示，从这个仓库克隆出新仓库，也可以吧一个已有的本地仓库与之关联，然后，把本地仓库的内容推送到GitHub仓库。

![mark](http://oo3g995ih.bkt.clouddn.com/blog/180626/350EHFg7gD.png?imageslim)

```
$ git remote add origin git@github.com:username/zmays-snps.git
```

在这行命令在中，需要注意的是我们不仅仅指定了远程仓库的地址`git@github.com:username/zmays-snps.git`，而且还为他指定了一个名称`origin`。现在，试试输入`git remote -v`，可以看到你的本地库与哪些远程库进行了关联。

```
$ git remote -v
origin git@github.com:username/zmays-snps.git (fetch)
origin git@github.com:username/zmays-snps.git (push)
```

如果你需要改变与本地库关联的远程仓库可以先删除再加入：

```
git remote rm origin
git remote add origin git@github.com:username/zmays-snps.git
```

### 使用`git push`提交你的`commit`至远程库

当我们添加了远程库后，就可以开始把我们的`commit`上传到远程仓库了（我们直接在之前创建的`Project`文件夹操作即可）。用Git进行合作的核心就是不断地上传你的工作，让你的同事可以看到并修改，然后再将这些改变传输到本地。

```
$ git push origin master
Counting objects: 15, done.
Delta compression using up to 8 threads.
Compressing objects: 100% (10/10), done.
Writing objects: 100% (15/15), 1.23 KiB | 631.00 KiB/s, done.
Total 15 (delta 2), reused 0 (delta 0)
remote: Resolving deltas: 100% (2/2), done.
To https://github.com/zwbao/zmayssnps.git
 * [new branch]      master -> master
```

 推送成功后，可以在GitHub页面中看到远程库的内容已经和本地的一模一样。

从现在起，只要本地作了提交，就可以通过命令：

```
$ git push origin master
```

把本地`master`分支的最新修改推送至GitHub。

### 使用`git pull`从远程库拉取`commit`

当你将新的提交推送到远程库时，这时候你同事的版本就已经过期了。在继续工作之前，她需要将这些提交拉取到本地。Git上的协作就像是一种来回交换，其中一个人将其最新的版本提交到远程库，其他协作者将这些更改下载到其本地库，进行自己的更改和提交，然后将这些提交推送到远程库供其他人查看和修改。

举个例子，现在，我们将把我们自己的仓库克隆到不同的目录中，模仿不同的成员进行项目合作。首先将远程存储库克隆到名为`zmay-snps-barbara` 的本地目录。 这个目录名反映出这个本地存储库是我们的同事Barbara的库。

```
$ git clone git@github.com:zwbao/zmayssnps.git
Cloning into 'zmayssnps'...
The authenticity of host 'github.com (13.229.188.59)' can't be established.
RSA key fingerprint is SHA256:nThbg6kXUpJWGl7E1IGOCspRomTxdCARLviKw6E5SY8.
Are you sure you want to continue connecting (yes/no)? yes
Warning: Permanently added 'github.com,13.229.188.59' (RSA) to the list of known hosts.
remote: Counting objects: 15, done.
remote: Compressing objects: 100% (8/8), done.
remote: Total 15 (delta 2), reused 15 (delta 2), pack-reused 0
Receiving objects: 100% (15/15), done.
Resolving deltas: 100% (2/2), done.
```

> - SSH警告
> 
> 当你第一次使用Git的`clone`或者`push`命令连接GitHub时，会得到一个警告：
> 
> ```
> The authenticity of host 'github.com (xx.xx.xx.xx)' can't be established.
> RSA key fingerprint is xx.xx.xx.xx.xx.
> Are you sure you want to continue connecting (yes/no)?
> ```
> 
> 这是因为Git使用SSH连接，而SSH连接在第一次验证GitHub服务器的Key时，需要你确认GitHub的Key的指纹信息是否真的来自GitHub的服务器，输入`yes`回车即可。
> 
> Git会输出一个警告，告诉你已经把GitHub的Key添加到本机的一个信任列表里了：
> 
> ```
> Warning: Permanently added 'github.com' (RSA) to the list of known hosts.
> ```
> 
> 这个警告只会出现一次，后面的操作就不会有任何警告了。
> 
> 如果你实在担心有人冒充GitHub服务器，输入`yes`前可以对照[GitHub的RSA Key的指纹信息](https://help.github.com/articles/what-are-github-s-ssh-key-fingerprints/)是否与SSH连接给出的一致。

现在这两个仓库已经有了相同的`commit`，你也可以用`git log`来确认一下。让我们返回`Project`目录进行一些修改，然后提交至远程库：

```
$  echo "Samples expected from sequencing core 2018-06-26" >> README.md
$ git commit -a -m "added information about samples"
warning: LF will be replaced by CRLF in README.md.
The file will have its original line endings in your working directory.
[master 1be82cf] added information about samples
 1 file changed, 1 insertion(+)
$ git push origin master
Counting objects: 3, done.
Delta compression using up to 8 threads.
Compressing objects: 100% (3/3), done.
Writing objects: 100% (3/3), 367 bytes | 367.00 KiB/s, done.
Total 3 (delta 1), reused 0 (delta 0)
remote: Resolving deltas: 100% (1/1), completed with 1 local object.
To https://github.com/zwbao/zmayssnps.git
   77eae85..1be82cf  master -> master
```

现在`zmay-snps-barbara` 目录下的版本需要更新，所以Barbara需要这样做：

```
$ # in zmays-snps-barbara/zmayssnps/
$ git pull origin master
remote: Counting objects: 3, done.
remote: Compressing objects: 100% (2/2), done.
remote: Total 3 (delta 1), reused 3 (delta 1), pack-reused 0
Unpacking objects: 100% (3/3), done.
From github.com:zwbao/zmayssnps
 * branch            master     -> FETCH_HEAD
   77eae85..1be82cf  master     -> origin/master
Updating 77eae85..1be82cf
Fast-forward
 README.md | 1 +
 1 file changed, 1 insertion(+)
```

在`zmay-snps-barbara` 目录下，看看`git log`的结果，可以看到相同的`commit`：

```
$ # in zmays-snps-barbara/zmayssnps/
$ git log --pretty=oneline --abbrev-commit
1be82cf (HEAD -> master, origin/master, origin/HEAD) added information about samples
77eae85 added .gitignore
13770d8 remove test.txt
043a28d add test.txt
630fbbb Add a line
d217683 Add README
```

### 使用Pushing 和 Pulling 来进行合作

Barbara对这个项目的`README.md`文件做了一些修改：

```
$ # in zmays-snps-barbara/zmayssnps/  -- Barbara's version
$ echo "\n\nMaize reference genome version: refgen3" >> README.md
$ git commit -a -m "added reference genome info"
warning: LF will be replaced by CRLF in README.md.
The file will have its original line endings in your working directory.
[master 734e105] added reference genome info
 1 file changed, 1 insertion(+)
$ git push origin master
Counting objects: 3, done.
Delta compression using up to 8 threads.
Compressing objects: 100% (3/3), done.
Writing objects: 100% (3/3), 357 bytes | 357.00 KiB/s, done.
Total 3 (delta 1), reused 0 (delta 0)
remote: Resolving deltas: 100% (1/1), completed with 1 local object.
To github.com:zwbao/zmayssnps.git
   1be82cf..734e105  master -> master
```

现在，我们的版本已经落后了，所以切换到我们的`Project/`目录，把这些更新pull下来：

```
$ # in Project/ -- our version
$ git pull origin master
remote: Counting objects: 3, done.
remote: Compressing objects: 100% (2/2), done.
remote: Total 3 (delta 1), reused 3 (delta 1), pack-reused 0
Unpacking objects: 100% (3/3), done.
From https://github.com/zwbao/zmayssnps
 * branch            master     -> FETCH_HEAD
   1be82cf..734e105  master     -> origin/master
Updating 1be82cf..734e105
Fast-forward
 README.md | 1 +
 1 file changed, 1 insertion(+)
 
$ cat README.md
Git is a version control system.
Git is free software.Git is free software distributed under the GPL.
Project started 2018-06-11
Samples expected from sequencing core 2018-06-26
\n\nMaize reference genome version: refgen3
```

查看最近的两次提交记录：

```
$ git log -n 2
commit 269aa09418b0d47645c5d077369686ff04b16393
Author: Barbara <barbara@barbarasmaize.com>
Date: Sat Sep 28 22:58:55 2013 -0700

added reference genome info

commit 46f0781e9e081c6c9ee08b2d83a8464e9a26ae1f
Author: Vince Buffalo <vsbuffaloAAAAAA@gmail.com>
Date: Tue Sep 24 00:31:31 2013 -0700

added information about samples
```

这就是用Git来合作的记录，可以看到哪些人做了哪些修改。

### 合并冲突

有时候当你`pull`到本地的时候，Git会警告你发生了一个合并错误。这时候，需要你手动输入一些东西来告诉Git如何来处理这个冲突。处理这些错误的方法总是相似的：

1. 首先输入`git status`看看哪些文件冲突了；
2. 打开并编辑这些文件，手动解决冲突；
3. 用`git add`递交这些冲突的文件到暂存区；
4. 用`git status`来检查一下，再用`git commit`提交；之后，你就可以马上`push`到远程仓库，这样你的合作者可以看到你解决了这个冲突，并在更新的版本上继续他们的工作。

> 本文主要参考了《Bioinformatics Data Skills》第五章以及[廖雪峰Git教程](https://www.liaoxuefeng.com/wiki/0013739516305929606dd18361248578c67b8067c8c017b000 )
