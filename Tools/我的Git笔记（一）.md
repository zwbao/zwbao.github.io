# 我的Git笔记（一）

## Git 是什么

Git是目前世界上最先进的分布式版本控制系统（没有之一）。版本控制（Revision control）是一种在开发的过程中用于管理我们对文件、目录或工程等内容的修改历史，方便查看更改历史记录，备份以便恢复以前的版本的软件工程技术。在处理生物信息数据的过程中，随着内容的增多，我们面临着如何管理不同的文件版本的挑战，同时，我们也希望可以有一种更好的方式来团队协作。对于普通用户来说例如Dropbox，Google Drive或坚果云可能是个不错的选择，它可以备份数据并且允许我们恢复我们所需要的某一个版本的文件。然而，这种专门的文件版本管理系统并不能很好地适应复杂的生物信息学项目，将整个生物信息学项目放入共享目录也不太现实，因为它可能包含千兆或者更多的数据，这些数据太大也无法在整个网络中共享。

碰巧，软件工程师们在协作的过程中也遇到了相同的问题，Linus编写了Git来管理Linux，这是一个包含了数千个协作者同时更改和处理文件的大型代码库。Git非常适合用于项目版本控制和协作工作。诚然，Git刚开始学习起来可能非常棘手，但只要勤加练习，在苦苦挣扎之后，你会发现它会远远超出你的期望。

![分布式版本控制](http://oo3g995ih.bkt.clouddn.com/blog/180702/HK6Il2d22m.png?imageslim)

## 为什么Git对于管理生物信息数据是必不可少的

### Git允许你创建项目的快照

使用版本控制系统，你可以在开发中的特定点处创建当前项目的快照。如果出现任何错误，你可以回退到项目的过去时间点（一个时间点即为一次`commit`）并恢复文件，所以在生物信息学的工作中，这样的安全措施非常有用。例如，假设你正在进行SNP数据的分析，你发现你的SNPs中有14％落在染色体的某一段编码区中，接着你也在你的论文中引用了这个百分比；两个月后，你早已忘记了关于这个分析的一些细节，但当你重新运行这些分析代码时，数据竟然变成了26％！如果你通过Git来跟踪项目的开发，那么当结果发生变化时，你可以很容易地查看项目的整个历史记录。`commit`允许你轻松再现和回滚到过去的分析版本。查看每一次`commit`，甚至比较任何两次`commit`之间的差异。

### Git帮助你跟踪代码的重要变化

随着新功能的添加或错误的修复，大多数软件或脚本都会随着时间而变化。Git在帮助你跟踪代码的变化方面非常有帮助。假设一位聪明的生物信息学家，编写了一个脚本，用于从`reads`中删去质量较差的区域。这位生物信息学家随后将其分发给他的所有实验室成员。一个月后，生物信息学家发现有一个bug会导致某些结果不正确。Git便可以很容易地跟踪软件的更改并下载最新的版本。此外，比如GitHub和Bitbucket这类在Web上托管的Git仓库也使得代码得以实现共享和协作。

### Git帮助人们保持软件的组织性和可用性

无序的代码会影响自己和其他实验室成员的不便，代码丢失更是会导致实验结果无法再现，并影响之后的实验。Git有助于保持工作的连续性和项目历史的完整记录。将整个项目集中到一个仓库中可以保持其组织性。Git存储了每一次提交的变更，因此即使主开发人员离开，项目的整个历史记录都依旧可用。 由于Git能够回滚到过去的版本，因此修改项目的风险也更小，从而更容易构建现有的工作。

## 安装和使用Git

### 安装Git

最早Git是在Linux上开发的，很长一段时间内，Git也只能在Linux和Unix系统上跑。不过，慢慢地有人把它移植到了Windows上。现在，Git可以在Linux、Unix、Mac和Windows这几大平台上正常运行了。

首先你可以在命令行输入`git`，看看系统有没有安装Git：

```
$ git
The program 'git' is currently not installed. You can install it by typing:
sudo apt-get install git
```

如果你使用OSX，可以使用Homebrew安装（比如`brew install git`）；而在Debian或Ubuntu Linux上则可以使用`apt-get`来安装。

- 更多方式请参见：https://git-scm.com/downloads

### Git理论基础

#### 工作区域

Git本地有三个工作区域：工作目录（Working Directory）、暂存区(Stage/Index)、资源库(Repository或Git Directory)。如果在加上远程的git仓库(Remote Directory)就可以分为四个工作区域。文件在这四个区域之间的转换关系如下：

![工作区域](http://oo3g995ih.bkt.clouddn.com/blog/180702/AjDA374Ff4.png?imageslim)

- Workspace：工作区，就是你平时存放项目代码的地方；
- Index / Stage：暂存区，用于临时存放你的改动，保存即将提交到文件列表信息；
- Repository：仓库区（或本地仓库），就是安全存放数据的位置，这里面有你提交到所有版本的数据。其中HEAD指向最新放入仓库的版本；
- Remote：远程仓库，托管代码的服务器，可以简单的认为是你项目组中的一台电脑用于远程数据交换。

#### 工作流程

git的工作流程一般是这样的：
1. 在工作目录中添加、修改文件；
2. 将需要进行版本管理的文件放入暂存区域；
3. 将暂存区域的文件提交到git仓库。

因此，git管理的文件有三种状态：已修改（modified）,已暂存（staged）,已提交(committed)。

![工作流程](http://oo3g995ih.bkt.clouddn.com/blog/180702/CL2l0h0EIl.png?imageslim)

### Git入门：创建仓库，跟踪文件和提交修改

#### 告诉Git你是谁

当你安装Git后要做的第一件事就是设置你的用户名称和e-mail地址。这是非常重要的，因为每次Git提交都会使用该信息。

```
$ git config --global user.name "zwbao"  # 名称
$ git config --global user.email "shinningbzw@foxmail.com" # 邮箱
```

我们可以通过`git <subcommand>`的子命令与Git进行交互。 Git有许多的子命令，但在日常工作中你只需要一些。另一个十分有用的Git设置是启用的终端颜色， 这样可以更直观地显示一些更改（例如，红色表示删除，绿色表示新的文件或修改）。

```
$ git config --global color.ui true
```

#### 创建仓库

要开始使用Git，我们首先需要将目录初始化为Git仓库。Git仓库是一个受版本控制的目录。它包含你当前的工作文件和项目在特定时间点的快照。在版本控制中，这些快照被称为`commit`，而使用Git无非就是创建和操作这些`commit`。

创建仓库有两种主要的方式：从一个现有目录初始化一个，或者是克隆别人的仓库。这次我们在`Project/`目录下新建一个仓库：

```
$ git init
Initialized empty Git repository in C:/Users/zwbao/Desktop/Project/.git/
```

`git init`命令使当前目录下多了一个`.git`隐藏目录（你可以用`ls -a`来查看它），这个目录是Git来跟踪管理版本库的，没事千万不要手动修改这个目录里面的文件，不然改乱了，就把Git仓库给破坏了。

另一种创建仓库的方法是通过克隆现有的库。你可以从任何位置克隆一个仓库：您的文件系统，本地网络或网络上的其他位置。现在，人们使用GitHub和Bitbucket等仓库托管服务，所以从网络上克隆Git仓库是也是很常见的。

让开始我们练习如何从GitHub克隆仓库吧。在这个例子中，我们将克隆Heng Li大神的GitHub中的Seqtk代码。Seqtk是SEQuence ToolKit的缩写，包含了一套用于处理FASTQ和FASTA文件的工具。首先，让我们访问GitHub仓库（https://github.com/lh3/seqtk），在这个页面上，右侧是复制URL的地方。

![](http://oo3g995ih.bkt.clouddn.com/_1527761422_385.png)

现在，让我们切换到`Project/`之外的其他目录，并从该目录运行：

```
$ git clone https://github.com/lh3/seqtk.git
Cloning into 'seqtk'...
remote: Counting objects: 319, done.
remote: Total 319 (delta 0), reused 0 (delta 0), pack-reused 319
Receiving objects: 100% (319/319), 141.97 KiB | 147.00 KiB/s, done.
Resolving deltas: 100% (183/183), done.
```

`git clone`命令将`seqtk`克隆到你的本地目录中。

#### 用Git跟踪文件

尽管您已经将`Project/`作为仓库初始化，但Git不会自动开始跟踪此目录中的每个文件。你需要使用`git add`命令告诉Git将哪些文件进行跟踪。这也是十分合理的，因为对于生物信息学数据而言，在一个项目中包含了许多我们不想跟踪的文件，比如一些大的数据文件，中间结果或任何可以通过命令轻松重新生成的文件（之后我们可以通过`.gitignore`来告诉Git自动忽略这些文件）。

首先创建一个`readme.txt`文件，内容如下：

```
Git is a version control system.
Git is free software.
```

>需要注意的是，所有的版本控制系统，其实只能跟踪文本文件的改动，比如TXT文件，网页，所有的程序代码等等，Git也不例外。版本控制系统可以告诉你每次的改动，比如在第5行加了一个单词“Linux”，在第8行删了一个单词“Windows”。而图片、视频这些二进制文件，虽然也能由版本控制系统管理，但没法跟踪文件的变化，只能把二进制文件每次改动串起来，也就是只知道图片从100KB改成了120KB，但到底改了啥，版本控制系统不知道，也没法知道。
不幸的是，Microsoft的Word格式是二进制格式，因此，版本控制系统是没法跟踪Word文件的改动的。如果要真正使用版本控制系统，就要以纯文本方式编写文件。

在跟踪文件之前，让我们使用命令`git status`查看Git仓库中的文件状态（如果你在其他地方，请切换到`Project/`目录）：

```
$ git status
On branch master

No commits yet

Untracked files:
  (use "git add <file>..." to include in what will be committed)

        readme.txt

nothing added to commit but untracked files present (use "git add" to track)
```

`git status`命令告诉我们：

- 我们正在`master`分支，也是Git的默认分支。分支允许我们在不同版本的项目之间转换，现在我们先在`master`分支上工作，在之后的章节中我们会学习更多的分支。
- 我们现在有一些`Untracked files`，现在是包括了在目录中的所有文件（其实就一个`readme.txt`），因为啥都还没提交嘛。

`git status`将会成为你最常用的命令之一，它会告诉你哪些文件被改动了，哪些文件已经准备好提交，哪些文件没有被跟踪。

接下来让我们用`git add`告诉Git将`readme.txt`进行跟踪：

```
$ git add readme.txt
```

现在我们已经跟踪了`readme.txt`，用`git status`看看：

```
$ git status
On branch master

No commits yet

Changes to be committed:
  (use "git rm --cached <file>..." to unstage)

        new file:   readme.txt

```

现在，我们试着给`readme.txt`文件加一行：

```
$ echo "Git is free software distributed under the GPL." >>readme.txt
```

接下来使用`git commit -m "your commit message"`来提交他们，`-m`选项后可以添加一些注释。

- 提交：
```
$ git commit -m "Add README"
[master (root-commit) d05a690] Add README
 1 file changed, 2 insertions(+)
 create mode 100644 readme.txt
```

提交后再看看状态：

```
$ git status
On branch master
Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

        modified:   readme.txt

no changes added to commit (use "git add" and/or "git commit -a")
```

可以看到第二次的修改并没有被提交，这是因为使用`git add`后，在工作区的第一次修改被放入暂存区，准备提交，但是，在工作区的第二次修改并没有放入暂存区，所以`git commit`只负责把暂存区的修改提交了，也就是第一次的修改被提交了，第二次的修改不会被提交。

要是想提交第二次的修改，可以继续`git add`再`git commit`，也可以别着急提交第一次修改，先`git add`第二次修改，再`git commit`，就相当于把两次修改合并后一块提交了。

## 小结

本文主要介绍了：

- Git在处理生物信息学项目中的用途；
- Git的一些概念；
- 简单介绍了Git的流程。

在接下来的教程中小编会介绍Git更多实用的命令和功能，另外后台回复`cheatsheet`有惊喜~

本文主要参考了《Bioinformatics Data Skills》第五章以及廖雪峰Git教程（https://www.liaoxuefeng.com/wiki/0013739516305929606dd18361248578c67b8067c8c017b000）
