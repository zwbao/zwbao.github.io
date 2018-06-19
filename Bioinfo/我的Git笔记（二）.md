# 我的Git笔记（二）

## 使用`git diff`查看文件之间的不同

现在，我们可以用`Git`将文件上传到暂存区并且进行`commit`，也可以用`git status`命令来查看哪些文件被跟踪，哪些文件有更改，哪些文件正在进行下一次提交。除此之外，另一个命令在这个过程中也非常有用，就是`git diff`。如果没有任何参数，`git diff`将比较工作区与暂存区。例如，如果我们新添加一行到`README.md`并运行`git diff`：

```
$  echo "Project started 2018-06-11" >> README.md
$ git diff
diff --git a/README.md b/README.md
index 97440ca..f95be81 100644
--- a/README.md
+++ b/README.md
@@ -1,2 +1,3 @@
 Git is a version control system.
 Git is free software.Git is free software distributed under the GPL.
+Project started 2018-06-11
```

```
--- a/README.md
+++ b/README.md
```

这两行表示了现在是`a`和`b`两个文件在进行比较；`---`表示我们最近一次提交的版本，`+++`则表示已经发生改变的版本

```
 Git is a version control system.
 Git is free software.Git is free software distributed under the GPL.
+Project started 2018-06-11
```

以空格开头的行表示没有被改变的内容，以`+`开头的是被添加的行，`-`开头的是被删除的行（这里我们只添加了一行所以没有这个符号）。

当我们将改变过的文件`git add`后再运行`git diff`就不会有任何输出了，因为`git diff`只是比较工作区与暂存区。

```
$ git add README.md
$ git diff # shows nothing
```

要是我们想比较暂存区的文件和最近一次提交的文件，可以使用` git diff --staged`。

```
$  git diff --staged
diff --git a/README.md b/README.md
index 97440ca..f95be81 100644
--- a/README.md
+++ b/README.md
@@ -1,2 +1,3 @@
 Git is a version control system.
 Git is free software.Git is free software distributed under the GPL.
+Project started 2018-06-11
```

## 使用`git log`查看你的`commit`记录

Git中的`commit`就是在某个时间点的快照，并且每个提交（除第一个之外）指向其父提交; 这个提交链就是一组连接的快照，展示了你的库是如何演变的。

![commit](http://oo3g995ih.bkt.clouddn.com/blog/180611/6Fj0b23j1j.png?imageslim)

我们可以使用`git log`来可视化`commit`记录：

```
$ git log
commit 630fbbb693aaeb46a980a3a7d983cf6dafe2adaf (HEAD -> master)
Author: zwbao <shinningbzw@foxmail.com>
Date:   Mon Jun 11 19:09:18 2018 +0800

    Add a line

commit d217683121237715e1f96f177c0cc11f01857b6a
Author: zwbao <shinningbzw@foxmail.com>
Date:   Mon Jun 11 17:39:05 2018 +0800

    Add README

```

随着我们的提交，这条命令的输出也会随之变长。如果你想要看看一个很长的`git log`输出不妨看看我们之前克隆的`seqtk`，在那个文件夹试试这条命令。

## 使用`git rm`删除文件

我们先添加一个新文件`test.txt`到Git并且提交：

```
$ git add test.txt
$ git commit -m "add test.txt"
[master 043a28d] add test.txt
 2 files changed, 1 insertion(+)
 create mode 100644 test.txt
```

一般情况下，当你不需要`test.txt`时，可以直接删除它或者使用`rm`命令：

```
$ rm test.txt
```

这时，Git也会发现你删除了这个文件，输入`git status`：

```
$ git status
On branch master
Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

        deleted:    test.txt

no changes added to commit (use "git add" and/or "git commit -a")
```

当你确实要删除这个文件时，就用`git rm`删掉，并进行`git commit`：

```
$ git rm test.txt
rm 'test.txt'

$ git commit -m "remove test.txt"
[master 13770d8] remove test.txt
 1 file changed, 0 insertions(+), 0 deletions(-)
 delete mode 100644 test.txt

```

现在，这个文件就已经完全从版本库中删除了。

## 使用`.gitignore`告诉Git哪些文件需要忽略

你也许已经发现`git status`会列出当前还没有被跟踪的文件，随着文件数目的增多，这个列表也会越来越长。但这些文件中，有些我们确实不想跟踪比如一些测序数据，因为他们太大了（当你提交这些大文件后，当别人克隆你的仓库时就不得不也下载如此大的数据了，之后我们会讲解如何管理这些数据，现在就先把他们忽略吧）。

如果我们希望忽略`data/seqs/`文件夹下的所有`FATSQ`文件，那么就下创建并编辑一个`.gitignore`文件，输入：

```
data/seqs/*.fastq
```

现在，输入`git status`看看：

```
$ git status
On branch master
Untracked files:
  (use "git add <file>..." to include in what will be committed)

        .gitignore

nothing added to commit but untracked files present (use "git add" to track)
```

接着，我们继续添加和提交这个文件：

```
$ git add .gitignore

$ git commit -m "added .gitignore"
[master 77eae85] added .gitignore
 1 file changed, 1 insertion(+)
 create mode 100644 .gitignore
```

### 那么哪些文件是我们应该告诉`.gitignore`将他们忽略的呢？

- 大文件
这些文件应该被忽略，并通过其他方式进行管理。 大文件会减慢创建，推送和提交的速度。当别人克隆你的仓库时，这可能会导致相当大的麻烦。

- 中间文件
生物信息学项目往往充满了中间文件。 例如，你将read回帖到基因组上时，则会创建`SAM`或`BAM`文件。 即使这些不是大文件，这些也应该被忽略。因为如果一个数据文件可以通过重新运行命令（或脚本）轻松地重新生成，我们通常只需存储它的创建方式。 

- 文本编辑器的临时文件
Emacs和Vim等文本编辑器有时会在你的目录中创建一些临时文件。这些可能是像`textfle.txt〜`或`＃textfle.txt＃`。在Git中跟踪这些内容也是没有意义的，所以，这些文件应该添加到`.gitignore`中。幸运的是，`.gitignore`也支持通配符，所以这些文件可以用`*〜`和`\＃* \＃`来忽略。

- 临时的代码文件
一些编程语言的编译器（例如`Python`）通常也会产生一些临时文件（例如`overlap.pyc`），需要将他们忽略。

## 使用`git reset`撤销暂存区的修改

如果你不小心把一个错误的修改`git add`到暂存区了，你可以使用`git reset`来取消它。 例如：

```
$ echo "TODO: ask sequencing center about adapters" >> README.md
$ git add README.md
$ git status
# On branch master
# Changes to be committed:
# (use "git reset HEAD <file>..." to unstage)
# #
new file: README.md
#
```

用`git status`查看后，我们可以发现我们对于`README.md`的改变已经包含在下一次的提交中，为了撤销这次修改，根据`git status`的提示输入`git reset HEAD README.md`：

```
$ git reset HEAD README.md
$ git status
# On branch master
# Changes not staged for commit:
# (use "git add <file>..." to update what will be committed)
# (use "git checkout -- <file>..." to discard changes in working
directory)
# #
modified: README.md
#
```


在`Git`中，我们用`HEAD`来表示当前版本，也就是最新的提交，上一个版本就是`HEAD^`，上上一个版本就是`HEAD^^`，当然往上100个版本写100个`^`比较容易数不过来，所以写成`HEAD~100`。

现在我们虽然已经撤回了这次的`git add`，但是`README.md`文件中已经添加了`TODO: ask sequencing center about adapters`这一行，倘若我们要丢弃工作区的修改呢？

## 使用`git checkout`丢弃工作区的修改

用`git checkout -- file`命令可以丢弃工作区的修改：

```
$ git checkout -- README.md
$ git status
On branch master
nothing to commit, working tree clean
```
