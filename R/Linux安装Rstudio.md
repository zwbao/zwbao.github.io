# 在 Linux 安装 Rstudio

1. 使用`conda`创建`R`环境：
```R
conda create -n r35 R=3.5
source activate r35
```
2. 安装必要的`R`包：
```R
conda install r-essentials
```
3. 安装`Rstudio`:
```R
conda install -c r rstudio
```
4. 使用`Rstudio`:
```bash
rstudio
```

如果使用 mobaxterm 登陆服务器，Rstudio 会直接弹出窗口，和 Windows 的没区别。

![](http://ww1.sinaimg.cn/mw690/c5d7b0ebly1g1249pqaj5j21ig0sin3z.jpg)

## 在 Jupyter Notebook 中使用 R

需要安装以下`R`包：
```R
install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'devtools', 'uuid', 'digest'))
devtools::install_github('IRkernel/IRkernel')
# 只在当前用户下安装
IRkernel::installspec()
```
等待执行完毕，打开`jupyter`就可以新建`R`的 notebook 了。