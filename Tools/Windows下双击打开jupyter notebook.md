# Windows下双击打开 jupyter notebook

1. 打开`cmd`，使用`pip`下载并安装`nbopen`：

```
pip install nbopen
python -m nbopen.install_win
```

2. 安装后需要关联`.ipynb`和`jupyter notebook`在命令行输入：

```
assoc .whl=jupyter& ftype jupyter=cmd.exe /c jupyter-notebook "%1"
```

缺点：这种方法只能是在已经打开了一个 jupyter notebook 后，再通过双击打开**相同目录**下的文件。

> 原文参见： http://axil.github.io/how-to-open-ipynb-file-with-one-doubleclick-on-windows.html 
