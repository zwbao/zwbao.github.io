# 16s 分析实战（QIIME2）

主要参考了 QIIME2文档的英文版（ https://docs.qiime2.org/2018.4/ ）以及中文版（ https://forum.qiime2.org/t/qiime2-chinese-manual/838 ）。

> 导入二代测序数据方法参见：https://docs.qiime2.org/2018.6/tutorials/importing/#fastq-manifest-formats

## 分析流程

[pipeline](https://www.processon.com/embed/5b6c116ce4b067df5a0216c0 ':include :type=iframe width=100% height=800px')

## 代码

```
source activate qiime2-2018.4
# import
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Samples_list.csv --output-path paired-end-demux.qza --source-format PairedEndFastqManifestPhred33
# 对拆分样品的结果和质量进行统计
qiime demux summarize --i-data paired-end-demux.qza --o-visualization paired-end-demux.qzv
# 去噪并生成Feature表和代表性序列，这步花了5个小时，时间太长，可选用 Deblur
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza  \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --p-n-threads 0 \
  --o-denoising-stats denoising-stats.qza
# At this stage, you will have artifacts containing the feature table and corresponding feature sequences. You can generate summaries of those as follows.
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
# As well, you can visualize the denoising stats by running:
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
# 多序列比对
qiime alignment mafft \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza
# 移除高变区
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza
# 建树
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza
# 无根树转换为有根树
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
# 查看table.qzv 确定深度（选最小即可）77753
# 计算多样性(包括所有常用的Alpha和Beta多样性方法)，输入有根树、Feature表、样本重采样深度(一般为最小样本数据量，或覆盖绝大多数样品的数据量)
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 77753 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results
  
# 输出结果包括多种多样性结果，文件列表和解释如下：
# beta多样性bray_curtis距离矩阵 bray_curtis_distance_matrix.qza 
# alpha多样性evenness(均匀度，考虑物种和丰度)指数 evenness_vector.qza
# alpha多样性faith_pd(考虑物种间进化关系)指数 faith_pd_vector.qza
# beta多样性jaccard距离矩阵 jaccard_distance_matrix.qza
# alpha多样性observed_otus(OTU数量)指数 observed_otus_vector.qza
# alpha多样性香农熵(考虑物种和丰度)指数 shannon_vector.qza
# beta多样性unweighted_unifrac距离矩阵，不考虑丰度 unweighted_unifrac_distance_matrix.qza
# beta多样性unweighted_unifrac距离矩阵，考虑丰度 weighted_unifrac_distance_matrix.qza

# 统计faith_pd算法Alpha多样性组间差异是否显著，输入多样性值、实验设计，输出统计结果
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

# 统计evenness组间差异是否显著
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

# 网页展示结果，只要是qzv的文件，均可用qiime tools view查看或在线https://view.qiime2.org/查看
qiime tools view evenness-group-significance.qzv
```

Option 2: Deblur

```
# 11m
qiime quality-filter q-score \
 --i-demux paired-end-demux.qza \
 --o-filtered-sequences demux-filtered.qza \
 --o-filter-stats demux-filter-stats.qza
# 44m 比dada2快了许多
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 150 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza
#
qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv
qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv
# 改名
mv rep-seqs-deblur.qza rep-seqs.qza
mv table-deblur.qza table.qza
# 可视化
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
# 多序列比对
qiime alignment mafft \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza
# 移除高变区
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza
# 建树
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza
# 无根树转换为有根树
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
# 查看table.qzv 确定深度（选最小即可）259966
# 计算多样性(包括所有常用的Alpha和Beta多样性方法)，输入有根树、Feature表、样本重采样深度(一般为最小样本数据量，或覆盖绝大多数样品的数据量)
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 259966 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results

# 输出结果包括多种多样性结果，文件列表和解释如下：
# beta多样性bray_curtis距离矩阵 bray_curtis_distance_matrix.qza 
# alpha多样性evenness(均匀度，考虑物种和丰度)指数 evenness_vector.qza
# alpha多样性faith_pd(考虑物种间进化关系)指数 faith_pd_vector.qza
# beta多样性jaccard距离矩阵 jaccard_distance_matrix.qza
# alpha多样性observed_otus(OTU数量)指数 observed_otus_vector.qza
# alpha多样性香农熵(考虑物种和丰度)指数 shannon_vector.qza
# beta多样性unweighted_unifrac距离矩阵，不考虑丰度 unweighted_unifrac_distance_matrix.qza
# beta多样性unweighted_unifrac距离矩阵，考虑丰度 weighted_unifrac_distance_matrix.qza

# 统计faith_pd算法Alpha多样性组间差异是否显著，输入多样性值、实验设计，输出统计结果
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

# 统计evenness组间差异是否显著
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

# 网页展示结果，只要是qzv的文件，均可用qiime tools view查看或在线https://view.qiime2.org/查看，以后不再赘述
qiime tools view evenness-group-significance.qzv
```

## 训练数据

```
# import
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 99_otus.fasta \
  --output-path 99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --source-format HeaderlessTSVTaxonomyFormat \
  --input-path 99_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza
# Extract reference reads
qiime feature-classifier extract-reads \
  --i-sequences 99_otus.qza \
  --p-f-primer ATTACCGCGGCKGCTGG  \
  --p-r-primer CCTACGGGNGGCWGCAG\
  --p-trunc-len 150 \
  --o-reads ref-seqs.qza
# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza
# Test the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

```

报错：

```

Plugin error from feature-classifier:

  No matches found

Debug info has been saved to /tmp/qiime2-q2cli-err-3tgbhs0v.log
```

训练出错，直接用：
- [Silva 119 99% OTUs full-length sequences](https://data.qiime2.org/2018.4/common/silva-119-99-nb-classifier.qza) (MD5: `ea31f08bd07ff510f179fa7e8af8014e`)
- [Silva 119 99% OTUs from 515F/806R region of sequences](https://data.qiime2.org/2018.4/common/silva-119-99-515-806-nb-classifier.qza) (MD5: `aeb2f3144ba8b64b118167d3b7c8c52c`)
- [Greengenes 13_8 99% OTUs full-length sequences](https://data.qiime2.org/2018.4/common/gg-13-8-99-nb-classifier.qza) (MD5: `bb72a9e3f1a4c810dd50bceef3508105`)
- [Greengenes 13_8 99% OTUs from 515F/806R region of sequences](https://data.qiime2.org/2018.4/common/gg-13-8-99-515-806-nb-classifier.qza) (MD5: `919f2a946433b4ff2f9b554d3a99f6ac`)

```
# 40s
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
# OTU表添加count，因为ANCOM不允许有零
qiime composition add-pseudocount \
  --i-table table.qza \
  --o-composition-table comp-table.qza

# 采用ancon，按BodySite分组进行差异统计
qiime composition ancom \
  --i-table comp-table.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column Subject \
  --o-visualization ancom-Subject.qzv

# 查看结果
qiime tools view ancom-BodySite.qzv
```

# 导入数据遇到的坑

以“Fastq manifest” formats 的格式导入数据，明明格式正确却还是报错：

![](_v_images/_1531373990_22184.png)

```
$ cat se-33-manifest
sample-id,absolute-filepath,direction
S029,/data/baozw/Shanxitest_18.7.11/data/S029_zy91-A_1_AHNJWYCCXY_S29_L002_R1_001.fastq.gz,forward
S029,/data/baozw/Shanxitest_18.7.11/data/S029_zy91-A_1_AHNJWYCCXY_S29_L002_R2_001.fastq.gz,reverse
S030,/data/baozw/Shanxitest_18.7.11/data/S030_zy91-A_2_AHNJWYCCXY_S30_L002_R1_001.fastq.gz,forward
S030,/data/baozw/Shanxitest_18.7.11/data/S030_zy91-A_2_AHNJWYCCXY_S30_L002_R2_001.fastq.gz,reverse
```

**问题出在文件编码格式！！！**

> create the manifest as UTF-8 (or ASCII, a subset of UTF-8).

用以下命令查看文件编码：
```
file -i file_name
```

