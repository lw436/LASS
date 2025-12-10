# 清空环境，加载依赖包
rm(list = ls())
options(stringsAsFactors = F)

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(harmony)

# 设置地址
x = getwd()
setwd(x)

# 加载上一步数据，
load("./step1_data_reading.Rdata")
sce
phe[1:4,]

# 计算线粒体基因比例
mito_genes = rownames(sce)[grep("^MT-", rownames(sce),ignore.case = T)] 
print(mito_genes) #可能是13个线粒体基因，小鼠数据基因名为小写"^mt-"
#sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
sce = PercentageFeatureSet(sce, features = mito_genes, col.name = "percent_mito")
fivenum(sce@meta.data$percent_mito)

# 计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
print(ribo_genes)
sce.all=PercentageFeatureSet(sce.all,  features = ribo_genes, col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)

# 计算红血细胞基因比例
Hb_genes=rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
print(Hb_genes)
sce.all=PercentageFeatureSet(sce.all,  features = Hb_genes,col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
head(sce.all@meta.data)

# 可视化
# "nFeature_RNA"、"nCount_RNA"可视化
feats <- c("nFeature_RNA","nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1
w = length(unique(sce.all$orig.ident))/3+5;w
ggsave(filename="./outcome/step2.vlnplot1.pdf",plot=p1,width = w,height = 5)

# "percent_mito"、"percent_ribo"、"percent_hb"可视化
feats <- c("percent_mito","percent_ribo","percent_hb")
p2 = VlnPlot(sce.all,group.by = "orig.ident",features = feats,pt.size = 0,ncol = 3,same.y.lims=T)+ 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2
w = length(unique(sce.all$orig.ident))/2+5;w
ggsave(filename="./outcome/step2.vlnplot2.pdf",plot=p2,width = w,height = 5)

# counts 和 features 相关性
p3 = FeatureScatter(sce.all,"nCount_RNA","nFeature_RNA",group.by = "orig.ident",pt.size = 0.5)
p3
ggsave(filename="./outcome/step2.scatterplot.pdf",plot=p3)

p = p1+p2+p3
ggsave(filename="./outcome/step2.filter.all.pdf",plot=p,width =21,height = 7)

# 过滤 counts 和 features
selected_c <- WhichCells(sce.all, expression = (nFeature_RNA > 200 & nFeature_RNA < 8000))
selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA$counts > 0) > 3] # 只留下总表达量大于 3 的基因
sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
dim(sce.all) 
dim(sce.all.filt) 

# 过滤线粒体和核糖体
selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 10)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = h < 1 )

sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)

dim(sce.all.filt)
table(sce.all.filt$orig.ident)

# 可视化
# "nFeature_RNA"、"nCount_RNA"可视化
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
w=length(unique(sce.all.filt$orig.ident))/3+5;w 
ggsave(filename="./outcome/step2.vlnplot1.filtered.pdf",plot=p1_filtered,width = w,height = 5)

# "percent_mito"、"percent_ribo"、"percent_hb"可视化
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
w=length(unique(sce.all.filt$orig.ident))/2+5;w 
ggsave(filename="./outcome/step2.vlnplot2.filtered.pdf",plot=p2_filtered,width = w,height = 5) 

p_filtered = p1_filtered+p2_filtered
ggsave(filename="./outcome/step2.vlnplot.filtered.pdf",plot=p_filtered,width = 14,height = 7) 

# 数据标准化
sce.all.filt <- NormalizeData(sce.all.filt, 
                             normalization.method = "LogNormalize",
                             scale.factor = 1e4) 

# 高变基因（目前用不到）
sce.all.filt <- FindVariableFeatures(sce.all.filt)
p4 <- VariableFeaturePlot(sce.all.filt) 
p4
ggsave(filename="./outcome/step2.VariableFeaturePlot.pdf",plot=p4,width = 7,height = 7) 

# 归一化
sce.all.filt <- ScaleData(sce.all.filt)

# PCA 降维
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))

p5 = VizDimLoadings(sce.all.filt, dims = 1:2, reduction = "pca")
p5
ggsave(filename="./outcome/step2.VariableFeaturePlot.pdf",plot=p5,width = 7,height = 7) 

p6 = DimPlot(sce.all.filt, reduction = "pca") + NoLegend()
p6
ggsave(filename="./outcome/step2.PCAplot.pdf",plot=p6,width = 7,height = 7) 

p7 = DimHeatmap(sce.all.filt, dims = 1:12, cells = 500, balanced = TRUE)
p7
ggsave(filename="./outcome/step2.Heapmap.pdf",plot=p7,width = 7,height = 7) 

p8 = ElbowPlot(sce.all.filt, ndims = 50)
p8
ggsave(filename="./outcome/step2.ElbowPlot.pdf",plot=p8,width = 7,height = 7) 

seuratObj = sce.all.filt
names(seuratObj@reductions)

# umap(pca)
seuratObj = RunUMAP(seuratObj,dims = 1:20,reduction = "pca")
DimPlot(seuratObj,reduction = "umap",label=F ) 
p9_pca = DimPlot(seuratObj,reduction = "umap",label=F)
ggsave(filename="./outcome/step2.umap.pca.pdf",plot = p9_pca,width = 7,height = 7)

# tsne(pca)
seuratObj = RunTSNE(seuratObj,dims = 1:20,eduction = "pca")
DimPlot(seuratObj,reduction = "tsne",label=F)
p10_pca = DimPlot(seuratObj,reduction = "tsne",label=F)
ggsave(filename="./outcome/step2.tsne.pca.pdf",plot = p10_pca,width = 7,height = 7)

# Harmony去批次
seuratObj <- RunHarmony(sce.all.filt, "orig.ident")
names(seuratObj@reductions)

# umap(harmony)
seuratObj = RunUMAP(seuratObj,dims = 1:20,reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=F) 
p9_harmony = DimPlot(seuratObj,reduction = "umap",label=F)
ggsave(filename="./outcome/step2.umap.harmony.pdf",plot = p9_harmony)

# tsne(harmony)
seuratObj = RunTSNE(seuratObj,dims = 1:20,reduction = "harmony")
DimPlot(seuratObj,reduction = "tsne",label=F)
p10_harmony = DimPlot(seuratObj,reduction = "tsne",label=F)
ggsave(filename="./outcome/step2.tsne.harmony.pdf",plot = p10_harmony)

p9 = p9_pca + p9_harmony
p9
ggsave(filename="./outcome/step2.umap.pdf",plot = p9)

p10 = p10_pca + p10_harmony
p10
ggsave(filename="./outcome/step2.tsne.pdf",plot = p10)

names(seuratObj@reductions)

# FindNeighbors
sce.all.filt = seuratObj
sce.all.filt = FindNeighbors(sce.all.filt,reduction = "harmony",dims = 1:20) 
sce.all.filt.all=sce.all.filt

# FindNeighbors 设置不同的分辨率，观察分群效果
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all.filt.all=FindClusters(sce.all.filt.all, #graph.name = "CCA_snn", 
                             resolution = res, algorithm = 1)
}
colnames(sce.all.filt.all@meta.data)
apply(sce.all.filt.all@meta.data[,grep("RNA_snn",colnames(sce.all.filt.all@meta.data))],2,table)

p11_dim_low = plot_grid(ncol = 3, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
ggsave(plot = p11_dim_low, filename="./outcome/step2.Dimplot_diff_resolution_low.pdf",width = 21)
  
p11_dim_high = plot_grid(ncol = 3, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
ggsave(plot = p11_dim_high, filename="./outcome/step2.Dimplot_diff_resolution_high.pdf",width = 21)

p12_tree = clustree(sce.all.filt.all@meta.data, prefix = "RNA_snn_res.")
p12_tree
ggsave(plot = p12_tree, filename="./outcome/step2.Tree_diff_resolution.pdf")
table(sce.all.filt.all@active.ident) 

# 选择 0.2 作为分群的标准
p13 = DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2")
p13
ggsave(plot = p13, filename="./outcome/step2.tsne.0.2.pdf")

DimPlot(sce.all.filt.all, reduction = "tsne", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2")

sce.all.int = sce.all.filt.all
save(sce.all.int,phe,file = "./outcome/step2_data_filtering.Rdata")