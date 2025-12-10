# 清空环境，加载依赖包
rm(list = ls())
options(stringsAsFactors = F) 

.libPaths(c("/home/zhangy/anaconda3/envs/hdWGCNA/lib/R/library", .libPaths()))

library(Seurat)
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)

# 设置地址

setwd("/home/jiangxijie/UVM/step1/") 
dir = "./data/GSE139829_RAW"

# 批量读取样本
fs = list.files(dir,'^GSM')
samples = fs[grepl("barcodes",fs)]
samples = str_split(samples,"_",simplify = T)[,1]
samples = str_replace(samples,"barcodes.tsv.gz",'')

# 给每个样本创建文件夹
lapply(samples,function(x){
  #x = "GSM4147091"
  y = fs[grepl(x,fs)]
  dir = paste0("./data/GSE139829_RAW/",x)
  dir.create(dir)
  
  file.rename(paste0("./data/GSE139829_RAW/",y[1]),file.path(dir,"barcodes.tsv.gz"))
  file.rename(paste0("./data/GSE139829_RAW/",y[2]),file.path(dir,"features.tsv.gz"))
  file.rename(paste0("./data/GSE139829_RAW/",y[3]),file.path(dir,"matrix.mtx.gz"))
})

# 返回创建的 Seurat 对象，将其存储在sceList中。
samples = list.files(dir)
sceList = lapply(samples,function(pro){ 
    print(pro)  
    tmp = Read10X(file.path(dir,pro)) 
    if(length(tmp)==2){
    ct = tmp[[1]] 
    }else{ct = tmp}
    sce =CreateSeuratObject(counts =  ct ,
                            project =  pro  ,
                            min.cells = 5,
                            min.features = 300 )
    return(sce)
}) 
View(sceList)
do.call(rbind,lapply(sceList, dim))# 观察每个样本的细胞数量和基因数量（左为基因数，右为细胞数）

# 将每个样本 merge 在一起，变成一个 seurat 对象（sce.all）
sce.all = merge(x = sceList[[1]],
                y = sceList[-1],
                add.cell.ids = samples)
sce.all #看一下总细胞数，layers，基因数等信息

# 用JoinLayers函数对layers进行合并(多个 layer 对象变为一个 layer 对象)
sce.all <- JoinLayers(sce.all)
sce.all #看看新的sce对象

table(sce.all$orig.ident) #看看每个样本的细胞量
as.data.frame(sce.all@assays$RNA$counts[1:4, 1:3]) #看看表达矩阵
head(sce.all@meta.data, 4) #meta 信息，细胞过滤的依据（orig.ident样品名、nCount_RNA每个细胞的总counts数目、nFeature_RNA每个细胞中检测到的基因数量）

# 添加分组信息
phe = sce.all@meta.data
phe$group = str_split(rownames(phe),'[_]',simplify = T)[,1] 
#phe$sample = ifelse(phe$group %in% c("GSM4147091","GSM4147092","GSM4147100"),"metastatic","primary")
phe$sample = "primary"
phe[1:4,]
table(phe$sample)

save(sce.all,phe,file = "./outcome/step1_data_reading.Rdata")