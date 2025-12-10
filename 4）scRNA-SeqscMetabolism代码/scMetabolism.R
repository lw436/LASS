# 清空环境，加载依赖包
rm(list = ls())
options(stringsAsFactors = F,future.globals.maxSize = 170 * 1024^3)

.libPaths(c(
  "/home/zhangy/anaconda3/envs/single/lib/R/library",
  "/home/zhangy/anaconda3/envs/hdWGCNA/lib/R/library",
  "/home/luwei/R/x86_64-pc-linux-gnu-library/4.4",
  "/usr/local/lib/R/site-library",
  "/usr/lib/R/site-library",
  "/usr/lib/R/library",
  "/home/zhangy/R/x86_64-conda-linux-gnu-library/4.4"
))
library(AUCell)
library(Seurat)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(dplyr)
library(gridExtra)
library(stringr)
library(future)
library(tidyverse)
library(ComplexHeatmap)
library(wesanderson)

availableCores() #12 #查看几个核可用
plan(multisession, workers=36) #观察htop的CPU，发现中后期确实使用了多核
nbrOfWorkers() #4 当前可用的核有多少个

# 设置地址
x = getwd()
setwd(x)
set.seed(123)

# 加载上一步数据，
sce = readRDS("/home/jiangxijie/UVM/step4/step4.tumor.new.rds")
sce

result <- sce@meta.data %>%
  group_by(celltype) %>%
  summarise(
    Proportion = mean(Top_LA_Cells == "top25%", na.rm = TRUE),
    Percentage = round(Proportion * 100, 1)
  ) %>%
  arrange(desc(Percentage))  # 按百分比降序排序

print(result)

sce@meta.data$group = factor(ifelse(sce$celltype %in% c("cluster 7","cluster 8","cluster 2"),"HighScore",
                     ifelse(sce$celltype %in% c("cluster 3","cluster 0","cluster 1","cluster 9"),"LowScore","Medium")),
                      levels = c("LowScore","Medium","HighScore"))

Idents(sce) <- factor(sce@meta.data$group,
                      levels = c("LowScore","Medium","HighScore"))

sce[["RNA"]] <- as(sce[["RNA"]], Class = "Assay")

countexp.Seurat <- sc.metabolism.Seurat(obj = sce,  #Seuratde单细胞object
                                      method = "AUCell", 
                                      imputation = F, 
                                      ncores = 36, 
                                      metabolism.type = "KEGG")


#提取score结果
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
#将score中barcode的点转为下划线
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))  #基因ID不规范会报错,下划线替换-
#确定细胞barcode椅子
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
#TRUE
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change))

#可以直接使用Seurat的相关函数
pathways = rownames(countexp.Seurat@assays$METABOLISM$score)

pdf('/home/luwei/UVM_output_new/new_KEGG_scMetabolism.PDF', width = 8, height = 20)

DotPlot.metabolism(countexp.Seurat, 
                   pathway = pathways, 
                   phenotype = "group",  #更改phenotype 参数
                   norm = "y")+
scale_x_discrete(limits = c("LowScore", "Medium", "HighScore")) 

dev.off()

countexp.Seurat@meta.data$group = factor(countexp.Seurat@meta.data$group,levels = c("LowScore","Medium","HighScore"))

#展示特定通路
input.pathway <- c(
  "Riboflavin metabolism",
  "Pyruvate metabolism",
  "Glycolysis / Gluconeogenesis",
  "Fatty acid elongation",
  "Arginine and proline metabolism",
  "Arachidonic acid metabolism",
  "Retinol metabolism",
  "Lysine degradation"
)



pdf('/home/luwei/UVM_output_new/KEGG_scMetabolism.PDF', width = 4.8, height = 10)
DotPlot.metabolism(countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "group",  #更改phenotype 参数
                   norm = "y")+
scale_x_discrete(limits = c("LowScore", "Medium", "HighScore")) 
dev.off()

sce@meta.data$group = factor(ifelse(sce$celltype %in% c("cluster 7","cluster 8","cluster 2"),"HighScore",
                     ifelse(sce$celltype %in% c("cluster 3","cluster 0","cluster 1","cluster 9"),"LowScore","Medium")),
                      levels = c("LowScore","Medium","HighScore"))

Idents(sce) <- factor(sce@meta.data$group,
                      levels = c("LowScore","Medium","HighScore"))

sce[["RNA"]] <- as(sce[["RNA"]], Class = "Assay")

countexp.Seurat <- sc.metabolism.Seurat(obj = sce,  #Seuratde单细胞object
                                      method = "AUCell", 
                                      imputation = F, 
                                      ncores = 36, 
                                      metabolism.type = "REACTOME")


#提取score结果
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
#将score中barcode的点转为下划线
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))  #基因ID不规范会报错,下划线替换-
#确定细胞barcode椅子
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
#TRUE
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )

#可以直接使用Seurat的相关函数
pathways = rownames(countexp.Seurat@assays$METABOLISM$score)

pdf('/home/luwei/UVM_output_new/new_RACTOME_scMetabolism.PDF', width = 8, height = 20)

DotPlot.metabolism(countexp.Seurat, 
                   pathway = pathways, 
                   phenotype = "group",  #更改phenotype 参数
                   norm = "y")+
scale_x_discrete(limits = c("LowScore", "Medium", "HighScore")) 

dev.off()

#展示特定通路
input.pathway <- c(
  "Pyruvate metabolism",
  "Glycogen metabolism",
  "Glucose metabolism",
  "Arachidonic acid metabolism"
)

pdf('/home/luwei/UVM_output_new/REACTOME_scMetabolism.PDF', width = 4.2, height = 8)
DotPlot.metabolism(countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "group",  #更改phenotype 参数
                   norm = "y")+
scale_x_discrete(limits = c("LowScore", "Medium", "HighScore")) 
dev.off()

