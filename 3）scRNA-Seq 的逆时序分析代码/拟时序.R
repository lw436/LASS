# 准备步骤
# 清空环境，加载依赖包
rm(list = ls())
options(stringsAsFactors = F,future.globals.maxSize = 128000 * 1024^3)
.libPaths(c(
  "/home/zhangy/anaconda3/envs/single/lib/R/library",
  "/home/zhangy/anaconda3/envs/hdWGCNA/lib/R/library",
  "/home/luwei/R/x86_64-pc-linux-gnu-library/4.4",
  "/usr/local/lib/R/site-library",
  "/usr/lib/R/site-library",
  "/usr/lib/R/library",
  "/home/zhangy/R/x86_64-conda-linux-gnu-library/4.4"
))

library(future)
library(tidyverse)
library(tinyarray)
library(data.table) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(monocle)
library(ggpubr)


# 设置地址
setwd("/home/luwei/UVM_output_new")
set.seed(123)

availableCores() #12 #查看几个核可用
plan(multisession, workers=36) #观察htop的CPU，发现中后期确实使用了多核
nbrOfWorkers() #4 当前可用的核有多少个

mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD')

setwd("/home/jiangxijie/UVM/step6")
Mono.cds = readRDS("./process_documentation/Mono.cds.DDRtree.rds")
lacScore = readRDS("/home/jiangxijie/UVM/step6/lacScore_new.rds")
pData(Mono.cds) = pData(Mono.cds)[,-18]
pd <- pData(Mono.cds)
pd$lacScore <- lacScore@meta.data[rownames(pd), "lacScore"]  # 按行名对齐
pData(Mono.cds) <- pd

# 展示Cluster/Pseudotime轨迹分布图：
p1 = plot_cell_trajectory(Mono.cds,color_by="celltype", size=1,show_backbone=TRUE)
p2 = plot_cell_trajectory(Mono.cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 

p2

ggsave(plot = p2, filename="/home/luwei/UVM_output_new/step6.Pseudotime.pdf",width = 7,height = 7)

# 风险三分群图
pData(Mono.cds)$celltype <- factor(
  pData(Mono.cds)$celltype,
  levels = c("cluster 7","cluster 8","cluster 2","cluster 4","cluster 6",
             "cluster 5","cluster 3","cluster 0","cluster 1","cluster 9")
)

# 设定 risk
pData(Mono.cds)$risk <- ifelse(
  pData(Mono.cds)$celltype %in% c("cluster 7","cluster 8","cluster 2"), 
  "High Lactylation Cluster",
  ifelse(
    pData(Mono.cds)$celltype %in% c("cluster 4","cluster 6","cluster 5"),
    "Medium Lactylation Cluster",
    "Low Lactylation Cluster"
  )
)

p <- plot_cell_trajectory(Mono.cds, color_by = "celltype") +
     facet_wrap("~risk", nrow = 1)

ggsave("/home/luwei/UVM_output_new/step6.lac.cluster.pdf",
       plot = p, width = 21, height = 7)

p

# 转移/未转移
pData(Mono.cds)$orig.ident <- factor(
  pData(Mono.cds)$orig.ident,
  levels = c("GSM4147093","GSM4147094","GSM4147095","GSM4147096",
             "GSM4147097","GSM4147098","GSM4147099","GSM4147101"))

# 设定 risk
pData(Mono.cds)$class <- ifelse(
  pData(Mono.cds)$orig.ident %in% c("GSM4147095","GSM4147098"), 
  "class1","class2")

p <- plot_cell_trajectory(Mono.cds, color_by = "class")

ggsave("/home/luwei/UVM_output_new/step6.lac.class.pdf",
       plot = p, width = 7, height = 7)

p

# 筛选基因,这里可以根据自己的需要筛选特定的基因：
disp_table <- dispersionTable(Mono.cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
Mono.cds <- setOrderingFilter(Mono.cds, unsup_clustering_genes$gene_id)

# 拟时相关基因聚类热图：
#高变基因
start_time <- Sys.time() 

if(F){
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)

diff_test <- differentialGeneTest(Mono.cds[disp.genes,], cores = 36, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
    
saveRDS(diff_test, "./process_documentation/diff_test.rds") 
}

end_time <- Sys.time() 
end_time - start_time

diff_test = readRDS("./process_documentation/diff_test.rds")

gene_to_cluster = diff_test %>% 
  arrange(qval) %>% 
  head(200) %>% 
  pull(gene_short_name)
head(gene_to_cluster)

p = plot_pseudotime_heatmap(Mono.cds[gene_to_cluster,], 
                            num_clusters=3,
                            show_rownames=T, 
                            return_heatmap=T,
                            trend_formula = "~sm.ns(Pseudotime, df=3)")
pdf("/home/luwei/UVM_output_new/pseudotime_heatmap_200.pdf", width = 10, height = 12)
p
dev.off()