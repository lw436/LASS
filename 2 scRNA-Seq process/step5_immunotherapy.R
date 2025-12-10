# 清空环境，加载依赖包
rm(list = ls())
options(stringsAsFactors = F,future.globals.maxSize = 172 * 1024^3)

.libPaths(c(
  "/home/zhangy/anaconda3/envs/single/lib/R/library",
  "/home/zhangy/anaconda3/envs/hdWGCNA/lib/R/library",
  "/home/luwei/R/x86_64-pc-linux-gnu-library/4.4",
  "/usr/local/lib/R/site-library",
  "/usr/lib/R/site-library",
  "/usr/lib/R/library",
  "/home/zhangy/R/x86_64-conda-linux-gnu-library/4.4"
))

library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(ggnewscale)
library(UCell)
library(reshape2)
library(stringr)
library(clusterProfiler)
library(ggpubr)
library(ggthemes)
library(ggnewscale)
library(readr)
library(data.table)
library(stringr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(survival)
library(survminer)
library(cowplot)
library(patchwork)
library(future)
library(plyr)


# 设置地址
x = getwd()
setwd(x)
set.seed(123)


mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD')

# 加载上一步数据
load("./data/step3_cell_annoation.Rdata")
sce.all = (scRNA)
table(sce.all$celltype)

sce1 = sce.all[, sce.all$celltype %in% c('Tumor')]
sce1 = JoinLayers(sce1)

as.data.frame(sce1@assays$RNA$counts[1:4,1:2])
head(sce1@meta.data,4)
table(sce1$orig.ident)

sce.tumor = sce1

# 提取counts矩阵以及需要的metadata信息，重新走降维聚类分群流程
sce.tumor = CreateSeuratObject(counts = sce.tumor@assays$RNA$counts,
                         meta.data = sce.tumor@meta.data[,1:4])
sce.tumor = NormalizeData(sce.tumor,
                    normalization.method = "LogNormalize",
                    scale.factor = 1e4)
sce.tumor = FindVariableFeatures(sce.tumor,selection.method = "vst",
                           nfeatures = 2000)
sce.tumor = ScaleData(sce.tumor)
sce.tumor = RunPCA(object = sce.tumor,pc.genes = VariableFeatures(sce.tumor))

DimHeatmap(sce.tumor,dims = 1:12,cells=100,balanced = TRUE)
ElbowPlot(sce.tumor,ndims = 50)

sce.tumor = FindNeighbors(sce.tumor,dims = 1:20)
for(res in c(0.01,0.05,0.1,0.2,0.3,0.5,0.8,1)){
    sce.tumor = FindClusters(sce.tumor,
                       resolution=res,
                       algorithm=1)
}
p2_tree = clustree(sce.tumor@meta.data,prefix="RNA_snn_res.")
p2_tree
table(sce.tumor@meta.data$RNA_snn_res.0.2)

sce.tumor = FindNeighbors(sce.tumor,dims = 1:20)
sce.tumor = FindClusters(sce.tumor,resolution = 0.2)

sce.tumor = RunTSNE(object = sce.tumor,dims = 1:20,do.fast=TRUE)
sce.tumor = RunUMAP(object = sce.tumor,dims = 1:20,do.fast=TRUE)

# 打分
# 命名
celltype = data.frame(ClusterID = 0:10, celltype = 0:10)
celltype$celltype = paste("cluster",0:10)

sce.tumor@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.tumor@meta.data[which(sce.tumor@meta.data$RNA_snn_res.0.05 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.tumor@meta.data$celltype)

# 载入乳酸化模型数据
new_model <- data.frame(
  gene = c("TCEB1","PSMA7","PSMB8","PSMB9","PSMC5","PSMD11","PSME1",
           "PSME2","RBX1","GOT2","PYCRL","GNPDA1","HK1","HK3","NUP85",
           "PFKP","PGAM1","NME1","MAF1","DZIP1","ZFHX4","EBF3","ZNF831","RBPJ"),
  weight = c(0.0766,-0.0415,0.0007,0.0082,0.0251,0.0283,0.1058,
             0.0003,-0.0155,-0.1831,0.0464,0.259,0.0413,0.0341,0.275,
             0.0316,0.1368,0.1453,0.2088,0.0679,0.0462,0.1018,0.0142,0.2082)
)

# 找到基因在sce.tumor中的行号
set = match(new_model$gene, rownames(sce.tumor@assays$RNA$counts))
missing_genes <- new_model$gene[is.na(set)]
missing_genes

# 初始化score数据框
score = data.frame(cellname = colnames(sce.tumor@assays$RNA$counts),
                   cellscore = 1:ncol(sce.tumor@assays$RNA$counts))

valid_idx <- which(!is.na(set))
new_model <- new_model[valid_idx, ]
set <- set[valid_idx]

# 计算每个基因加权值
for (i in 1:nrow(new_model)) {
  x = set[i]                # 基因位置
  y = new_model$weight[i]   # 对应权重
  z = new_model$gene[i]     # 基因名
  score[,z] = sce.tumor@assays$RNA$counts[x,]*y
}

# 总分
score$cellscore = rowSums(score[,3:(2+nrow(new_model))])
head(score)

# 计算各个细胞群的乳酸化水平
# 将乳酸化分数写入 meta.data，并归一化
sce.tumor@meta.data$lacScore = (score$cellscore-min(score$cellscore))/(max(score$cellscore)-min(score$cellscore))

# 细胞分群可视化
# 提取tsne坐标
tsne_df <- as.data.frame(sce.tumor@reductions$tsne@cell.embeddings)
tsne_df$cluster <- as.factor(sce.tumor$celltype)

# 提取基因表达数据并与tsne坐标合并
gene_df <- as.data.frame(sce.tumor@meta.data$lacScore)
merged_df <- cbind(tsne_df, gene_df)
colnames(merged_df) = c('tSNE_1','tSNE_2','cluster','lacScore')
head(merged_df)

p = ggplot(merged_df, vars = c("tSNE_1", "tSNE_2", 'lacScore'), aes(x = tSNE_1, y = tSNE_2, colour = lacScore)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("white", "#930404"), limits = c(0, 1), oob = scales::squish) +
  new_scale_color()+
  theme_classic()
p

x = aggregate(lacScore ~ cluster, merged_df, median)
x1 = arrange(x, lacScore)

x1

Tumor_names_all = rownames(subset(sce.all@meta.data,celltype %in% "Tumor"))
Tumor_names_tumor = rownames(sce.tumor@meta.data)
identical(Tumor_names_all,Tumor_names_tumor)

start_time <- Sys.time() 
sce.all.filter = sce.all
for (i in 1:ncol(sce.all.filter)){
    name = rownames(sce.all.filter@meta.data[i,])
    if(sce.all.filter@meta.data[i,]$celltype %in% "Tumor"){
        sce.all.filter@meta.data[i,]$celltype = sce.tumor@meta.data[name,]$celltype        
    }  
}
end_time <- Sys.time() 
end_time - start_time

saveRDS(object = sce.all.filter, file = "/home/jiangxijie/UVM/step5/step5.single.cell.tumor.new.rds")
start_time <- Sys.time() 
sce.all.filter = sce.all
for (i in 1:ncol(sce.all.filter)){
    name = rownames(sce.all.filter@meta.data[i,])
    if(sce.all.filter@meta.data[i,]$celltype %in% "Tumor"){
        sce.all.filter@meta.data[i,]$celltype = sce.tumor@meta.data[name,]$celltype        
    }  
}
end_time <- Sys.time() 
end_time - start_time

saveRDS(object = sce.all.filter, file = "/home/jiangxijie/UVM/step5/step5.single.cell.tumor.new.rds")

set.seed(123)
# 根据每一个 cluster随机抽样，每个 cluster 抽取 2/3 的细胞（大约只够 20000 个左右的细胞）
colnames(sce.all.filter@meta.data)
Idents(sce.all.filter) = "celltype"

exp_matrix = GetAssayData(sce.all.filter, assay = "RNA", layer ="counts") #表达矩阵 sce.all.filter@assays$RNA@layers$counts
metadata = sce.all.filter@meta.data #metadata
row.names(sce.all.filter@meta.data)[1:5]
dim(exp_matrix)
identical(row.names(sce.all.filter@meta.data),colnames(exp_matrix))#TRUE

##开始抽样
cellnames = as.character(names(table(sce.all.filter$celltype)))
sample_name = c()
#分层抽样
for (i in 1:length(cellnames)) {
    newmetadata = metadata[metadata$celltype == cellnames[i],]
    nsample = floor(table(sce.all.filter$celltype)[i]/2)
    sample_name = c(sample_name,
                    row.names(newmetadata[sample(1:nrow(newmetadata),nsample,replace = F),]))
    sample_name = gsub("\\.[0-9]","",sample_name)
    rm(newmetadata)
}

# 选择抽出来的样本
df = as.data.frame(GetAssayData(sce.all.filter, assay = "RNA", layer ="counts"))
df = df[,sample_name]
df[1:6,1:3]
dim(df)

idents = as.data.frame(t(as.data.frame(Idents(object = sce.all.filter))))
idents = idents[,sample_name]
identical(colnames(df), colnames(idents))
colnames(df) = idents[1, ]
x = data.frame(GeneSymbol = rownames(df))
rownames (df) = NULL 
final = cbind(x, df)
final[1:6,1:6]

keep <- rowSums(final>0) >= floor(0.5*ncol(final))
table(keep)
filter_count <- final[keep,] 
filter_count[1:6,1:6]
dim(filter_count)

load("./data/step5.SKCM.clin.Rdata")

# 反卷积完成后读取数据
cell_prop <- fread("./data/step5_CIBERSORTx_skcm_result.txt")#job15是反卷积到每个 cluster 上的结果
colnames(cell_prop)[1] = "sample"

cell_prop = merge(cell_prop,clin,by =  "sample",all = F)

cell_prop = as.data.frame(cell_prop)
rownames(cell_prop) = cell_prop$sample

cell_prop$top_clusters_proportion = cell_prop$`cluster 7`+cell_prop$`cluster 8`+cell_prop$`cluster 2`

dim(cell_prop)
cell_prop

cell_prop$dataset = factor(cell_prop$dataset,levels = c("GSE158403_1","GSE115821_2","GSE91061_3","GSE78220_4","GSE168204_5"))

table(cell_prop$dataset)

x = cell_prop[cell_prop$dataset %in% "GSE91061_3",]
write.csv(x,file = "./outcome/GSE91061_3_clin.csv")

names(head(cell_prop))

cell_prop2 = cell_prop[,c("dataset","timepoint","resp","top_clusters_proportion")] 
cell_prop2 = cell_prop2[!is.na(cell_prop2$resp),]
cell_prop2 = cell_prop2[cell_prop2$dataset %in% "GSE91061_3",]
cell_prop2 = cell_prop2[cell_prop2$timepoint %in% "on",]
head(cell_prop2)
dim(cell_prop2)

cell_prop3 = cell_prop2

m = mean(cell_prop2$top_clusters_proportion)
cell_prop2$risk = ifelse(cell_prop2$top_clusters_proportion >= m,"HighRisk","LowRisk")


p = ggplot(data = cell_prop2,aes(x = risk,y = 1,fill = resp)) +
  geom_bar(stat = "identity",
           position = "fill") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  scale_fill_manual(values=c('#4DBBD4','#E64A35'))+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("")+
  xlab("")
p

contingency_table <- table(cell_prop2$risk, cell_prop2$resp)

# 执行卡方检验
chi_result <- chisq.test(contingency_table)
print(chi_result)

 
#R 函数 fisher.test() 用于执行 Fisher's exact test，详情 ?fisher.test
contingency_table
fisher.test(contingency_table, alternative = 'two.sided', conf.level = 0.95)
stack_data = data.frame(cell_prop2[,c("risk","resp")])
a <- data.frame(table(stack_data$risk,stack_data$resp))
a <- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
a

# 绘制两组免疫治疗响应比例
a$Var1 = factor(a$Var1,levels = c("HighRisk","LowRisk"))

pvalue = fisher.test(contingency_table, alternative = 'two.sided', conf.level = 0.95)$p.value
squared = fisher.test(contingency_table, alternative = 'two.sided', conf.level = 0.95)$estimate

p = ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack())+
  scale_fill_manual(values = c("#DB423E","#008ECA"),label=c("SD/PD","CR/PR"))+
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = scales::percent_format(scale = 1))+ #百分比y轴
  labs(x="Group",y="Percent Weidght",
       fill="")+
  annotate(geom = "text",
           cex=4,
           x=1.5, y=105, # 根据自己的数据调节p value的位置
           label=paste0("P ", ifelse(pvalue<0.001, "< 0.0001", paste0("= ",round(pvalue,3)))), # 添加P值
           color="black")+

  annotate(geom = "text",
           cex=4,
           x=1.5, y=110, # 根据自己的数据调节p value的位置
           label=paste0("Odds Ratio = ",round(squared,3),sep = ""), # 添加P值
           color="black")+

  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10))+
  annotate("text", x = 2 , y = 19,label = a$label[4],colour="black")+
  annotate("text", x = 2 , y = 69,label = a$label[3],colour="black")+
  annotate("text", x = 1 , y = 52,label = a$label[1],colour="black")+
  annotate("text", x = 1 , y = 2,label = a$label[2],colour="black")+
  scale_x_discrete("", labels = c("N" = "Nomal","H" = "Hypertension"))

p

ggsave(plot = p, filename="/home/luwei/UVM_output_new/step5.barplot.percentage.pdf",width = 5,height = 7)

# 免疫治疗数据按均值分组
cell_prop2 = cell_prop[!is.na(cell_prop$OS.time),]
cell_prop2 = cell_prop2[cell_prop2$dataset %in% "GSE91061_3",]
cell_prop2 = cell_prop2[cell_prop2$timepoint %in% "on",]
head(cell_prop2)
dim(cell_prop2)

m = mean(cell_prop2$top_clusters_proportion)
cell_prop2$risk = ifelse(cell_prop2$top_clusters_proportion > m,"HighRisk","LowRisk")
cell_prop2$risk = factor(cell_prop2$risk,levels = c("HighRisk","LowRisk"))

fit <- survfit(Surv(OS.time,OS) ~ risk, data = cell_prop2)
print(fit)


pdf("/home/luwei/UVM_output_new/step5.melanoma.immu.threapy.OS(on).pdf",width = 5,height = 4)
p = ggsurvplot(fit,
               pval = TRUE, 
               conf.int = TRUE,
               surv.median.line = "hv",
               xlab = "Overall Survival (day)",
               legend.title = "",
               legend.labs = c("High Score", "Low Score"),
               palette = c('#E64A35', '#4DBBD4'),
               break.x.by = 300,
               font.x = c(14), # 设置x轴字体大小、格式和颜色
               font.y = c(14),# 设置y轴字体大小、格式和颜色
               risk.table = TRUE, # 设置添加风险因子表
               tables.height = 0.2, # 设置风险表的高度
               tables.theme = theme_cleantable(), # 设置风险表的主题
               ggtheme = theme_bw())
p1 = p$plot+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(axis.text.y=element_text(size=12,vjust=.5),
        axis.text.x=element_text(size=12,vjust=.5))


p1

dev.off()
p1