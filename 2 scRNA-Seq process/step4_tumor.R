# 清空环境，加载依赖包
rm(list = ls())
options(stringsAsFactors = F)

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
library(SCP)
options(reticulate.conda_binary = "/home/zhangy/anaconda3/condabin/conda", SCP_env_name="SCP_env")
SCP::PrepareEnv()

# 设置地址
x = getwd()

setwd("/home/jiangxijie/UVM/step4/") 
set.seed(123)

# 加载上一步数据，
load("./data/step3_cell_annoation.Rdata")
lacGenes = read.gmt("./data/HP_INCREASED_SERUM_LACTATE.v2023.2.Hs.gmt")
sce.all = (scRNA)

table(sce.all$celltype)

# 挑选 tumor 细胞出来
table(sce.all$celltype)

sce1 = sce.all[, sce.all$celltype %in% c('Tumor')]
sce1 = JoinLayers(sce1)

as.data.frame(sce1@assays$RNA$counts[1:4,1:2])
head(sce1@meta.data,4)
table(sce1$orig.ident)

# 提取counts矩阵以及需要的metadata信息，重新走降维聚类分群流程
sce = sce1
sce = CreateSeuratObject(counts = sce@assays$RNA$counts,
                         meta.data = sce@meta.data[,1:4])
sce = NormalizeData(sce,
                    normalization.method = "LogNormalize",
                    scale.factor = 1e4)
sce = FindVariableFeatures(sce,selection.method = "vst",
                           nfeatures = 2000)
sce = ScaleData(sce)
sce = RunPCA(object = sce,pc.genes = VariableFeatures(sce))

DimHeatmap(sce,dims = 1:12,cells=100,balanced = TRUE)
ElbowPlot(sce,ndims = 50)

sce = FindNeighbors(sce,dims = 1:20)
for(res in c(0.01,0.05,0.1,0.2,0.3,0.5,0.8,1)){
    sce = FindClusters(sce,
                       resolution=res,
                       algorithm=1)
}
p2_tree = clustree(sce@meta.data,prefix="RNA_snn_res.")
p2_tree

ggsave(filename="./outcome/tree.pdf",plot = p2_tree)

table(sce@meta.data$RNA_snn_res.0.2)

sce = FindNeighbors(sce,dims = 1:20)
sce = FindClusters(sce,resolution = 0.2)

sce = RunTSNE(object = sce,dims = 1:20,do.fast=TRUE)
sce = RunUMAP(object = sce,dims = 1:20,do.fast=TRUE)

setwd("/home/jiangxijie/UVM/step4/") 
# 可视化
mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
p = DimPlot(sce, reduction = "tsne",cols = mycolors,pt.size = 0.8,
                  group.by = "RNA_snn_res.0.05",label = T,label.box = T)
p
ggsave(filename="/home/luwei/UVM_output_new/step4.tumor.tsne.pdf",plot = p)

FeaturePlot(sce,features ="EGF",reduction = "tsne",cols = c("lightgrey" ,"#DE1F1F"),raster=FALSE)

# 命名
celltype = data.frame(ClusterID = 0:10, celltype = 0:10)
celltype$celltype = paste("cluster",0:10)

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.05 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)

rownames(sce@assays$RNA$counts)

# 载入乳酸化模型数据
new_model <- data.frame(
  gene = c("TCEB1","PSMA7","PSMB8","PSMB9","PSMC5","PSMD11","PSME1",
           "PSME2","RBX1","GOT2","PYCRL","GNPDA1","HK1","HK3","NUP85",
           "PFKP","PGAM1","NME1","MAF1","DZIP1","ZFHX4","EBF3","ZNF831","RBPJ"),
  weight = c(0.0766,-0.0415,0.0007,0.0082,0.0251,0.0283,0.1058,
             0.0003,-0.0155,-0.1831,0.0464,0.259,0.0413,0.0341,0.275,
             0.0316,0.1368,0.1453,0.2088,0.0679,0.0462,0.1018,0.0142,0.2082)
)

# 找到基因在sce中的行号
set = match(new_model$gene, rownames(sce@assays$RNA$counts))
missing_genes <- new_model$gene[is.na(set)]
missing_genes

# 初始化score数据框
score = data.frame(cellname = colnames(sce@assays$RNA$counts),
                   cellscore = 1:ncol(sce@assays$RNA$counts))

valid_idx <- which(!is.na(set))
new_model <- new_model[valid_idx, ]
set <- set[valid_idx]

# 计算每个基因加权值
for (i in 1:nrow(new_model)) {
  x = set[i]                # 基因位置
  y = new_model$weight[i]   # 对应权重
  z = new_model$gene[i]     # 基因名
  score[,z] = sce@assays$RNA$counts[x,]*y
}

# 总分
score$cellscore = rowSums(score[,3:(2+nrow(new_model))])
head(score)

# 计算各个细胞群的乳酸化水平
# 将乳酸化分数写入 meta.data，并归一化
sce@meta.data$lacScore = (score$cellscore-min(score$cellscore))/(max(score$cellscore)-min(score$cellscore))

# 细胞分群可视化
# 提取tsne坐标
tsne_df <- as.data.frame(sce@reductions$tsne@cell.embeddings)
tsne_df$cluster <- as.factor(sce$celltype)

# 提取基因表达数据并与tsne坐标合并
gene_df <- as.data.frame(sce@meta.data$lacScore)
merged_df <- cbind(tsne_df, gene_df)
colnames(merged_df) = c('tSNE_1','tSNE_2','cluster','lacScore')
print(head(merged_df))

p = ggplot(merged_df, vars = c("tSNE_1", "tSNE_2", 'lacScore'), aes(x = tSNE_1, y = tSNE_2, colour = lacScore)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("white", "#930404"), limits = c(0, 1), oob = scales::squish) +
  new_scale_color()+
  theme_classic()
p
ggsave(filename="/home/luwei/UVM_output_new/step4.tumor.tsne.lacScore.pdf",plot = p)

saveRDS(sce,file = "/home/jiangxijie/UVM/step6/lacScore_new.rds")

print(head(merged_df))
sce@meta.data$group = ifelse(sce@meta.data$orig.ident %in% c("GSM4147091","GSM4147092","GSM4147100"),
                            "metastatic","primary")
print(head(sce@meta.data))


aggregate(lacScore ~ group, sce@meta.data, mean)
aggregate(lacScore ~ orig.ident, sce@meta.data, mean)
aggregate(lacScore ~ cluster, merged_df, mean)

x = aggregate(lacScore ~ cluster, merged_df, median)
x1 = arrange(x, lacScore)

# 对模型的乳酸化打分画直方图，判断大概位置和区分
p = gghistogram(merged_df, x = "lacScore",
            fill = "#3bc9db", # 设置填充色
            add = "mean", # 添加均值线
            rug = TRUE, # 添加轴须线
            binwidth = 0.02
)
p
ggsave(filename="/home/luwei/UVM_output_new/step4.histogram.pdf",plot = p,width = 7,height = 7)

floor(length(merged_df$lacScore)*0.2)
m = (merged_df$lacScore)
merged_df2 = merged_df
merged_df3 = merged_df


pdf("/home/luwei/UVM_output_new/step_arrangement.pdf")
for (i in seq(0.05, 0.5, by=0.05)){
m_top = sort(merged_df$lacScore,decreasing = T)[floor(length(merged_df$lacScore)*i)]
x = i*100
merged_df2$high_lac_cells = ifelse(merged_df$lacScore>m_top,paste("top",x,"%",sep = ""),"others")
merged_df3 = cbind(merged_df3,merged_df2$high_lac_cells)
colnames(merged_df3)[dim(merged_df3)[2]] = paste("top",x,"cells",sep = "")
    
p1 = ggplot(merged_df2, vars = c("tSNE_1", "tSNE_2", 'mean'), aes(x = tSNE_1, y = tSNE_2, colour = high_lac_cells)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_manual(values = c("lightgrey","#930404")) +
  new_scale_color()+
  theme_classic()+
  ggtitle(paste("lacScore divided by ",x,"%",sep = ""))
print(p1)
}

m_top = sort(merged_df$lacScore,decreasing = T)[floor(length(merged_df$lacScore)*0.25)]
x = m_top*100
merged_df2$high_lac_cells = ifelse(merged_df$lacScore>m_top,paste(">","mean",sep = ""),"others")
merged_df3 = cbind(merged_df3,merged_df2$high_lac_cells)
    
p1 = ggplot(merged_df2, vars = c("tSNE_1", "tSNE_2", 'mean'), aes(x = tSNE_1, y = tSNE_2, colour = high_lac_cells)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_manual(values = c("#930404","lightgrey")) +
  new_scale_color()+
  theme_classic()+
  ggtitle(paste("lacScore divided by ","mean",sep = ""))
print(p1)
dev.off()

p1

# 根据样本原本的分类，写入 class1 和 class2 的样本
sce@meta.data$class = ifelse(sce@meta.data$orig.ident %in% c("GSM4147091","GSM4147095","GSM4147098"),"class1","class2")
px1 = DimPlot(sce, reduction = "tsne",cols = mycolors,pt.size = 0.8,
                  group.by = "RNA_snn_res.0.05",label = T,label.box = T)

    

m_top = sort(merged_df$lacScore,decreasing = T)[floor(length(merged_df$lacScore)*0.25)]
x = m_top*100
merged_df2$high_lac_cells = ifelse(merged_df$lacScore>m_top,paste(">","25%",sep = ""),"others")
merged_df3 = cbind(merged_df3,merged_df2$high_lac_cells)
colnames(merged_df3)[dim(merged_df3)[2]] = "mean"
    
sce@meta.data$Top_LA_Cells = merged_df3$top25cells
px2 = DimPlot(sce, reduction = "tsne",cols = c("lightgrey","#930404"),pt.size = 0.8,
                  group.by = "Top_LA_Cells",label = F,label.box = F)+
   labs(title = "Top 25% LA Cells")



px3 = DimPlot(sce, reduction = "tsne",cols = mycolors,pt.size = 0.8,
                  group.by = "class",label = T,label.box = T)



print(px1)
print(px2)
print(px3)


p = px1+px2+px3
ggsave(filename="/home/luwei/UVM_output_new/step4.class.pdf",plot = p,width = 21,height = 7)
p = CellDimPlot(sce,
            group.by = "RNA_snn_res.0.05", 
            reduction = "tsne", 
            theme_use = "theme_blank",
            palcolor = mycolors,
            theme_args = list(base_size = 16))
p
ggsave(plot = p, filename="/home/luwei/UVM_output_new/1.pdf",width = 7,height = 7)

p = CellDimPlot(sce,
            group.by = "Top_LA_Cells", 
            reduction = "tsne", 
            theme_use = "theme_blank",
            palcolor = c("#930404","lightgrey"),
            theme_args = list(base_size = 16))
p
ggsave(plot = p, filename="/home/luwei/UVM_output_new/2.pdf",width = 7,height = 7)
p = CellDimPlot(sce,
            group.by = "class", 
            reduction = "tsne", 
            theme_use = "theme_blank",
            palcolor = c('#4DBBD4','#E64A35'),
            theme_args = list(base_size = 16))
p
ggsave(plot = p, filename="/home/luwei/UVM_output_new/3.pdf",width = 7,height = 7)

saveRDS(sce,file = "/home/luwei/UVM_output_new/step4.tumor.rds")

mycolors <-c('#4DBBD4','#6BD66B','#01A187','#3C5588' ,'#F29F80','#E64A35',
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD')

merged_df4 = merged_df3
colnames(merged_df4)
merged_df4 <- merged_df4[, !duplicated(colnames(merged_df4))]

q = quantile(merged_df4$lacScore, probs = rev(seq(0,1,0.25)))

merged_df4$cluster = factor(merged_df4$cluster,levels = rev(c("cluster 9","cluster 1","cluster 0","cluster 3","cluster 5","cluster 6",
"cluster 4","cluster 2","cluster 8","cluster 7")))

merged_df4$cells = ifelse(merged_df4$lacScore > q[[2]],"0~25%",
                          ifelse(merged_df4$lacScore < q[[4]],"76~100%","26~75%"))
merged_df4$cells = factor(merged_df4$cells,levels = c("76~100%","26~75%","0~25%"))






p = ggplot(data = merged_df4,aes(x = cluster,y = 1,fill = cells)) +
  geom_bar(stat = "identity",
           position = "fill") + 
#  coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  scale_fill_manual(values=c("#4DBBD4","#FEA82F",'#E64A35' ))+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1),  # 0到1之间每0.1一个刻度

                    expand = c(0, 0))+
  ylab("")+
  xlab("")+
  labs(fill = "Cell Rank")+
  geom_hline(aes(yintercept=0.6), colour="#2F4F4F", linetype="dashed")
p

ggsave(filename="/home/luwei/UVM_output_new/step4.barplot.pdf",plot = p,width = 5,height = 7)

start_time <- Sys.time() 
set.seed(123)
# 根据每一个 cluster随机抽样，每个 cluster 抽取 2/3 的细胞（大约只够 20000 个左右的细胞）
colnames(sce@meta.data)
Idents(sce) = "celltype"

exp_matrix = GetAssayData(sce, assay = "RNA", layer ="counts") #表达矩阵 sce@assays$RNA@layers$counts
metadata = sce@meta.data #metadata
row.names(sce@meta.data)[1:5]
dim(exp_matrix)
identical(row.names(sce@meta.data),colnames(exp_matrix))#TRUE

##开始抽样
cellnames = as.character(names(table(sce$celltype)))
sample_name = c()
#分层抽样
for (i in 1:length(cellnames)) {
    newmetadata = metadata[metadata$celltype == cellnames[i],]
    nsample = floor(table(sce$celltype)[i]/3*2)
    sample_name = c(sample_name,
                    row.names(newmetadata[sample(1:nrow(newmetadata),nsample,replace = F),]))
    sample_name = gsub("\\.[0-9]","",sample_name)
    rm(newmetadata)
}

# 选择抽出来的样本
df = as.data.frame(GetAssayData(sce, assay = "RNA", layer ="counts"))
df = df[,sample_name]
df[1:6,1:3]
dim(df)

idents = as.data.frame(t(as.data.frame(Idents(object = sce))))
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

## 输出single cell文件
#fwrite(final,file = "./outcome/step4.CibersortX.scRNAassay.txt",sep = "\t",row.names = F)
end_time <- Sys.time() 
end_time - start_time

# 读取 bulk 的数据，处理成 cibersortx 的输入样本
dat = read.table("./data/combat-dat.txt")
clin = read.table("./data/combat-clin.txt",head = T)
dat[1:6,1:6]
clin[1:6,1:4]

# 处理表型数据
## 给表型输入样本名称信息
clin$samples = colnames(dat)
identical(colnames(dat),clin$sample)

dat = rownames_to_column(dat,var = "gene")
head(dat)

## 去掉不必要的数据 
clin = clin[,-1]

head(clin)
head(dat)

# 输出 bulk 数据文件
write.table(dat,file = "./outcome/step4.CibersortX.dat.txt",sep = "\t",row.names = F)
write.table(clin,file = "./outcome/step4.CibersortX.clinical.txt",sep = "\t",row.names = F)

# 反卷积完成后读取数据（这一步设置了打分依据）
cell_prop <- fread("./data/CIBERSORTx_all_cluster_result.txt")#job15是反卷积到每个 cluster 上的结果

colnames(cell_prop)[1] = "samples"

cell_prop = merge(cell_prop,clin,by =  "samples",all = F)
cell_prop = as.data.frame(cell_prop)
rownames(cell_prop) = cell_prop$samples

cell_prop$top_clusters_proportion = cell_prop$`cluster 7`+cell_prop$`cluster 8`+cell_prop$`cluster 2`

dim(cell_prop)
cell_prop

# "E-MTAB-4097","GSE84976"按均值分组(OS)
cell_prop2 = subset(cell_prop,group %in% c("E-MTAB-4097","GSE84976"))

m = mean(cell_prop2$top_clusters_proportion)
cell_prop2$risk = ifelse(cell_prop2$top_clusters_proportion >= m,"HighRisk","LowRisk")
cell_prop2$risk = factor(cell_prop2$risk,levels = c("HighRisk","LowRisk"))

fit <- survfit(Surv(time,status) ~ risk, data = cell_prop2)
print(fit)

res.sum <- surv_summary(fit,data = cell_prop2)
head(res.sum[res.sum$risk=="HighRisk",])
head(res.sum[res.sum$risk=="LowRisk",])

surv_diff <- survdiff(Surv(time, status) ~ risk, data = cell_prop2)
p.val = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p.val

pdf("/home/luwei/UVM_output_new/step4.combat.OS.pdf",width = 5,height = 4)
p = ggsurvplot(fit,
               pval = TRUE, 
               conf.int = TRUE,
               surv.median.line = "hv",
               xlab = "Overall Survival (day)",
               legend.title = "",
               legend.labs = c("High Score", "Low Score"),
               palette = c('#E64A35', '#4DBBD4'),
               break.x.by = 1000,
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

# "GSE22138","GSE39717"按均值分组(MFS)
cell_prop2 = subset(cell_prop,group %in% c("GSE22138","GSE39717"))

m = mean(cell_prop2$top_clusters_proportion)
cell_prop2$risk = ifelse(cell_prop2$top_clusters_proportion >= m,"HighRisk","LowRisk")
cell_prop2$risk = factor(cell_prop2$risk,levels = c("HighRisk","LowRisk"))

fit <- survfit(Surv(time,status) ~ risk, data = cell_prop2)
print(fit)

res.sum <- surv_summary(fit,data = cell_prop2)
head(res.sum[res.sum$risk=="HighRisk",])
head(res.sum[res.sum$risk=="LowRisk",])

surv_diff <- survdiff(Surv(time, status) ~ risk, data = cell_prop2)
p.val = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p.val

pdf(paste0("/home/luwei/UVM_output_new/step4.combat.MFS.pdf"),width = 5,height = 4)
p = ggsurvplot(fit,
               pval = TRUE, 
               conf.int = TRUE,
               surv.median.line = "hv",
               xlab = "Metastasis Free Survival (day)",
               legend.title = "",
               legend.labs = c("High Score", "Low Score"),
               palette = c('#E64A35', '#4DBBD4'),
               break.x.by = 1000,
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

# 读取 bulk 的数据，处理成 cibersortx 的输入样本
dat = read.table("./data/TCGA-UVM-dat.txt")

dat[1:6,1:6]

dat = rownames_to_column(dat,var = "gene")
head(dat)

# 输出 bulk 数据文件
write.table(dat,file = "/home/luwei/UVM_output_new/step4.CibersortX.TCGA.dat.txt",sep = "\t",row.names = F)

# 反卷积完成后读取数据
cell_prop <- fread("./data/CIBERSORTx_all_cluster_TCGA_result.txt")#job15是反卷积到每个 cluster 上的结果
clin = fread("./data/TCGA_clin.txt")
cell_prop
clin

colnames(cell_prop)[1] = "samples"
colnames(clin)[1] = "samples"


cell_prop = merge(cell_prop,clin,by =  "samples",all = F)
cell_prop = as.data.frame(cell_prop)
rownames(cell_prop) = cell_prop$samples

cell_prop$top_clusters_proportion = cell_prop$`cluster 7`+cell_prop$`cluster 8`+cell_prop$`cluster 2`

dim(cell_prop)
cell_prop

# TCGA按均值分组(OS)
cell_prop2 = cell_prop

m = mean(cell_prop2$top_clusters_proportion)
cell_prop2$risk = ifelse(cell_prop2$top_clusters_proportion >= m,"HighRisk","LowRisk")
cell_prop2$risk = factor(cell_prop2$risk,levels = c("HighRisk","LowRisk"))

fit <- survfit(Surv(OS.time,OS) ~ risk, data = cell_prop2)
print(fit)

res.sum <- surv_summary(fit,data = cell_prop2)
head(res.sum[res.sum$risk=="HighRisk",])
head(res.sum[res.sum$risk=="LowRisk",])

surv_diff <- survdiff(Surv(OS.time, OS) ~ risk, data = cell_prop2)
p.val = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p.val

pdf(paste0("/home/luwei/UVM_output_new/step4.TCGA.OS.pdf"),width = 5,height = 4)
p = ggsurvplot(fit,
               pval = TRUE, 
               conf.int = TRUE,
               surv.median.line = "hv",
               xlab = "Overall Survival (day)",
               legend.title = "",
               legend.labs = c("High Score", "Low Score"),
               palette = c('#E64A35', '#4DBBD4'),
               break.x.by = 500,
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

# TCGA按均值分组(MFS)
cell_prop2 = cell_prop

m = mean(cell_prop2$top_clusters_proportion)
cell_prop2$risk = ifelse(cell_prop2$top_clusters_proportion >= m,"HighRisk","LowRisk")
cell_prop2$risk = factor(cell_prop2$risk,levels = c("HighRisk","LowRisk"))

fit <- survfit(Surv(PFI.time,PFI) ~ risk, data = cell_prop2)
print(fit)

res.sum <- surv_summary(fit,data = cell_prop2)
head(res.sum[res.sum$risk=="HighRisk",])
head(res.sum[res.sum$risk=="LowRisk",])

surv_diff <- survdiff(Surv(PFI.time, PFI) ~ risk, data = cell_prop2)
p.val = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p.val

pdf(paste0("/home/luwei/UVM_output_new/step4.TCGA.PFI.pdf"),width = 5,height = 4)
p = ggsurvplot(fit,
               pval = TRUE, 
               conf.int = TRUE,
               surv.median.line = "hv",
               xlab = "Progression Free Interval (day)",
               legend.title = "",
               legend.labs = c("High Score", "Low Score"),
               palette = c('#E64A35', '#4DBBD4'),
               break.x.by = 500,
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

# TCGA按均值分组(DSS)
cell_prop2 = cell_prop

m = median(cell_prop2$top_clusters_proportion)
cell_prop2$risk = ifelse(cell_prop2$top_clusters_proportion >= m,"HighRisk","LowRisk")
cell_prop2$risk = factor(cell_prop2$risk,levels = c("HighRisk","LowRisk"))

fit <- survfit(Surv(DSS.time,PFI) ~ risk, data = cell_prop2)
print(fit)

res.sum <- surv_summary(fit,data = cell_prop2)
head(res.sum[res.sum$risk=="HighRisk",])
head(res.sum[res.sum$risk=="LowRisk",])

surv_diff <- survdiff(Surv(DSS.time, DSS) ~ risk, data = cell_prop2)
p.val = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p.val

pdf(paste0("/home/luwei/UVM_output_new/step4.TCGA.DSS.pdf"),width = 5,height = 4)
p = ggsurvplot(fit,
               pval = TRUE, 
               conf.int = TRUE,
               surv.median.line = "hv",
               xlab = "Disease Specific Survival (day)",
               legend.title = "",
               legend.labs = c("High Score", "Low Score"),
               palette = c('#E64A35', '#4DBBD4'),
               break.x.by = 500,
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

saveRDS(sce, "/home/jiangxijie/UVM/step4/step4.single.cell.tumor.new.rds")