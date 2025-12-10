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
library(data.table)
library(dplyr)
library(harmony)
library(patchwork)
library(plotrix)
library(ggsci)
library(celldex)
library(singleseqgset)
library(devtools)
library(tidyr)
library(reshape2)
library(gplots)
library(ggthemes)
library(ggnewscale)
library(ggpubr)



# 设置地址
setwd("/home/jiangxijie/UVM/step3/")


# 加载上一步数据，
load("./data/step2_data_filtering.Rdata")
sce.all.int
head(phe)

# 选择 0.2
sel.clust = "RNA_snn_res.0.2"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)

table(sce.all.int@active.ident) #看看每一分群的细胞数量
colnames(sce.all.int@meta.data) 
scRNA = sce.all.int

# 添加注释信息
Bcells <- c("PTPRC", "MS4A1", "CD79A", "CD19" )
Dendritic <- c("PTPRC", "IL3RA", "IRF7", "IRF8", "GZMB", "CD4", "CLEC4C", "PTGDS", "JCHAIN", "PLAC8", "PLD4", "TCF4", "BCL11A", "GPR183", "CCDC50", "LILRA4", "TSPAN13", "CLIC3", "MPEG1")
Endothelial <- c("CLEC14A", "PECAM1", "VWF", "CAV1", "EMCN", "CDH5", "MCAM", "IL3RA", "IGFBP7", "COL4A1", "COL4A2", "COL15A1", "SPARCL1", "A2M", "HSPG2", "PLVAP", "AQP1", "ENG", "RAMP2", "GNG11", "EGFL7", "CLDN5", "INSR")
Fibroblast <- c("COL1A1", "COL3A1", "WT1", "ACTA2", "CAV1", "COL1A2", "DCN", "SPARC", "COL6A1", "CCDC80", "LUM", "COL6A2", "COL6A3", "CALD1", "RARRES2", "MGP", "CTHRC1", "AEBP1", "POSTN", "COL5A2", "FBLN1", "TAGLN", "C1S", "C1R", "NNMT", "MMP2", "IGFBP5", "TIMP1", "FN1", "IGFBP7", "C3", "COL5A1", "LGALS1")
Myeloid <- c("PTPRC", "CD14", "FCER1G", "FCGR3A", "LYZ", "CTSS", "CD33", "CD68", "CD163", "ITGAX", "ITGAM", "CD4", "MRC1", "VSIG4", "SPP1", "APOE", "C1QA", "C1QB", "C1QC", "APOC1", "FTL", "S100A9", "TYROBP", "AIF1", "CD74", "PSAP", "CTSB")
Epi <- c("WFDC2", "CD24", "CLDN3", "KRT7", "KRT8", "KRT17", "KRT18", "KRT19", "EPCAM", "WT1", "CLDN4", "MSLN", "FOLR1", "MUC1")
Plasma <- c("PTPRC", "IGKC", "IGHG1", "CD79A", "IGHG2", "IGLC2", "IGLC3", "IGHG3", "IGHG4", "JCHAIN", "MZB1", "XBP1")
Tcells <- c("PTPRC", "CD2", "CD3D", "TRAC", "GZMA", "NKG7", "CD3E", "CD3G", "CD4", "TCF7", "CD8A", "PRF1", "GZMB", "CCL5", "CCL4", "IL32", "CD52")
Mast <- c("PTPRC", "KIT", "CPA3", "CTSG", "MS4A2", "TPSAB1", "TPSB2", "HPGD", "HPGDS", "GATA2")
SMC=c('NOTCH3','RGS5','NDUFA4L2','MYH11','COX4I2','PLN')
cycle=c('RRM2','MKI67','BIRC5','UBE2C','TOP2A','AURKB')

cells = list(Bcells,Dendritic,Endothelial,Fibroblast,Myeloid,Epi,Plasma,Tcells,Mast,SMC,cycle)
names(cells) = c("Bcells","Dendritic","Endothelial","Fibroblast","Myeloid","Epi","Plasma","Tcells","Mast","SMC","cycle")



pdf("./outcome/step3.cellMarkers.pdf")
for (i in 1:length(cells)){
    p = DotPlot(scRNA, features = unique(cells[[i]]),
                assay='RNA'  )  + coord_flip() + ggtitle(names(cells)[i])
    print(p)

}
dev.off()




genes_to_check = c("MLANA", "MITF","DCT","PRAME","GEP",#Tumor
                   "CD3D", "CD3E", "CD8A",#T
                   "CD19", "MS4A1","CD79A",#B
                   "IGHG1", "MZB1", "SDC1",#plasma
                   "CD68", "CD163", "CD14",#mono M
                   "FGFBP2", "FCGR3A", "CX3CR1",#NK
                   "RPE65",#Retinal pigment epithelium
                   "RCVRN",#Photoreceptor cells
                   "FGF7",#Fibroblasts
                   "PECAM1", "VWF"#Endothelial cells
                  )  
p1 = DotPlot(scRNA, features = unique(genes_to_check),
                assay='RNA'  )  + coord_flip()
p1
ggsave(plot = p1, filename="./outcome/step3.dotplot.pdf",width = 7,height = 7)

# 2 类肿瘤和 1 类的区别
genes_to_check = c("CDH1","ECM1","HTR2B","RAB31",#up
                   "EIF1B","FXR1","ID2","LMCD1","LTA4H","MTUS1","ROBO1","SATB1",#down
                   "FZD6","PHLDA1","ENPP2","BAP1"
                  )  
p = DotPlot(scRNA, features = unique(genes_to_check),
                assay='RNA'  )  + coord_flip()
p

# 重新聚类
mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
tsne =DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 0.8,
                  group.by = "RNA_snn_res.0.2",label = T,label.box = T)
tsne

ggsave(plot = tsne, filename="./outcome/step3.tsne.0.2.pdf",width = 7,height = 7)

# 重新分组
num = length(table(scRNA@meta.data$RNA_snn_res.0.2))-1
celltype=data.frame(ClusterID=0:num, celltype= 0:num) 

celltype[celltype$ClusterID %in% c(0,1,2,4,6,8),2]='Tumor'
celltype[celltype$ClusterID %in% c(3,9),2]='T cell'
celltype[celltype$ClusterID %in% c(5),2]='Macrophage'
celltype[celltype$ClusterID %in% c(7),2]='RPE cell'
celltype[celltype$ClusterID %in% c(10),2]='Endothelial'

# 将细胞类型信息添加到meta.data
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.2 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype)

scRNA@meta.data$sample = ifelse(scRNA@meta.data$orig.ident %in% c("GSM4147091","GSM4147092","GSM4147100"),"metastatic","primary")

# 可视化
celltype_tsne = DimPlot(scRNA, 
                       reduction = "tsne",
                       cols = mycolors,
                       pt.size = 1,
                       group.by = "celltype",
                       label = T)

patient_tsne = DimPlot(scRNA, reduction = "tsne",
                      cols = mycolors,
                      pt.size = 0.2,
                      group.by = "orig.ident") 

sample_tsne =DimPlot(scRNA, reduction = "tsne",
                     cols = mycolors,
                     pt.size = 0.2,
                     group.by = "RNA_snn_res.0.2") 
patient_tsne + sample_tsne + celltype_tsne

# 保存图片
p2 = patient_tsne + sample_tsne + celltype_tsne
p2
ggsave(plot = p2, filename="./outcome/step3.tsne.group.pdf",width = 21,height = 7)

library(SCP)
options(reticulate.conda_binary = "/home/zhangy/anaconda3/condabin/conda", SCP_env_name="SCP_env")
SCP::PrepareEnv()

p = CellDimPlot(scRNA,
            group.by = "orig.ident", 
            reduction = "tsne", 
            theme_use = "theme_blank",
            palcolor = mycolors,
            theme_args = list(base_size = 16))
ggsave(plot = p, filename="./outcome/1.pdf",width = 7,height = 7)

p = CellDimPlot(scRNA,
            group.by = "RNA_snn_res.0.2", 
            reduction = "tsne", 
            theme_use = "theme_blank",
            palcolor = mycolors,
            theme_args = list(base_size = 16))
ggsave(plot = p, filename="./outcome/2.pdf",width = 7,height = 7)

p = CellDimPlot(scRNA,
            group.by = "celltype", 
            reduction = "tsne", 
            theme_use = "theme_blank",
            palcolor = c('#6BD66B','#3C5588','#4DBBD4' ,'#E64A35' ,'#01A187'),
            theme_args = list(base_size = 16))
ggsave(plot = p, filename="./outcome/3.pdf",width = 7,height = 7)

head(scRNA@meta.data$RNA_snn_res.0.2)

# 看看感兴趣的基因在群中的分布
FeaturePlot(scRNA,
            features = c("FZD6","PHLDA1","ENPP2","BAP1"),
            reduction = "tsne",
            cols = c("lightgrey" ,"#DE1F1F"),
            ncol=2,
            raster=FALSE)

tb = table(scRNA$sample, scRNA$celltype)
head(tb)

pdf("./outcome/step3.cell.number.pdf",width = 7,height = 7)
balloonplot(tb)
dev.off()

bar_data <- as.data.frame(tb)
bar_per <- bar_data %>% 
  group_by(Var1) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
write.csv(bar_per,file = "./outcome/celltype_by_group_percent.csv")
col =c("#3176B7","#F78000","#3FA116","#CE2820","#9265C1",
       "#885649","#DD76C5","#BBBE00","#41BED1")
colnames(bar_per)

p4 <- ggplot(bar_per, aes(y = percent, x = Var1)) +
  geom_bar(aes(fill = Var2) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Relative proportion(%)')+
  scale_fill_manual(values=col)+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4
ggsave(plot = p4, filename="./outcome/step3.tsne.cell.precentage.pdf",width = 7,height = 7)

new_model <- data.frame(
  gene = c("TCEB1","PSMA7","PSMB8","PSMB9","PSMC5","PSMD11","PSME1",
           "PSME2","RBX1","GOT2","PYCRL","GNPDA1","HK1","HK3","NUP85",
           "PFKP","PGAM1","NME1","MAF1","DZIP1","ZFHX4","EBF3","ZNF831","RBPJ"),
  weight = c(0.0766,-0.0415,0.0007,0.0082,0.0251,0.0283,0.1058,
             0.0003,-0.0155,-0.1831,0.0464,0.259,0.0413,0.0341,0.275,
             0.0316,0.1368,0.1453,0.2088,0.0679,0.0462,0.1018,0.0142,0.2082)
)

# 找到基因在sce中的行号
set = match(new_model$gene, rownames(scRNA@assays$RNA$counts))

# 初始化score数据框
score = data.frame(cellname = colnames(scRNA@assays$RNA$counts),
                   cellscore = 1:ncol(scRNA@assays$RNA$counts))

valid_idx <- which(!is.na(set))
new_model <- new_model[valid_idx, ]
set <- set[valid_idx]

# 计算每个基因加权值
for (i in 1:nrow(new_model)) {
  x = set[i]                # 基因位置
  y = new_model$weight[i]   # 对应权重
  z = new_model$gene[i]     # 基因名
  score[,z] = scRNA@assays$RNA$counts[x,]*y
}

# 总分
score$cellscore = rowSums(score[,3:(2+nrow(new_model))])
head(score)

# 计算各个细胞群的乳酸化水平
# 将乳酸化分数写入 meta.data，并归一化
scRNA@meta.data$lacScore = (score$cellscore-min(score$cellscore))/(max(score$cellscore)-min(score$cellscore))
scRNA@meta.data$lacScore = score$cellscore


# 细胞分群可视化
# 提取tsne坐标
tsne_df <- as.data.frame(scRNA@reductions$tsne@cell.embeddings)
tsne_df$cluster <- as.factor(scRNA$celltype)

# 提取基因表达数据并与tsne坐标合并
gene_df <- as.data.frame(scRNA@meta.data$lacScore)
merged_df <- cbind(tsne_df, gene_df)
colnames(merged_df) = c('tSNE_1','tSNE_2','cluster','lacScore')
head(merged_df)

#head(merged_df)
#head(scRNA@meta.data)

#aggregate(lacScore ~ celltype, scRNA@meta.data, median)
#aggregate(lacScore ~ sample, scRNA@meta.data, median)
#aggregate(lacScore ~ orig.ident, scRNA@meta.data, mean)

lac_list = aggregate(lacScore ~ celltype, scRNA@meta.data, median)
x = aggregate(lacScore ~ celltype+sample, scRNA@meta.data, median)
y = dcast(x,celltype~sample,value.var = "lacScore")
lac_list = merge(lac_list,y,by = "celltype")
lac_list
write.csv(lac_list,file = "/home/luwei/UVM_output_new/step3.cell.lacScore.median.csv")

lac_list = aggregate(lacScore ~ celltype, scRNA@meta.data, mean)
x = aggregate(lacScore ~ celltype+sample, scRNA@meta.data, mean)
y = dcast(x,celltype~sample,value.var = "lacScore")
lac_list = merge(lac_list,y,by = "celltype")
lac_list
write.csv(lac_list,file = "/home/luwei/UVM_output_new/step3.cell.lacScore.mean.csv")

p5 = ggplot(merged_df, vars = c("tSNE_1", "tSNE_2", 'lacScore'), aes(x = tSNE_1, y = tSNE_2, colour = lacScore)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("lightgrey", "#930404"), limits = c(0, 10), oob = scales::squish) +
  new_scale_color()+
  theme_classic()
p5
ggsave(filename="/home/luwei/UVM_output_new/step3.lacScore.pdf",plot = p5,width = 7,height = 7)

save(scRNA,file = "./outcome/step3_cell_annoation.Rdata")

# 加载包
library(dplyr)
library(ggplot2)
library(ggsignif)

plot_df <- merged_df
colnames(plot_df)[3] <- "celltype"

# 按中位数排序
cell_medians <- plot_df %>%
  group_by(celltype) %>%
  summarize(median_score = median(lacScore, na.rm = TRUE)) %>%
  arrange(desc(median_score))

plot_df$celltype <- factor(plot_df$celltype, levels = cell_medians$celltype)

# Wilcoxon 两两显著性检验
stat_df <- pairwise.wilcox.test(plot_df$lacScore, plot_df$celltype, p.adjust.method = "BH")$p.value
sig_pairs <- as.data.frame(as.table(stat_df)) %>%
  filter(!is.na(Freq) & Freq < 0.05)

# 只保留 Tumor vs 其他组
sig_pairs <- sig_pairs %>%
  filter(Var1 == "Tumor" | Var2 == "Tumor")

# 确保 comparisons 顺序与 factor levels 一致
sig_pairs_list <- apply(sig_pairs[,1:2], 1, function(x){
  factor_levels <- levels(plot_df$celltype)
  if(match(x[1], factor_levels) > match(x[2], factor_levels)) x <- rev(x)
  as.character(x)
}, simplify = FALSE)

sig_labels <- ifelse(sig_pairs$Freq < 0.001, "***",
                     ifelse(sig_pairs$Freq < 0.01, "**",
                            ifelse(sig_pairs$Freq < 0.05, "*", "ns")))

# 计算 Tumor 与其他组上须
whisker_df <- plot_df %>%
  filter(celltype == "Tumor" | celltype %in% c(sig_pairs$Var1, sig_pairs$Var2)) %>%
  group_by(celltype) %>%
  summarize(
    Q3 = quantile(lacScore, 0.75, na.rm = TRUE),
    IQR = IQR(lacScore, na.rm = TRUE),
    upper_whisker = min(max(lacScore, na.rm = TRUE), Q3 + 1.5 * IQR)
  )

y_min <- min(plot_df$lacScore, na.rm = TRUE)
y_max <- max(whisker_df$upper_whisker) * 1.8  # 放大 50%，避免被截断

# 绘制箱线图
p_box <- ggplot(plot_df, aes(x = celltype, y = lacScore, fill = celltype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")) +
  labs(x = "Cell Type", y = "Lactate Score") +
  coord_cartesian(ylim = c(y_min, y_max)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  )

# 添加显著性标注
if(nrow(sig_pairs) > 0){
  whisker_lookup <- setNames(whisker_df$upper_whisker, whisker_df$celltype)
  y_positions <- sapply(1:nrow(sig_pairs), function(i){
    group1 <- as.character(sig_pairs[i, "Var1"])
    group2 <- as.character(sig_pairs[i, "Var2"])
    max(whisker_lookup[group1], whisker_lookup[group2]) * 1.05 + (i-1) * (y_max - y_min)*0.05
  })
  
  p_box <- p_box +
    geom_signif(
      comparisons = sig_pairs_list,
      annotations = sig_labels,
      map_signif_level = TRUE,
      textsize = 6,
      tip_length = 0.02,
      y_position = y_positions
    )
}

# 显示与保存
print(p_box)

ggsave("/home/luwei/UVM_output_new/boxplot_tumor_vs_others.pdf",
       p_box, width = 8, height = 6)

pairwise.wilcox.test(plot_df$lacScore, plot_df$celltype, p.adjust.method = "BH")