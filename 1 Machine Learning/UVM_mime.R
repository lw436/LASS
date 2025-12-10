# 准备工作
rm(list = ls())
setwd("/Users/luwei/UVM")
library(Mime1)
library(randomForestSRC)
library(survminer)
library(Cairo)
library(snowfall)
load("/Users/luwei/UVM/lac_candidate_genes_184.Rdata")
genelist <- lac_candidate_genes_184
length(unique(genelist))

# 数据预处理
ArrayExpress_df <- read.table("ArrayExpress-dat.txt")
ArrayExpress_t <- as.data.frame(t(ArrayExpress_df))
ArrayExpress_cli <- read.table("ArrayExpress-clin.txt")
ArrayExpress_t$ID <- row.names(ArrayExpress_t)
ArrayExpress_t <- merge(ArrayExpress_t, ArrayExpress_cli, by.x = "ID", by.y = "sample")
n <- ncol(ArrayExpress_t)
Dataset1 <- ArrayExpress_t[, c(
  1,                       # 第1列保留在最前
  (n-1):n,                 # 最后两列移到第2、3列
  2:(n-2)                  # 其余列顺次往后排
)]
colnames(Dataset1)[2:3] <- c("OS.time", "OS")

TCGA_UVM_df <- read.table("TCGA-UVM-dat.txt")
TCGA_UVM_t <- as.data.frame(t(TCGA_UVM_df))
TCGA_UVM_cli <- read.table("TCGA-UVM-clin.txt")
TCGA_UVM_t$ID <- row.names(TCGA_UVM_t)
TCGA_UVM_t <- merge(TCGA_UVM_t, TCGA_UVM_cli, by.x = "ID", by.y = "sample")
n <- ncol(TCGA_UVM_t)
Dataset2 <- TCGA_UVM_t[, c(
  1,                       # 第1列保留在最前
  (n-1):n,                 # 最后两列移到第2、3列
  2:(n-2)                  # 其余列顺次往后排
)]
colnames(Dataset2)[2:3] <- c("OS.time", "OS")

Dataset1 <- Dataset1[Dataset1[, 2] != 0, ]
Dataset2 <- Dataset2[Dataset2[, 2] != 0, ]


list_train_vali_Data <- list(ArrayExpress = Dataset1,
                             TCGA_UVM = Dataset2)

# 构建模型
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$ArrayExpress,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = genelist,
                       mode = 'all',
                       nodesize =5,seed = 5201314)    

write.csv(res$Cindex.res, "result_C_index.csv", row.names = F)

save.image(file = "Mime_model_result.Rdata")
dir.create("./result_new")
setwd("./result_new")

# 绘制所有模型的C-index热图
CairoPDF('1_Heatmap_C_index.pdf', width = 7, height = 12)
cindex_dis_all(
  res,                                  # 模型结果对象
  validate_set = names(list_train_vali_Data)[-1],  # 验证集名称（除第一个之外）
  order = names(list_train_vali_Data),               # 排序顺序
  width = 0.35                         # 热图单元格宽度
)
dev.off()

# 绘制最佳模型的C-index
CairoPDF('2_Best_model_C_index.pdf', width = 8, height = 2.33)
cindex_dis_select(
  res,
  model = "StepCox[forward] + Ridge",   # 选择的最佳模型
  order = names(list_train_vali_Data)     # 排序顺序
)
dev.off()

# 绘制生存曲线
survplot <- vector("list", 2)
for (i in 1:2) {
  # 为每个数据集生成生存曲线
  print(
    survplot[[i]] <- rs_sur(
      res,
      model_name = "StepCox[forward] + Ridge",  # 使用的模型
      dataset = names(list_train_vali_Data)[i],   # 数据集名称
      median.line = "hv",                         # 中位生存线类型
      cutoff = 0.5,                              # 风险分组的截断值
      conf.int = TRUE,                           # 显示置信区间
      xlab = "Day",                             # X轴标签
      pval.coord = c(1000, 0.9)                  # p值位置
    )
  )
}
# 合并生存曲线图
CairoPDF('3_survival_curve.pdf', width = 12, height = 6)
aplot::plot_list(gglist = survplot, ncol = 2)
dev.off()

# 计算各模型在不同时间点的AUC
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["ArrayExpress"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["ArrayExpress"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["ArrayExpress"]],
                             inputmatrix.list = list_train_vali_Data,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")

# 1年AUC绘图
CairoPDF('4a_1year_AUC.pdf', width = 7.2, height = 12)
auc_dis_all(all.auc.1y, dataset = names(list_train_vali_Data), validate_set = names(list_train_vali_Data)[-1],
            order = names(list_train_vali_Data), width = 0.35, year = 1)
dev.off()

# 绘制最佳模型的ROC曲线（1年）
CairoPDF('4b_1year_AUC_by_best_model.pdf', width = 7.08, height = 6)
roc_vis(all.auc.1y, model_name = "StepCox[forward] + Ridge", dataset = names(list_train_vali_Data),
        order = names(list_train_vali_Data), anno_position = c(0.65, 0.55), year = 1)
dev.off()

# 3年AUC绘图
CairoPDF('5a_3year_AUC.pdf', width = 7.2, height = 12)
auc_dis_all(all.auc.3y, dataset = names(list_train_vali_Data), validate_set = names(list_train_vali_Data)[-1],
            order = names(list_train_vali_Data), width = 0.35, year = 3)
dev.off()

# 绘制最佳模型的ROC曲线（3年）
CairoPDF('5b_3year_AUC_by_best_model.pdf', width = 7.08, height = 6)
roc_vis(all.auc.3y, model_name = "StepCox[forward] + Ridge", dataset = names(list_train_vali_Data),
        order = names(list_train_vali_Data), anno_position = c(0.65, 0.55), year = 3)
dev.off()

# 5年AUC绘图
CairoPDF('6a_5year_AUC.pdf', width = 7.2, height = 12)
auc_dis_all(all.auc.5y, dataset = names(list_train_vali_Data), validate_set = names(list_train_vali_Data)[-1],
            order = names(list_train_vali_Data), width = 0.35, year = 5)
dev.off()

# 5年最佳模型ROC曲线
CairoPDF('6b_5year_AUC_by_best_model.pdf', width = 7.08, height = 6)
roc_vis(all.auc.5y, model_name = "StepCox[forward] + Ridge", dataset = names(list_train_vali_Data),
        order = names(list_train_vali_Data), anno_position = c(0.65, 0.55), year = 5)
dev.off()


# 绘制最佳模型在不同时间点的AUC比较
pdf('7_123year_AUC_by_best_model.pdf', width = 12, height = 2)
auc_dis_select(list(all.auc.1y, all.auc.2y, all.auc.3y),  # 两个时间点的结果
               model_name = "StepCox[forward] + Lasso",    # 模型名称
               dataset = names(list_train_vali_Data),       # 数据集
               order = names(list_train_vali_Data),         # 排序
               year = c(1, 2, 3))                            # 三个时间点
dev.off()

# 计算单因素Cox回归结果
unicox.rs.res <- cal_unicox_ml_res(
  res.by.ML.Dev.Prog.Sig = res,                   # 模型结果
  optimal.model = "StepCox[forward] + Lasso",               # 选择的最佳模型
  type = 'categorical'                            # 分类变量（高/低风险组）
)

# 进行meta分析
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)

# 调整数据框列顺序
metamodel$data <- metamodel$data[, c(1:6, 8, 9, 7)]

# 可视化meta分析结果
png('8_Cox_regression_by_best_model.png', width = 11, height = 4.5, units = 'in', res = 300)
meta_unicox_vis(
  metamodel,                                     # meta分析结果
  dataset = names(list_train_vali_Data)          # 数据集名称
)
dev.off()