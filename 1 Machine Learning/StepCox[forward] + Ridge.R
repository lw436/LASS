library(ggplot2)
library(survival)
library(survminer)

# 提取 glmnet 模型对象
glmnet_fit <- res$ml.res$`StepCox[forward] + Ridge`$glmnet.fit

# 最优 lambda
best_lambda <- res$ml.res$`StepCox[forward] + Ridge`$lambda.min

# 提取对应的系数 (包括截距和基因)
coef_mat <- as.matrix(coef(glmnet_fit, s = best_lambda))

# 去掉截距（如果有的话），只保留基因和对应系数
coef_df <- data.frame(
  gene = rownames(coef_mat),
  coef = coef_mat[,1],
  stringsAsFactors = FALSE
)
coef_df <- coef_df[coef_df$coef != 0, ]  # 只保留非零系数
head(coef_df)

# 打印基因数量
gene_used <- rownames(coef_df)
num_genes <- nrow(coef_df)
cat("模型中使用的基因数量:", num_genes, "\n")
write.table(gene_used, file = "./result/used_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 打印公式
formula_str <- paste0(round(coef_df$coef, 4), "*", coef_df$gene, collapse = " + ")
cat("RiskScore = ", formula_str, "\n")

out_dir <- "./Additional_Plots"
dir.create(out_dir, showWarnings = FALSE)

# =========================
# 1. Ridge Coefficient Path & CV Deviance
# =========================
ridge_model <- res$ml.res[["StepCox[forward] + Ridge"]]
best_lambda <- ridge_model$lambda.min

pdf(file.path(out_dir, "Ridge_Coefficient_Path.pdf"), width = 6, height = 6)
plot(ridge_model$glmnet.fit, xvar = "lambda", label = TRUE, main = "Coefficient Path")
dev.off()

pdf(file.path(out_dir, "Ridge_CV_Deviance.pdf"), width = 6, height = 6)
plot(ridge_model, xvar = "lambda", label = TRUE, main = "CV Deviance")
abline(v = log(best_lambda), col = "red", lty = 2, lwd = 2)
dev.off()
# =========================
# 2. 绘制 ROC 曲线
# =========================
plot_multi_year_roc <- function(dataset_name, years = c(1,3,5),
                                model_name = "StepCox[forward] + Ridge",
                                colors = c("#124f7b","#e89b38","#46a5c7"),
                                out_dir = ".") {
  out_file <- file.path(out_dir, paste0(dataset_name,"_MultiYear_ROC.pdf"))
  
  # 确保无论如何都会关闭设备
  pdf(out_file, width=6, height=6)
  on.exit({
    try(dev.off(), silent = TRUE)
  }, add = TRUE)
  
  # 空坐标系
  plot(NA, NA, xlim=c(0,1), ylim=c(0,1),
       xlab="False Positive Rate", ylab="True Positive Rate",
       main=paste0(dataset_name," Multi-Year ROC"))
  
  first <- TRUE
  auc_text <- c()
  colors <- rep(colors, length.out = length(years))
  plotted_any <- FALSE
  
  for(i in seq_along(years)){
    year <- years[i]
    auc_list_name <- paste0("all.auc.", year, "y")
    
    if(!exists(auc_list_name, inherits = TRUE)) next
    auc_list <- get(auc_list_name, inherits = TRUE)
    
    if(!(model_name %in% names(auc_list))) next
    model_entry <- auc_list[[model_name]]
    if(!(dataset_name %in% names(model_entry))) next
    
    df <- model_entry[[dataset_name]]
    if(nrow(df)==0) next
    
    # 清理 FP/TP
    df <- df[is.finite(df$FP) & is.finite(df$TP), ]
    if(nrow(df)==0) next
    
    # 按 FP 排序
    df <- df[order(df$FP), ]
    FP <- df$FP
    TP <- df$TP
    
    # 单调化，保证 ROC 不往回走
    TP <- cummax(TP)
    FP <- cummax(FP)
    
    # 绘制阶梯线
    lines(c(0, FP, 1), c(0, TP, 1), type="s", col=colors[i], lwd=2)
    
    # 计算 AUC（梯形法）
    auc_val <- sum(diff(c(0, FP, 1)) * (c(0, TP, 1)[-length(TP)+0:1] + c(0, TP, 1)[-1])/2)
    auc_text <- c(auc_text, paste0("AUC@", year,"yr: ", round(auc_val,3)))
    
    plotted_any <- TRUE
  }
  
  # 画对角线和图例
  if(plotted_any){
    abline(0,1,lty=2,col="gray")
    legend("bottomright", legend = auc_text, col=colors[1:length(auc_text)],
           lty=1, lwd=2, bty="n")
  } else {
    text(0.5, 0.5, "No valid ROC data", cex=1)
  }
  
  invisible(out_file)
}

# 调用示例
plot_multi_year_roc("ArrayExpress", years=c(1,3,5))
plot_multi_year_roc("TCGA_UVM", years=c(1,3,5))