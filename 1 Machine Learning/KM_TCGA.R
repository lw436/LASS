rm(list=ls())

library(survival)
library(survminer)

setwd("/Users/luwei/UVM/survival")
TCGA_group <- read.table("/Users/luwei/UVM/TCGA_group.txt", header = T)
TCGA_surv <- read.table("/Users/luwei/UVM/TCGA_clin.txt", header = T, sep = "\t")
TCGA_group$sample <- gsub("-", "\\.", TCGA_group$sample)

TCGA_mean <- mean(TCGA_group$RiskScore)
TCGA_group$group <- ifelse(TCGA_group$RiskScore > TCGA_mean, "High Risk", "Low Risk")

colnames(TCGA_surv)[1] <- "sample"
TCGA_group <- TCGA_group[,-2]
TCGA_clin <- merge(TCGA_surv, TCGA_group, by = "sample")

# 总体生存分析 OS-----
# 构建 COX 模型计算 HR
cox_OS <- coxph(Surv(OS.time, OS) ~ group, data = TCGA_clin)
summary(cox_OS)

# 提取HR和置信区间
HR <- round(summary(cox_OS)$coefficients[,"exp(coef)"], 2)
HR.confint <- round(summary(cox_OS)$conf.int[,c("lower .95", "upper .95")], 2)
pvalue <- signif(summary(cox_OS)$coefficients[,"Pr(>|z|)"], 3)

HR_label <- paste0("HR = ", HR, " (95% CI: ", HR.confint[1], "-", HR.confint[2], ")\nP = ", pvalue)
HR_label

fit_OS <- survfit(Surv(OS.time, OS) ~ group, data = TCGA_clin)

p_OS <- ggsurvplot(
  fit_OS,
  data = TCGA_clin,
  pval = FALSE,
  conf.int = TRUE,
  risk.table = FALSE,   # 不显示风险表
  palette = c("#E64B35FF", "#4DBBD5FF"),
  ggtheme = theme_bw(),
  legend.title = "Group",
  legend.labs = c("High Risk", "Low Risk"),
  surv.median.line = "hv",
  #title = "TCGA-UVM Overall Survival",
  xlab = "Overall Survival (day)",
  ylab = "Survival probability"
)

# 添加HR标签（放在左下角）& 标题居中
p_OS$plot <- p_OS$plot +
  annotate("text",
           x = max(TCGA_clin$OS.time) * 0.005,   # ← 控制横向位置
           y = 0.2,                             # ← 控制纵向位置
           label = HR_label,
           size = 3.5,
           hjust = 0,
           color = "black") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # 居中标题
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "top"
  )

pdf("TCGA_UVM_OS_KM.pdf", width = 6, height = 5)
print(p_OS, newpage = FALSE)
dev.off()

# 疾病特异生存分析 DSS-----
# 构建 COX 模型计算 HR
cox_DSS <- coxph(Surv(DSS.time, DSS) ~ group, data = TCGA_clin)
summary(cox_DSS)

# 提取HR和置信区间
HR <- round(summary(cox_DSS)$coefficients[,"exp(coef)"], 2)
HR.confint <- round(summary(cox_DSS)$conf.int[,c("lower .95", "upper .95")], 2)
pvalue <- signif(summary(cox_DSS)$coefficients[,"Pr(>|z|)"], 3)

HR_label <- paste0("HR = ", HR, " (95% CI: ", HR.confint[1], "-", HR.confint[2], ")\nP = ", pvalue)
HR_label

fit_DSS <- survfit(Surv(DSS.time, DSS) ~ group, data = TCGA_clin)

p_DSS <- ggsurvplot(
  fit_DSS,
  data = TCGA_clin,
  pval = FALSE,
  conf.int = TRUE,
  risk.table = FALSE,   # 不显示风险表
  palette = c("#E64B35FF", "#4DBBD5FF"),
  ggtheme = theme_bw(),
  legend.title = "Group",
  legend.labs = c("High Risk", "Low Risk"),
  surv.median.line = "hv",
  #title = "TCGA-UVM Overall Survival",
  xlab = "Disease Specific Survival (day)",
  ylab = "Survival probability"
)

# 添加HR标签（放在左下角）& 标题居中
p_DSS$plot <- p_DSS$plot +
  annotate("text",
           x = max(TCGA_clin$DSS.time) * 0.05,   # ← 控制横向位置
           y = 0.2,                             # ← 控制纵向位置
           label = HR_label,
           size = 3.5,
           hjust = 0,
           color = "black") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # 居中标题
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "top"
  )

pdf("TCGA_UVM_DSS_KM.pdf", width = 6, height = 5)
print(p_DSS, newpage = FALSE)
dev.off()

# 进展自由生存分析 PFI-----
# 构建 COX 模型计算 HR
cox_PFI <- coxph(Surv(PFI.time, PFI) ~ group, data = TCGA_clin)
summary(cox_PFI)

# 提取HR和置信区间
HR <- round(summary(cox_PFI)$coefficients[,"exp(coef)"], 2)
HR.confint <- round(summary(cox_PFI)$conf.int[,c("lower .95", "upper .95")], 2)
pvalue <- signif(summary(cox_PFI)$coefficients[,"Pr(>|z|)"], 3)

HR_label <- paste0("HR = ", HR, " (95% CI: ", HR.confint[1], "-", HR.confint[2], ")\nP = ", pvalue)
HR_label

fit_PFI <- survfit(Surv(PFI.time, PFI) ~ group, data = TCGA_clin)

p_PFI <- ggsurvplot(
  fit_PFI,
  data = TCGA_clin,
  pval = FALSE,
  conf.int = TRUE,
  risk.table = FALSE,   # 不显示风险表
  palette = c("#E64B35FF", "#4DBBD5FF"),
  ggtheme = theme_bw(),
  legend.title = "Group",
  legend.labs = c("High Risk", "Low Risk"),
  surv.median.line = "hv",
  #title = "TCGA-UVM Overall Survival",
  xlab = "Progression Free Interval (day)",
  ylab = "Survival probability"
)

# 添加HR标签（放在左下角）& 标题居中
p_PFI$plot <- p_PFI$plot +
  annotate("text",
           x = 1,   # ← 控制横向位置
           y = 0.2,                             # ← 控制纵向位置
           label = HR_label,
           size = 3.5,
           hjust = 0,
           color = "black") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # 居中标题
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "top"
  )

pdf("TCGA_UVM_PFI_KM.pdf", width = 6, height = 5)
print(p_PFI, newpage = FALSE)
dev.off()
