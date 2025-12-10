rm(list=ls())

library(survival)
library(survminer)

setwd("/Users/luwei/UVM/survival")
ArrayExpress_group <- read.table("/Users/luwei/UVM/ArrayExpress_group.txt", header = T)
ArrayExpress_surv <- read.table("/Users/luwei/UVM/ArrayExpress-clin.txt", header = T)
ArrayExpress_mean <- mean(ArrayExpress_group$RiskScore)
ArrayExpress_group$group <- ifelse(ArrayExpress_group$RiskScore > ArrayExpress_mean, "High Risk", "Low Risk")

ArrayExpress_group <- ArrayExpress_group[,-2]
ArrayExpress_clin <- merge(ArrayExpress_surv, ArrayExpress_group, by = "sample")

# 总体生存分析 OS-----
# 构建 COX 模型计算 HR
cox_OS <- coxph(Surv(time, status) ~ group, data = ArrayExpress_clin)
summary(cox_OS)

# 提取HR和置信区间
HR <- round(summary(cox_OS)$coefficients[,"exp(coef)"], 2)
HR.confint <- round(summary(cox_OS)$conf.int[,c("lower .95", "upper .95")], 2)
pvalue <- signif(summary(cox_OS)$coefficients[,"Pr(>|z|)"], 3)

HR_label <- paste0("HR = ", HR, " (95% CI: ", HR.confint[1], "-", HR.confint[2], ")\nP = ", pvalue)
HR_label

fit_OS <- survfit(Surv(time, status) ~ group, data = ArrayExpress_clin)

p_OS <- ggsurvplot(
  fit_OS,
  data = ArrayExpress_clin,
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
           x = max(ArrayExpress_clin$time) * 0.05,   # ← 控制横向位置
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

pdf("ArrayExpress_UVM_OS_KM.pdf", width = 6, height = 5)
print(p_OS, newpage = FALSE)
dev.off()