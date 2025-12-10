library(dplyr)
library(ggplot2)

# =======================================================
# 1. 读取数据
# =======================================================
go_ready <- read.csv("C:/Users/LuWei/OneDrive/桌面/GO富集分析.csv")

# 改列名
go_ready <- go_ready %>%
  rename(
    Description = Name,
    P_value = p.value,
    Gene_Count = Hit.Count.in.Query.List,
    Group = Cluster
  ) %>%
  mutate(
    Group = paste0("cluster ", Group)  # cluster 1 / cluster 2 / cluster 3
  )

# 计算 -log10(p)
go_ready <- go_ready %>%
  mutate(NegLogP = -log10(P_value)) %>%
  arrange(Group, desc(NegLogP)) %>%
  mutate(Unique_ID = paste(Group, Description, sep = " - "))

go_ready$Plot_ID <- factor(go_ready$Unique_ID, levels = rev(go_ready$Unique_ID))
y_labels <- setNames(as.character(go_ready$Description), go_ready$Plot_ID)


# =======================================================
# 2. 配色方案 + 二轴转换
# =======================================================
sci_colors <- c(
  "cluster 1" = "#00BA38",
  "cluster 2" = "#F8766D",
  "cluster 3" = "#619CFF"
)

max_neglogp <- max(go_ready$NegLogP)
max_gene_count <- max(go_ready$Gene_Count)
conversion_factor <- max_neglogp / max_gene_count
grid_lines <- seq(1.5, length(unique(go_ready$Plot_ID)) - 0.5, 1)


# =======================================================
# 3. 绘图
# =======================================================
p <- ggplot(go_ready, aes(y = Plot_ID)) +
  
  geom_hline(yintercept = grid_lines, color = "gray92", size = 0.5) +
  
  # 柱状图：P-value
  geom_bar(
    aes(x = NegLogP, fill = Group),
    stat = "identity",
    width = 0.65,
    color = "black",
    size = 0.3
  ) +
  
  # shape 图例（P-value）
  geom_point(
    aes(shape = "P-value Metric"),
    x = -Inf,
    y = go_ready$Plot_ID[1],
    fill = sci_colors[1],
    color = "black",
    size = 4,
    data = go_ready[1,]
  ) +
  
  # Gene count 折线
  geom_path(
    aes(x = Gene_Count * conversion_factor, group = 1, color = "Gene Count"),
    size = 0.6
  ) +
  
  # Gene count 点
  geom_point(
    aes(x = Gene_Count * conversion_factor, color = "Gene Count"),
    size = 2.5,
    shape = 15
  ) +
  
  scale_x_continuous(
    name = expression(-log[10](italic(p))),
    expand = expansion(mult = c(0, 0.15)),
    sec.axis = sec_axis(~ . / conversion_factor, name = "Gene Count")
  ) +
  
  scale_y_discrete(labels = y_labels) +
  
  scale_fill_manual(
    name = NULL,
    values = sci_colors,
    guide = "none"   # 不显示填充色图例
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c("Gene Count" = "black"),
    labels = c("Gene Count" = "Gene Count (Line)")
  ) +
  
  scale_shape_manual(
    name = NULL,
    values = c("P-value Metric" = 22),
    labels = c("P-value Metric" = expression(-log[10](italic(p))))
  ) +
  
  guides(
    shape = guide_legend(
      order = 1,
      keywidth = unit(1.0, "cm"),
      keyheight = unit(0.5, "cm"),
      override.aes = list(
        fill = "#00BA38",
        size = 6,
        linetype = 0,
        color = "black"
      )
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(
        shape = 15, size = 3, linetype = 1, fill = NA
      )
    )
  ) +
  
  theme_classic() +
  theme(
    axis.line.x.bottom = element_line(color = "black", size = 0.6),
    axis.line.x.top = element_line(color = "black", size = 0.6),
    axis.line.y = element_line(color = "black", size = 0.6),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length.x = unit(-0.2, "cm"),
    axis.text.x.bottom = element_text(color = "black", margin = margin(t = 8)),
    axis.text.x.top = element_text(color = "black", margin = margin(b = 8)),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.box = "horizontal",
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.x = unit(0.5, "cm")
  ) +
  
  labs(title = "GO Enrichment Analysis")

print(p)

ggsave("GO_Enrichment_Analysis_Plot.pdf", p, width = 8, height = 6)
