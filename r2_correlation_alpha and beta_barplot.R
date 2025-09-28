setwd("C:\\Users\\True\\OneDrive\\桌面")
alpha.r2.data <- read.csv("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\alpha\\top vs sub\\shannon_LMM_perm test.csv")
beta.r2.data <- read.csv("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\beta\\beta_LMM_perm test.csv")

library(reshape2)
alpha.r2.data <- alpha.r2.data[4, ]
alpha.r2.data <- melt(alpha.r2.data, id.vars = "X")
alpha.r2.data$group <- "Alpha"
beta.r2.data <- beta.r2.data[4, ]
beta.r2.data <- melt(beta.r2.data, id.vars = "X")
beta.r2.data$group <- "Beta"
plot.data <- rbind(alpha.r2.data, beta.r2.data)
plot.data$group <- factor(plot.data$group, levels = c("Alpha", "Beta"))
plot.data$variable <- factor(plot.data$variable, levels = c("Topsoil", "Subsoil"))

library(ggplot2)
library(ggchicklet)
ggplot(plot.data, aes(x = variable, y = value, fill = variable)) +
  geom_chicklet(
    color = NA, # 边缘颜色设置
    width = 0.85,
    radius = unit(9, "pt"),
    position = position_dodge(width = 0.9),
    alpha = 0.8
  ) + # 设置圆角角度
  # 设置图案映射
  xlab("") +
  ylab("R2") +
  facet_wrap(~group, scales = "free_y", ncol = 2) +
  # 设置填充颜色
  scale_fill_manual(
    values = c("#CD2626", "#1874CD"),
    limits = c("Topsoil", "Subsoil")
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.title = element_text(
      size = 15,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      colour = "#000000"
    ),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12, face = "bold")
  )

ggsave("r2_correlation_alpha and beta.pdf", width = 5.26, height = 2.58, units = "in")
