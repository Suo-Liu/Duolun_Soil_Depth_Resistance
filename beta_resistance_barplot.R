setwd("C:\\Users\\True\\OneDrive\\桌面")

# plot ####
library(ggplot2)
prefix.m = "mixed"
resistance.data = read.csv(paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\beta\\resistance\\",prefix.m,"_distance for resistance.csv"),
                           header = T, row.names = 1, sep = ",", stringsAsFactors = F)
resistance.data$layer1 = ifelse(resistance.data$layer %in% c("L1","L2"), "Topsoil", "Subsoil")
resistance.data$resistance = 1 - resistance.data$dis
resistance.data = subset(resistance.data, treat == "Bray")

## all treatments & all steppes ####
library(Rmisc)
plot.data = resistance.data
plot.data = subset(plot.data, treat3 %in% c("RP","EP","W","WRP","WEP"))
plot.data = summarySE(data = plot.data, measurevar = "resistance", groupvars = c("layer1"))
plot.data$layer1 <- factor(plot.data$layer1, levels = c("Topsoil","Subsoil"))

ggplot(data = plot.data,
       aes(x = layer1, y = resistance, color = layer1, fill = layer1)) +
  geom_bar(alpha = 0.8, stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = resistance - se, ymax = resistance + se),
                color = "black", width = 0.5,
                position = position_dodge(.9), lwd = 0.7
  ) +
  scale_fill_manual(values = c("#1874CD", "#CD2626"), limits = c("Topsoil","Subsoil")) +
  scale_color_manual(values = c("#1874CD", "#CD2626"), limits =  c("Topsoil","Subsoil")) +
  theme_bw() +
  # remove the grid lines
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 取消刻度线：axis.ticks.x, axis.ticks.y
    # axis.ticks = element_blank(),
    plot.title = element_text(
      size = 20,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 13,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 12,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    )
    # ,
    # legend.title = element_text(
    #   size = 12,
    #   face = "bold"
    # ),
    # legend.text = element_text(
    #   size = 10,
    #   face = "bold"
    # )
  ) +
  labs(y="Resistance",x=NULL) +
  # facet_wrap(~treat, scales = "free_y", ncol = 3) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"))

ggsave(paste0(prefix.m,"_","all treatments","_","all steppes",".pdf"), width = 4.15, height = 4.15, units = "in")

## each treatment & all steppes ####
library(Rmisc)
plot.data = resistance.data
plot.data = subset(plot.data, treat3 %in% c("RP","EP","W","WRP","WEP"))
plot.data = summarySE(data = plot.data, measurevar = "resistance", groupvars = c("layer1","treat3"))
plot.data$layer1 <- factor(plot.data$layer1, levels = c("Topsoil","Subsoil"))

ggplot(data = plot.data,
       aes(x = treat3, y = resistance, color = layer1, fill = layer1)) +
  geom_bar(alpha = 0.8, stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = resistance - se, ymax = resistance + se),
                color = "black", width = 0.5,
                position = position_dodge(.9), lwd = 0.7
  ) +
  scale_fill_manual(values = c("#1874CD", "#CD2626"), limits = c("Topsoil","Subsoil")) +
  scale_color_manual(values = c("#1874CD", "#CD2626"), limits =  c("Topsoil","Subsoil")) +
  theme_bw() +
  # remove the grid lines
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 取消刻度线：axis.ticks.x, axis.ticks.y
    # axis.ticks = element_blank(),
    plot.title = element_text(
      size = 20,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 13,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 12,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    )
    # ,
    # legend.title = element_text(
    #   size = 12,
    #   face = "bold"
    # ),
    # legend.text = element_text(
    #   size = 10,
    #   face = "bold"
    # )
  ) +
  labs(y="Resistance",x=NULL) +
  # facet_wrap(~treat, scales = "free_y", ncol = 3) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"))

ggsave(paste0(prefix.m,"_","each treatment","_","all steppes",".pdf"), width = 7.54, height = 5.69, units = "in")
