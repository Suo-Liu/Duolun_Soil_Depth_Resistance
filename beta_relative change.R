setwd("C:/Users/True/OneDrive/桌面")
library(ieggr)
library(dplyr)
library(vegan)
library(Rmisc)
library(reshape2)

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\data for use\\treatment.csv"

library(ggplot2)
library(tidyverse)
library(gghalves)

Bac.data.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\beta\\relative change\\beta_relative_change", "_", "Bacteria", ".csv")
Fungi.data.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\beta\\relative change\\beta_relative_change", "_", "Fungi", ".csv")
Pro.data.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\beta\\relative change\\beta_relative_change", "_", "Protist", ".csv")

Bac.data <- read.csv(Bac.data.file, row.names = 1, header = T, sep = ",")
Fungi.data <- read.csv(Fungi.data.file, row.names = 1, header = T, sep = ",")
Pro.data <- read.csv(Pro.data.file, row.names = 1, header = T, sep = ",")

sum(Bac.data$block1 == Fungi.data$block1)
sum(Bac.data$treat3 == Fungi.data$treat3)
sum(Bac.data$treat == Fungi.data$treat)
sum(Bac.data$plant1 == Fungi.data$plant1)
sum(Bac.data$layer1 == Fungi.data$layer1)
sum(Bac.data$group == Fungi.data$group)

result <- cbind(
  Bac.data$relative_change, Fungi.data$relative_change, Pro.data$relative_change,
  Bac.data[, c("block1", "treat3", "treat", "plant1", "layer1")]
)
colnames(result) <- c("Bacteria", "Fungi", "Protist", "block", "treat", "method", "plant", "layer")

result$ID <- rownames(result)
library(reshape2)
library(tidyverse)
result <- melt(result, id = c("ID", "layer", "block", "plant", "method", "treat"))

# plot ####
plot.data <- subset(result, method == "Sorenson" & plant == "TS")

prefix.m <- "Sorenson"

plot.data$layer <- ifelse(plot.data$layer %in% c("L1", "L2"), "Topsoil", "Subsoil")
plot.data$layer <- factor(plot.data$layer, levels = c("Topsoil", "Subsoil"))

ggplot(
  data = plot.data,
  aes(x = layer, y = value, fill = layer)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.5) +
  geom_half_boxplot(
    side = "r", errorbar.draw = FALSE, width = 0.25,
    # position = position_dodge(width = 0.76),
    linewidth = 0.2
  ) +
  geom_half_point_panel(
    side = "l", shape = 21, size = 2,
    alpha = 0.7
    # position = position_dodge(width = 0.78)
    , color = "white"
  ) +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("Topsoil", "Subsoil")) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  labs(y = "RC_beta", x = NULL) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
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
    )
  )
ggsave(paste0("beta_relaive change","_", prefix.m, "_top_sub.pdf"), width = 7.17, height = 2.40, units = "in")

# statistics test ####
div.table <- data.frame(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
plot.data <- plot.data[order(plot.data$variable), ]
plot.data <- plot.data[order(plot.data$block), ]
plot.data <- plot.data[order(plot.data$treat), ]
plot.data <- plot.data[order(plot.data$layer), ]
plot.data$ID <- rep(1:2, nrow(plot.data) / 2)
plot.data %<>% pivot_wider(
  names_from = variable,
  values_from = value
)

divindex <- plot.data[, 7:9, drop = F]
treat.used <- plot.data[, -c(7:9)]
treat.used$layer <- ifelse(treat.used$layer == "Topsoil", 0, 1)
divindex <- scale(divindex)
library(car)
library(lme4)
divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ layer + (1 | block) + (1 | treat), data = div)
  presult <- car::Anova(fm, type = 2)
  coefs <- coef(summary(fm))[, "Estimate"] # four coefs
  names(coefs) <- paste0(names(coefs), ".mean")
  SEvalues <- coef(summary(fm))[, "Std. Error"] # standard errors
  names(SEvalues) <- paste0(names(SEvalues), ".se")
  tvalues <- coef(summary(fm))[, "t value"] # t values
  names(tvalues) <- paste0(names(tvalues), ".t")
  chisqP <- c(presult[, 1], presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
  chisqP <- c(presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".P"))
  result <- c(coefs, tvalues, SEvalues, chisqP)
})
colnames(divs1) <- colnames(divindex)
divs1 <- divs1[!grepl("Intercept", rownames(divs1)), ]
divs1 <- as.data.frame(divs1)
divs1$treat <- rownames(divs1)
divs1 <- separate(divs1, col = treat, sep = "\\.", into = c("treat", "type"))
divs1 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)
divs1$method <- "with treat"
div.table <- rbind(div.table, divs1)

divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ layer + (1 | block), data = div)
  presult <- car::Anova(fm, type = 2)
  coefs <- coef(summary(fm))[, "Estimate"] # four coefs
  names(coefs) <- paste0(names(coefs), ".mean")
  SEvalues <- coef(summary(fm))[, "Std. Error"] # standard errors
  names(SEvalues) <- paste0(names(SEvalues), ".se")
  tvalues <- coef(summary(fm))[, "t value"] # t values
  names(tvalues) <- paste0(names(tvalues), ".t")
  chisqP <- c(presult[, 1], presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
  chisqP <- c(presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".P"))
  result <- c(coefs, tvalues, SEvalues, chisqP)
})
colnames(divs1) <- colnames(divindex)
divs1 <- divs1[!grepl("Intercept", rownames(divs1)), ]
divs1 <- as.data.frame(divs1)
divs1$treat <- rownames(divs1)
divs1 <- separate(divs1, col = treat, sep = "\\.", into = c("treat", "type"))
divs1 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)
divs1$method <- "without treat"
div.table <- rbind(div.table, divs1)
write.csv(div.table, paste0("beta", "_", "relative change", "_", prefix.m, "_", "P value.csv"))
