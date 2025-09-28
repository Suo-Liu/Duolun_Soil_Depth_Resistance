setwd("C:\\Users\\True\\OneDrive\\桌面")
com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
prefix.m <- "Bacteria"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\ITS\\all samples\\Unoise\\zotutab220_max_2_resample_10000.txt"
prefix.m <- "Fungi"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"
prefix.m <- "Protist"

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

comm <- t(read.table(com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))

name.row <- rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] <- "YD3L3"
rownames(comm) <- name.row

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
treat <- subset(treat, plant.type == "TS")
library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, comm = comm))
dim(comm)
comm <- sampc$comm
treat <- sampc$treat
comm <- comm[, colSums(comm) > 0]
dim(comm)

library(spaa)
niche_width <- niche.width(comm, method = "shannon")

# 计算样本平均niche width ####
# weighted by abundance
niche_width_sample_weighted <-
  sapply(1:nrow(comm), function(i) {
    sum(comm[i, ] * niche_width) / sum(comm[i, ])
  })
names(niche_width_sample_weighted) <- rownames(comm)
niche_width_sample_weighted <- as.data.frame(niche_width_sample_weighted)
niche_width_sample_weighted$sample <- rownames(niche_width_sample_weighted)
colnames(niche_width_sample_weighted) <- c("niche_width", "sample")
niche_width_sample_weighted$Layer <- treat$Layer[match(niche_width_sample_weighted$sample, rownames(treat))]
niche_width_sample_weighted$treat <- treat$combined_treat1[match(niche_width_sample_weighted$sample, rownames(treat))]
# unweighted by abundance
niche_width_sample_unweighted <-
  sapply(1:nrow(comm), function(i) {
    sum(niche_width[comm[i, ] > 0]) / length(comm[i, comm[i, ] > 0])
  })
names(niche_width_sample_unweighted) <- rownames(comm)
niche_width_sample_unweighted <- as.data.frame(niche_width_sample_unweighted)
niche_width_sample_unweighted$sample <- rownames(niche_width_sample_unweighted)
colnames(niche_width_sample_unweighted) <- c("niche_width", "sample")
niche_width_sample_unweighted$Layer <- treat$Layer[match(niche_width_sample_unweighted$sample, rownames(treat))]
niche_width_sample_unweighted$treat <- treat$combined_treat1[match(niche_width_sample_unweighted$sample, rownames(treat))]

write.csv(niche_width_sample_weighted, paste0(prefix.m,".","niche_width_sample_weighted.csv"))
write.csv(niche_width_sample_unweighted, paste0(prefix.m,".","niche_width_sample_unweighted.csv"))

# plot ####
library(ggplot2)
library(tidyverse)
library(gghalves)
plot.data.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\niche\\top vs sub\\shannon\\niche_width_sample.csv"
plot.data <- read.csv(plot.data.file, header = T)
plot.data <- subset(plot.data, method == "unweighted")
plot.data$Layer <- ifelse(plot.data$Layer %in% c("L1", "L2"), "Topsoil", "Subsoil")
plot.data$Layer <- factor(plot.data$Layer, levels = c("Topsoil", "Subsoil"))

ggplot(
  data = plot.data,
  aes(x = Layer, y = niche_width, fill = Layer)
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
  labs(y = "Niche width", x = NULL) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 3) +
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

ggsave(paste0("shannon", "_", "unweighted", "_", "top vs sub", "_", "niche width.pdf"), width = 6.02, height = 2.17, units = "in")

# statistics test ####
plot.data$block = treat$block[match(plot.data$sample, rownames(treat))]
plot.data$Layer <- ifelse(plot.data$Layer == "Topsoil", 0, 1)
div.table <- data.frame(stringsAsFactors = F)
divindex = plot.data[,1,drop = F]
treat.used = plot.data[, -1]
divindex <- scale(divindex)
library(car)
library(lme4)
divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ Layer + (1 | block) + (1 | treat), data = div)
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
library(tidyverse)
divs1 <- separate(divs1, col = treat, sep = "\\.", into = c("treat", "type"))
divs1 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)
divs1$method <- "with treat"
div.table <- rbind(div.table, divs1)

divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ Layer + (1 | block), data = div)
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
write.csv(div.table, paste0(prefix.m, "_","shannon","_", "unweighted", "_", "niche width", "_", "P value.csv"))
