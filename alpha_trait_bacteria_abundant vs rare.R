setwd("C:/Users/True/OneDrive/桌面")
# abundant and rare species ####
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
prefix.m <- "Bacteria"

comm <- t(read.table(com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))
name.row <- rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] <- "YD3L3"
rownames(comm) <- name.row

pre.treat <- read.csv(treat.file, header = T, row.names = 1)

pre.treat <- subset(pre.treat, plant.type == "TS")
prefix.m <- paste0(prefix.m, "_", "TS")

treat <- pre.treat
library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, comm = comm), silent = T)
dim(comm)
comm <- sampc$comm
treat <- sampc$treat
comm <- comm[, colSums(comm) > 0]
dim(comm)

# 统计物种出现频率
library(dplyr)
fre.comm <- comm
fre.comm <- as.data.frame(fre.comm)
fre.comm$Layer <- treat$Layer
fre.comm$Layer <- ifelse(fre.comm$Layer %in% c("L1", "L2"),
  "Topsoils", "Subsoils"
)

fre.comm %<>% group_by(Layer) %>%
  summarise_all(sum)
fre.comm <- as.data.frame(fre.comm)
rownames(fre.comm) <- fre.comm[, 1]
fre.comm <- fre.comm[, -1]
# 识别交集中的物种和独有物种
inter.otu <- names(fre.comm)[which(colSums(fre.comm > 0) == nrow(fre.comm))]
# 识别丰富物种和稀有物种，定义源于Dai et al., NC, 2022，但是丰富物种的定义少了出现频率的限制，因为过多限制会导致最后识别的丰富物种太少
abun.otu <- which(colSums(comm / rowSums(comm) >= 0.001) > 0.5 & (colSums(comm > 0) > (0.5 * nrow(comm))))
# abun.otu = which(colSums(comm/rowSums(comm) >= 0.001) > 0.5)
abun.otu <- abun.otu[which(names(abun.otu) %in% inter.otu)]

rare.otu <- which(colSums(comm / rowSums(comm) < 0.001) == nrow(comm))

abun.comm <- comm[, colnames(comm) %in% names(abun.otu)]
rare.comm <- comm[, colnames(comm) %in% names(rare.otu)]

# 因为要用物种多样性和丰度数据预测抗性，而抗性是通过处理与对照间的群落差异计算得到的，
# 因此要match对照的物种丰度（应该不需要对照和处理的平均丰度？）
mean.abun.table <- data.frame(stringsAsFactors = F)
mean.rare.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(
    ct.treat = ct.treat, abun.comm = abun.comm,
    rare.comm = rare.comm
  ), , silent = T)
  ct.abun.comm <- sampc$abun.comm
  ct.rare.comm <- sampc$rare.comm
  ct.treat <- sampc$ct.treat

  mean.abun.comm <- ct.abun.comm
  mean.rare.comm <- ct.rare.comm

  mean.abun.comm <- as.data.frame(mean.abun.comm)
  mean.abun.comm$block <- ct.treat$block
  mean.abun.comm$Layer <- ct.treat$Layer
  mean.abun.comm$treat <- ct.treat$combined_treat1

  mean.abun.table <- rbind(mean.abun.table, mean.abun.comm)

  mean.rare.comm <- as.data.frame(mean.rare.comm)
  mean.rare.comm$block <- ct.treat$block
  mean.rare.comm$Layer <- ct.treat$Layer
  mean.rare.comm$treat <- ct.treat$combined_treat1
  mean.rare.table <- rbind(mean.rare.table, mean.rare.comm)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(
    ct.treat = ct.treat, abun.comm = abun.comm,
    rare.comm = rare.comm
  ), , silent = T)
  ct.abun.comm <- sampc$abun.comm
  ct.rare.comm <- sampc$rare.comm
  ct.treat <- sampc$ct.treat

  mean.abun.comm <- ct.abun.comm
  mean.rare.comm <- ct.rare.comm

  mean.abun.comm <- as.data.frame(mean.abun.comm)
  mean.abun.comm$block <- ct.treat$block
  mean.abun.comm$Layer <- ct.treat$Layer
  mean.abun.comm$treat <- ct.treat$combined_treat1
  mean.abun.comm$treat <- paste0(i, ".", mean.abun.comm$treat)
  mean.abun.table <- rbind(mean.abun.table, mean.abun.comm)

  mean.rare.comm <- as.data.frame(mean.rare.comm)
  mean.rare.comm$block <- ct.treat$block
  mean.rare.comm$Layer <- ct.treat$Layer
  mean.rare.comm$treat <- ct.treat$combined_treat1
  mean.rare.comm$treat <- paste0(i, ".", mean.rare.comm$treat)
  mean.rare.table <- rbind(mean.rare.table, mean.rare.comm)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(
    ct.treat = ct.treat, abun.comm = abun.comm,
    rare.comm = rare.comm
  ), silent = T)
  ct.abun.comm <- sampc$abun.comm
  ct.rare.comm <- sampc$rare.comm
  ct.treat <- sampc$ct.treat

  mean.abun.comm <- ct.abun.comm
  mean.rare.comm <- ct.rare.comm

  mean.abun.comm <- as.data.frame(mean.abun.comm)
  mean.abun.comm$block <- ct.treat$block
  mean.abun.comm$Layer <- ct.treat$Layer
  mean.abun.comm$treat <- ct.treat$combined_treat1
  mean.abun.comm$treat <- paste0(i, ".", "W")
  mean.abun.table <- rbind(mean.abun.table, mean.abun.comm)

  mean.rare.comm <- as.data.frame(mean.rare.comm)
  mean.rare.comm$block <- ct.treat$block
  mean.rare.comm$Layer <- ct.treat$Layer
  mean.rare.comm$treat <- ct.treat$combined_treat1
  mean.rare.comm$treat <- paste0(i, ".", "W")
  mean.rare.table <- rbind(mean.rare.table, mean.rare.comm)
}

mean.abun.table <- mean.abun.table[order(mean.abun.table$treat), ]
mean.abun.table <- mean.abun.table[order(mean.abun.table$block), ]
mean.abun.table <- mean.abun.table[order(mean.abun.table$Layer), ]

mean.rare.table <- mean.rare.table[order(mean.rare.table$treat), ]
mean.rare.table <- mean.rare.table[order(mean.rare.table$block), ]
mean.rare.table <- mean.rare.table[order(mean.rare.table$Layer), ]

mean.abun.table <- mean.abun.table[, -c((ncol(mean.abun.table) - 2):ncol(mean.abun.table))]
mean.rare.table <- mean.rare.table[, -c((ncol(mean.rare.table) - 2):ncol(mean.rare.table))]

# bacterial trait abundance ####
trait.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\predicted_traits.csv"

trait.data <- read.csv(trait.file, header = T, row.names = 1)
trait.data$cell_diameter <- 10^trait.data$cell_diameter
trait.data$cell_length <- 10^trait.data$cell_length
trait.data$doubling_h <- 10^trait.data$doubling_h
trait.data <- trait.data[, 1:11]

mean.abun.table <- t(mean.abun.table)
mean.rare.table <- t(mean.rare.table)

spc <- match.name(rn.list = list(mean.abun.table = mean.abun.table, trait.data = trait.data), silent = T)
mean.abun.table <- spc$mean.abun.table
abun.trait.data <- spc$trait.data

spc <- match.name(rn.list = list(mean.rare.table = mean.rare.table, trait.data = trait.data), silent = T)
mean.rare.table <- spc$mean.rare.table
rare.trait.data <- spc$trait.data

abun.list <- list()
for (i in 1:ncol(abun.trait.data)) {
  for (j in 1:ncol(mean.abun.table)) {
    mean.abun.table[, j] <- mean.abun.table[, j] * abun.trait.data[, i]
  }
  abun.list[[i]] <- mean.abun.table
  names(abun.list)[i] <- colnames(abun.trait.data)[i]
}

rare.list <- list()
for (i in 1:ncol(rare.trait.data)) {
  for (j in 1:ncol(mean.rare.table)) {
    mean.rare.table[, j] <- mean.rare.table[, j] * rare.trait.data[, i]
  }
  rare.list[[i]] <- mean.rare.table
  names(rare.list)[i] <- colnames(rare.trait.data)[i]
}

trait_names <- names(abun.trait.data)
abun.adjusted_mats <- list()
for (trait in trait_names) {
  adjusted_mat <- mean.abun.table / abun.trait.data[[trait]]
  adjusted_mat[is.infinite(adjusted_mat)] <- NA

  abun.adjusted_mats[[trait]] <- t(adjusted_mat)
}
abun.weighted.result <- data.frame(row.names = colnames(mean.abun.table))
for (trait in trait_names) {
  w.sample.aver <- colSums(mean.abun.table) / colSums(mean.abun.table / abun.trait.data[[trait]], na.rm = TRUE)
  abun.weighted.result[[trait]] <- w.sample.aver
}

trait_names <- names(rare.trait.data)
rare.adjusted_mats <- list()
for (trait in trait_names) {
  adjusted_mat <- mean.rare.table / rare.trait.data[[trait]]
  adjusted_mat[is.infinite(adjusted_mat)] <- NA

  rare.adjusted_mats[[trait]] <- t(adjusted_mat)
}
rare.weighted.result <- data.frame(row.names = colnames(mean.rare.table))
for (trait in trait_names) {
  w.sample.aver <- colSums(mean.rare.table) / colSums(mean.rare.table / rare.trait.data[[trait]], na.rm = TRUE)
  rare.weighted.result[[trait]] <- w.sample.aver
}

abun.unweighted.result <- data.frame(row.names = colnames(mean.abun.table))
for (trait in trait_names) {
  uw.sample.aver <- colSums(mean.abun.table > 0) / colSums((mean.abun.table > 0) / abun.trait.data[[trait]], na.rm = TRUE)
  abun.unweighted.result[[trait]] <- uw.sample.aver
}

rare.unweighted.result <- data.frame(row.names = colnames(mean.rare.table))
for (trait in trait_names) {
  uw.sample.aver <- colSums(mean.rare.table > 0) / colSums((mean.rare.table > 0) / rare.trait.data[[trait]], na.rm = TRUE)
  rare.unweighted.result[[trait]] <- uw.sample.aver
}

# 去除重复的sample
abun.weighted.result.unique <- abun.weighted.result[!duplicated(abun.weighted.result[, 1]), ]
abun.weighted.result.unique$name <- rownames(abun.weighted.result.unique)
abun.weighted.result.unique$group <- "abundant"
abun.weighted.result.unique$method <- "weighted"

rare.weighted.result.unique <- rare.weighted.result[!duplicated(rare.weighted.result[, 1]), ]
rare.weighted.result.unique$name <- rownames(rare.weighted.result.unique)
rare.weighted.result.unique$group <- "rare"
rare.weighted.result.unique$method <- "weighted"

abun.unweighted.result.unique <- abun.unweighted.result[!duplicated(abun.unweighted.result[, 1]), ]
abun.unweighted.result.unique$name <- rownames(abun.unweighted.result.unique)
abun.unweighted.result.unique$group <- "abundant"
abun.unweighted.result.unique$method <- "unweighted"

rare.unweighted.result.unique <- rare.unweighted.result[!duplicated(rare.unweighted.result[, 1]), ]
rare.unweighted.result.unique$name <- rownames(rare.unweighted.result.unique)
rare.unweighted.result.unique$group <- "rare"
rare.unweighted.result.unique$method <- "unweighted"

result.unique <- rbind(
  abun.weighted.result.unique, rare.weighted.result.unique,
  abun.unweighted.result.unique, rare.unweighted.result.unique
)
library(tidyr)
result.unique <- result.unique %>%
  separate(name, into = c("plot", "layer"), sep = "L", remove = FALSE)
result.unique$layer <- ifelse(result.unique$layer %in% 1:2, "Topsoil", "Subsoil")

# plot ####
plot.data <- result.unique
plot.data <- subset(plot.data, plot.data$method == "weighted")
plot.data <- plot.data[, -16]
library(reshape2)
plot.data <- melt(
  plot.data,
  id.vars = c("name", "plot", "layer", "group"),
  variable.name = "trait",
  value.name = "value"
)
plot.data <- plot.data[plot.data$trait != "doubling_h", ]
plot.data <- plot.data[plot.data$trait != "growth_tmp", ]
library(ggplot2)
library(gghalves)
plot.data$layer <- factor(plot.data$layer, levels = c("Topsoil", "Subsoil"))
plot.data$group <- factor(plot.data$group, levels = c("abundant", "rare"))

plot.data$trait = as.character(plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "cell_diameter", "Cell diameter (μm)",
                         plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "cell_length", "Cell length (μm)",
                         plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "genome_size", "Genome size (Mbp)",
                         plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "gc_content", "GC content (%)",
                         plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "optimum_ph", "Opt. pH",
                         plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "optimum_tmp", "Opt. temp. (°C)",
                         plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "coding_genes", "Coding gene",
                         plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "rRNA16S_genes", "rrn copy number",
                         plot.data$trait)
plot.data$trait = ifelse(plot.data$trait == "tRNA_genes", "tRNA gene",
                         plot.data$trait)
plot.data[plot.data$trait == "Genome size (Mbp)", "value"] <- plot.data[plot.data$trait == "Genome size (Mbp)", "value"] / 1e6

ggplot(
  data = plot.data,
  # aes(x = layer, y = value, fill = layer)
  aes(x = group, y = value, fill = group)
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
  # scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("Topsoil", "Subsoil")) +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("abundant", "rare")) +
  labs(y = "Value", x = NULL) +
  facet_wrap(~trait, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.position = "none",
    strip.text = element_text(size = 13.5),
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

ggsave(paste0("abundant vs rare", "_", "weighted", "_", "bacteria", ".pdf"), width = 7, height = 5.6, units = "in")

# statistics test ####
plot.data <- result.unique
plot.data <- subset(plot.data, plot.data$method == "unweighted")
plot.data <- plot.data[, -16]

div.table <- data.frame(stringsAsFactors = F)

divindex <- plot.data[, 1:11]
treat.used <- plot.data[, -c(1:11)]
treat.used$group <- ifelse(treat.used$group == "rare", 1, 0)
treat.used$block <- treat$block[match(treat.used$name, rownames(treat))]
treat.used$treat <- treat$combined_treat1[match(treat.used$name, rownames(treat))]

divindex <- scale(divindex)
library(car)
library(lme4)
divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ group + (1 | block) + (1 | treat) + (1 | layer), data = div)
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
  fm <- lmer(divtest ~ group + (1 | block) + (1 | layer), data = div)
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
write.csv(div.table, paste0(prefix.m, "_", "unweighted", "_", "abundant vs rare", "_", "P value.csv"))
