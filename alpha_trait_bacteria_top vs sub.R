setwd("C:/Users/True/OneDrive/桌面")
# trait calculation ####
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"

comm <- t(read.table(com.file,
                     header = TRUE, sep = "\t", row.names = 1,
                     as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE
))
name.row <- rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] <- "YD3L3"
rownames(comm) <- name.row

treat <- read.csv(treat.file, header = T, row.names = 1)

treat <- subset(treat, plant.type == "TS")

library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, comm = comm), silent = T)
dim(comm)
comm <- sampc$comm
treat <- sampc$treat
comm <- comm[, colSums(comm) > 0]
dim(comm)

trait.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\predicted_traits.csv"

trait.data <- read.csv(trait.file, header = T, row.names = 1)
trait.data$cell_diameter <- 10^trait.data$cell_diameter
trait.data$cell_length <- 10^trait.data$cell_length
trait.data$doubling_h <- 10^trait.data$doubling_h
trait.data <- trait.data[, 1:11]

otu.table <- t(comm)

spc <- match.name(rn.list = list(otu.table = otu.table, trait.data = trait.data))
otu.table <- spc$otu.table
trait.data <- spc$trait.data

data.list <- list()
for (i in 1:ncol(trait.data)) {
  for (j in 1:ncol(otu.table)) {
    otu.table[, j] <- otu.table[, j] * trait.data[, i]
  }
  data.list[[i]] <- otu.table
  names(data.list)[i] <- colnames(trait.data)[i]
}

trait_names <- names(trait.data)
used.adjusted_mats <- list()
for (trait in trait_names) {
  adjusted_mat <- otu.table / trait.data[[trait]]
  adjusted_mat[is.infinite(adjusted_mat)] <- NA
  
  used.adjusted_mats[[trait]] <- t(adjusted_mat)
}
weighted.result <- data.frame(row.names = colnames(otu.table))
for (trait in trait_names) {
  w.sample.aver <- colSums(otu.table) / colSums(otu.table / trait.data[[trait]], na.rm = TRUE)
  weighted.result[[trait]] <- w.sample.aver
}

unweighted.result <- data.frame(row.names = colnames(otu.table))
for (trait in trait_names) {
  uw.sample.aver <- colSums(otu.table > 0) / colSums((otu.table > 0) / trait.data[[trait]], na.rm = TRUE)
  unweighted.result[[trait]] <- uw.sample.aver
}

weighted.result$name <- rownames(weighted.result)
weighted.result$method <- "weighted"

unweighted.result$name <- rownames(unweighted.result)
unweighted.result$method <- "unweighted"

result <- rbind(weighted.result, unweighted.result)

result$Layer <- treat$Layer[match(result$name, rownames(treat))]
result$plot <- treat$plot[match(result$name, rownames(treat))]
result$block <- treat$block[match(result$name, rownames(treat))]
result$treat <- treat$combined_treat1[match(result$name, rownames(treat))]

# plot ####
plot.data <- result
plot.data <- subset(plot.data, plot.data$method == "weighted")
plot.data <- plot.data[, -13]
library(reshape2)
plot.data <- melt(
  plot.data,
  id.vars = c("name", "plot", "block", "Layer", "treat"),
  variable.name = "trait",
  value.name = "value"
)
plot.data <- plot.data[plot.data$trait != "doubling_h", ]
plot.data <- plot.data[plot.data$trait != "growth_tmp", ]
library(ggplot2)
library(gghalves)

plot.data$Layer <- ifelse(plot.data$Layer %in% c("L1", "L2"), "Topsoil", "Subsoil")
plot.data$Layer <- factor(plot.data$Layer, levels = c("Topsoil", "Subsoil"))

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
  aes(x = Layer, y = value, fill = Layer)
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

ggsave(paste0("top vs sub","_","weighted","_","bacteria",".pdf"), width = 7, height = 5.6, units = "in")

# statistics test ####
plot.data <- result
plot.data <- subset(plot.data, plot.data$method == "weighted")
plot.data <- plot.data[, -13]

div.table <- data.frame(stringsAsFactors = F)

divindex <- plot.data[, 1:11]
treat.used <- plot.data[, -c(1:11)]
treat.used$Layer <- ifelse(treat.used$Layer %in% c("L1", "L2"), "Topsoil", "Subsoil")
treat.used$Layer <- ifelse(treat.used$Layer == "Topsoil", 0, 1)
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
write.csv(div.table, paste0(prefix.m, "_", "weighted", "_", "all treatments", "_", "P value.csv"))
