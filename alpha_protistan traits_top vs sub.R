setwd("C:/Users/True/OneDrive/桌面")

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"
clas.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\18S_tax.txt"
prefix.m <- "Protist"

comm <- t(read.table(com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))
name.row <- rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] <- "YD3L3"
rownames(comm) <- name.row

clas <- read.table(clas.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
)
treat <- read.csv(treat.file, header = T, row.names = 1)

treat <- subset(treat, plant.type == "TS")
prefix.m <- paste0(prefix.m, "_", "TS")

library(ieggr)
sampc <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), silent = T)
comm <- sampc$comm
clas <- sampc$clas

sampc <- match.name(rn.list = list(treat = treat, comm = comm), silent = T)
dim(comm)
comm <- sampc$comm
treat <- sampc$treat
comm <- comm[, colSums(comm) > 0]
dim(comm)

# sporulation, cyst, and omnivore ####
sporu.trait <- clas[, c("sporulation"), drop = F]
sporu.trait <- sporu.trait[sporu.trait$sporulation != "", , drop = F]
sporu.trait$sporulation <- ifelse(sporu.trait$sporulation == "present", 1, 0)

cyst.trait <- clas[, c("cyst"), drop = F]
cyst.trait <- cyst.trait[cyst.trait$cyst != "", , drop = F]
cyst.trait$cyst <- ifelse(cyst.trait$cyst == "present", 1, 0)

omnivore.trait <- clas[, c("main_functional_class", "detailed_functional_class"), drop = F]
omnivore.trait <- omnivore.trait[omnivore.trait$main_functional_class %in% c("predator", "predator (add)"), , drop = F]
omnivore.trait$detailed_functional_class <- ifelse(omnivore.trait$detailed_functional_class %in% c("omnivore", "omnivore (add)"), 1, 0)
omnivore.trait <- omnivore.trait[, -1, drop = F]

sampc <- match.name(cn.list = list(comm = comm), rn.list = list(sporu.trait = sporu.trait), silent = T)
sporu.comm <- sampc$comm
sporu.trait <- sampc$sporu.trait

sampc <- match.name(cn.list = list(comm = comm), rn.list = list(cyst.trait = cyst.trait), silent = T)
cyst.comm <- sampc$comm
cyst.trait <- sampc$cyst.trait

sampc <- match.name(cn.list = list(comm = comm), rn.list = list(omnivore.trait = omnivore.trait), silent = T)
omnivore.comm <- sampc$comm
omnivore.trait <- sampc$omnivore.trait

sporu.otu.table <- t(sporu.comm)
cyst.otu.table <- t(cyst.comm)
omnivore.otu.table <- t(omnivore.comm)

sporu.weighted.result <- data.frame(row.names = colnames(sporu.otu.table))

w.sample.aver <- as.matrix(sporu.comm) %*% as.matrix(sporu.trait) / colSums(sporu.otu.table)
sporu.weighted.result <- w.sample.aver

sporu.unweighted.result <- data.frame(row.names = colnames(sporu.otu.table))

uw.sample.aver <- as.matrix(sporu.comm > 0) %*% as.matrix(sporu.trait) / colSums(sporu.otu.table > 0)
sporu.unweighted.result <- uw.sample.aver

sporu.weighted.result <- as.data.frame(sporu.weighted.result)
sporu.unweighted.result <- as.data.frame(sporu.unweighted.result)

sporu.weighted.result$name <- rownames(sporu.weighted.result)
sporu.weighted.result$method <- "weighted"

sporu.unweighted.result$name <- rownames(sporu.unweighted.result)
sporu.unweighted.result$method <- "unweighted"

sporu.result <- rbind(sporu.weighted.result, sporu.unweighted.result)

cyst.weighted.result <- data.frame(row.names = colnames(cyst.otu.table))

w.sample.aver <- as.matrix(cyst.comm) %*% as.matrix(cyst.trait) / colSums(cyst.otu.table)
cyst.weighted.result <- w.sample.aver

cyst.unweighted.result <- data.frame(row.names = colnames(cyst.otu.table))
uw.sample.aver <- as.matrix(cyst.comm > 0) %*% as.matrix(cyst.trait) / colSums(cyst.otu.table > 0)
cyst.unweighted.result <- uw.sample.aver

cyst.weighted.result <- as.data.frame(cyst.weighted.result)
cyst.unweighted.result <- as.data.frame(cyst.unweighted.result)
cyst.weighted.result$name <- rownames(cyst.weighted.result)
cyst.weighted.result$method <- "weighted"

cyst.unweighted.result$name <- rownames(cyst.unweighted.result)
cyst.unweighted.result$method <- "unweighted"

cyst.result <- rbind(cyst.weighted.result, cyst.unweighted.result)

omnivore.weighted.result <- data.frame(row.names = colnames(omnivore.otu.table))
w.sample.aver <- as.matrix(omnivore.comm) %*% as.matrix(omnivore.trait) / colSums(omnivore.otu.table)
omnivore.weighted.result <- w.sample.aver

omnivore.unweighted.result <- data.frame(row.names = colnames(omnivore.otu.table))
uw.sample.aver <- as.matrix(omnivore.comm > 0) %*% as.matrix(omnivore.trait) / colSums(omnivore.otu.table > 0)
omnivore.unweighted.result <- uw.sample.aver

omnivore.weighted.result <- as.data.frame(omnivore.weighted.result)
omnivore.unweighted.result <- as.data.frame(omnivore.unweighted.result)

omnivore.weighted.result$name <- rownames(omnivore.weighted.result)
omnivore.weighted.result$method <- "weighted"

omnivore.unweighted.result$name <- rownames(omnivore.unweighted.result)
omnivore.unweighted.result$method <- "unweighted"

omnivore.result <- rbind(omnivore.weighted.result, omnivore.unweighted.result)
colnames(omnivore.result)[1] <- "omnivore"

table(sporu.result$name == cyst.result$name)
table(sporu.result$name == omnivore.result$name)

# shell ####
shell.trait <- clas[, c("shell"), drop = F]
shell.trait <- shell.trait[shell.trait$shell != "", , drop = F]
shell.trait <- shell.trait[shell.trait$shell != "naked_or_silica_or_organic", , drop = F]
shell.trait$shell <- ifelse(shell.trait$shell == c("naked"), 0, 1)

sampc <- match.name(cn.list = list(comm = comm), rn.list = list(shell.trait = shell.trait), silent = T)
shell.comm <- sampc$comm
shell.trait <- sampc$shell.trait

shell.otu.table <- t(shell.comm)

shell.weighted.result <- data.frame(row.names = colnames(shell.otu.table))

w.sample.aver <- as.matrix(shell.comm) %*% as.matrix(shell.trait) / colSums(shell.otu.table)
shell.weighted.result <- w.sample.aver

shell.unweighted.result <- data.frame(row.names = colnames(shell.otu.table))

uw.sample.aver <- as.matrix(shell.comm > 0) %*% as.matrix(shell.trait) / colSums(shell.otu.table > 0)
shell.unweighted.result <- uw.sample.aver

shell.weighted.result <- as.data.frame(shell.weighted.result)
shell.unweighted.result <- as.data.frame(shell.unweighted.result)

shell.weighted.result$name <- rownames(shell.weighted.result)
shell.weighted.result$method <- "weighted"

shell.unweighted.result$name <- rownames(shell.unweighted.result)
shell.unweighted.result$method <- "unweighted"

shell.result <- rbind(shell.weighted.result, shell.unweighted.result)

# locomotion ####
locomotion.trait <- clas[, c("locomotion"), drop = F]
locomotion.trait <- locomotion.trait[locomotion.trait$locomotion != "", , drop = F]

locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion == c("non_motile"), 0, locomotion.trait$locomotion)
locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion == c("pseudopodia_and_flagella"), 2, locomotion.trait$locomotion)
locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion %in% c(
  "cilia", "flagella", "pseudopodia",
  "pseudopodia_and_flagella_or_pseudopodia_or_flagella"
), 1, locomotion.trait$locomotion)
# locomotion.trait$locomotion = ifelse(locomotion.trait$locomotion == c("pseudopodia_and_flagella"), 1, 0)

# locomotion.trait$locomotion = ifelse(locomotion.trait$locomotion == c("non_motile"), 0, 1)

locomotion.trait$locomotion <- as.numeric(locomotion.trait$locomotion)
sampc <- match.name(cn.list = list(comm = comm), rn.list = list(locomotion.trait = locomotion.trait), silent = T)
locomotion.comm <- sampc$comm
locomotion.trait <- sampc$locomotion.trait

locomotion.otu.table <- t(locomotion.comm)

locomotion.weighted.result <- data.frame(row.names = colnames(locomotion.otu.table))

w.sample.aver <- as.matrix(locomotion.comm) %*% as.matrix(locomotion.trait) / colSums(locomotion.otu.table)
locomotion.weighted.result <- w.sample.aver

locomotion.unweighted.result <- data.frame(row.names = colnames(locomotion.otu.table))

uw.sample.aver <- as.matrix(locomotion.comm > 0) %*% as.matrix(locomotion.trait) / colSums(locomotion.otu.table > 0)
locomotion.unweighted.result <- uw.sample.aver

locomotion.weighted.result <- as.data.frame(locomotion.weighted.result)
locomotion.unweighted.result <- as.data.frame(locomotion.unweighted.result)

locomotion.weighted.result$name <- rownames(locomotion.weighted.result)
locomotion.weighted.result$method <- "weighted"

locomotion.unweighted.result$name <- rownames(locomotion.unweighted.result)
locomotion.unweighted.result$method <- "unweighted"

locomotion.result <- rbind(locomotion.weighted.result, locomotion.unweighted.result)

# bacterivore ####
bacterivore.trait <- clas[, c("main_functional_class", "associated_organism"), drop = F]
bacterivore.trait <- bacterivore.trait[bacterivore.trait$main_functional_class %in% c("predator", "predator (add)"), , drop = F]
bacterivore.trait$associated_organism <- ifelse(grepl("bacteria", bacterivore.trait$associated_organism), 1, 0)
bacterivore.trait <- bacterivore.trait[, -1, drop = F]

sampc <- match.name(cn.list = list(comm = comm), rn.list = list(bacterivore.trait = bacterivore.trait), silent = T)
bacterivore.comm <- sampc$comm
bacterivore.trait <- sampc$bacterivore.trait

bacterivore.otu.table <- t(bacterivore.comm)

bacterivore.weighted.result <- data.frame(row.names = colnames(bacterivore.otu.table))
w.sample.aver <- as.matrix(bacterivore.comm) %*% as.matrix(bacterivore.trait) / colSums(bacterivore.otu.table)
bacterivore.weighted.result <- w.sample.aver

bacterivore.unweighted.result <- data.frame(row.names = colnames(bacterivore.otu.table))
uw.sample.aver <- as.matrix(bacterivore.comm > 0) %*% as.matrix(bacterivore.trait) / colSums(bacterivore.otu.table > 0)
bacterivore.unweighted.result <- uw.sample.aver

bacterivore.weighted.result <- as.data.frame(bacterivore.weighted.result)
bacterivore.unweighted.result <- as.data.frame(bacterivore.unweighted.result)

bacterivore.weighted.result$name <- rownames(bacterivore.weighted.result)
bacterivore.weighted.result$method <- "weighted"

bacterivore.unweighted.result$name <- rownames(bacterivore.unweighted.result)
bacterivore.unweighted.result$method <- "unweighted"

bacterivore.result <- rbind(bacterivore.weighted.result, bacterivore.unweighted.result)
colnames(bacterivore.result)[1] <- "bacterivore"

sum(sporu.result$name == cyst.result$name)
sum(sporu.result$name == omnivore.result$name)
sum(sporu.result$name == bacterivore.result$name)
sum(sporu.result$name == locomotion.result$name)
sum(sporu.result$name == shell.result$name)

trait.result <- cbind(
  sporu.result[, c(2, 3, 1)], cyst.result[, "cyst"], omnivore.result[, "omnivore"],
  bacterivore.result[, "bacterivore"], locomotion.result[, "locomotion"], shell.result[, "shell"]
)
colnames(trait.result)[3:8] <- c("Sporulation", "Cyst", "Omnivore", "Bacterivore", "Locomotion", "Shell")

trait.result$Layer <- treat$Layer[match(trait.result$name, rownames(treat))]
trait.result$plot <- treat$plot[match(trait.result$name, rownames(treat))]
trait.result$block <- treat$block[match(trait.result$name, rownames(treat))]
trait.result$treat <- treat$combined_treat1[match(trait.result$name, rownames(treat))]

# plot ####
plot.data <- trait.result
plot.data <- subset(plot.data, plot.data$method == "unweighted")
plot.data <- plot.data[, -2]
library(reshape2)
plot.data <- melt(
  plot.data,
  id.vars = c("name", "plot", "block", "Layer", "treat"),
  variable.name = "trait",
  value.name = "value"
)
library(ggplot2)
library(gghalves)
library(dplyr)
plot.data <- subset(plot.data, trait %in% c("Bacterivore", "Locomotion", "Shell"))
plot.data$Layer <- ifelse(plot.data$Layer %in% c("L1", "L2"), "Topsoil", "Subsoil")
plot.data$Layer <- factor(plot.data$Layer, levels = c("Topsoil", "Subsoil"))

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

ggsave(paste0("top vs sub", "_", "unweighted", "_", "protistan trait", ".pdf"), width = 6.99, height = 2.22, units = "in")

# statistics test ####
plot.data <- trait.result
plot.data <- subset(plot.data, plot.data$method == "unweighted")
plot.data <- plot.data[, -2]

div.table <- data.frame(stringsAsFactors = F)

divindex <- plot.data[, 2:7]
treat.used <- plot.data[, -c(2:7)]
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
write.csv(div.table, paste0("unweighted", "_", "protistan trait", "_", "P value.csv"))
