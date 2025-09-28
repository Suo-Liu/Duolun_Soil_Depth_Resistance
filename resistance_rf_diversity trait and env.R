setwd("C:\\Users\\True\\OneDrive\\桌面")
# calculation ####
library(randomForest)
library(rfPermute)
library(rfUtilities)
env.data.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\env_used.csv"
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"
treat <- read.csv(treat.file, header = T, row.names = 1)
treat <- subset(treat, plant.type == "TS")

index.name <- "beta"
library(ieggr)

j <- "Sorenson"
## bacteria ####
prefix.m <- "Bacteria"

divindex.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\alpha\\diversity index\\alpha index_", prefix.m, ".csv")
divindex.data <- read.csv(divindex.file, row.names = 1, header = T, sep = ",")

result.data.file <- paste0(
  "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\", index.name,
  "\\resistance\\", index.name, "_resistance_", prefix.m, ".csv"
)
result.data <- read.csv(result.data.file, header = T, row.names = 1)
prefix.m <- paste0(prefix.m, "_", index.name, "_", "TS")
result.data <- subset(result.data, treat == j &
  plant1 == "TS")

### bacteria traits ####
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

### others ####
env.data <- read.csv(env.data.file, header = T, row.names = 1)

env.data <- env.data[, c(1:4, 11)]
sames <- match.name(rn.list = list(
  env.data = env.data,
  treat = treat
), silent = T)
env.data <- sames$env.data
treat <- sames$treat

divindex.data <- divindex.data[rownames(divindex.data) %in% rownames(treat), ]

index.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), , silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  index.table <- rbind(index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), , silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", used.index$treat)
  index.table <- rbind(index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", "W")
  index.table <- rbind(index.table, used.index)
}
index.table <- index.table[order(index.table$treat), ]
index.table <- index.table[order(index.table$block), ]
index.table <- index.table[order(index.table$Layer), ]

env.table <- data.frame(stringsAsFactors = F)
# 同上，保留对照的环境因子数据
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1

  env.table <- rbind(env.table, mean.env)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1
  mean.env$treat <- paste0(i, ".", mean.env$treat)
  env.table <- rbind(env.table, mean.env)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1
  mean.env$treat <- paste0(i, ".", "W")
  env.table <- rbind(env.table, mean.env)
}
env.table <- env.table[order(env.table$treat), ]
env.table <- env.table[order(env.table$block), ]
env.table <- env.table[order(env.table$Layer), ]

trait.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, weighted.result = weighted.result), , silent = T)
  used.index <- sampc$weighted.result
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  trait.table <- rbind(trait.table, used.index)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, weighted.result = weighted.result), , silent = T)
  used.index <- sampc$weighted.result
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", used.index$treat)
  trait.table <- rbind(trait.table, used.index)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, weighted.result = weighted.result), silent = T)
  used.index <- sampc$weighted.result
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", "W")
  trait.table <- rbind(trait.table, used.index)
}
trait.table <- trait.table[order(trait.table$treat), ]
trait.table <- trait.table[order(trait.table$block), ]
trait.table <- trait.table[order(trait.table$Layer), ]

result.data$treat3 <- ifelse(result.data$treat3 %in% c("EP", "RP", "W", "WRP", "WEP"), "C",
  result.data$treat3
)
result.data <- result.data[order(result.data$treat3), ]
result.data <- result.data[order(result.data$block1), ]
result.data <- result.data[order(result.data$layer1), ]

used.index.table <- index.table[, c("richness"), drop = F]
used.env.table <- env.table[, c(
  "root.biomass",
  "NO3", "NH4", "moisture.one.nighbor.aver"
)]
used.trait.table <- trait.table[, c(
  "cell_diameter",
  "cell_length",
  "optimum_tmp",
  # "optimum_ph",
  "genome_size",
  "gc_content",
  "coding_genes",
  "rRNA16S_genes",
  "tRNA_genes"
), drop = F]
dat <- cbind(
  result.data$resistance,
  used.index.table,
  used.env.table, used.trait.table
)

prefix <- prefix.m

dat[is.na(dat)] <- 0
colnames(dat)[1] <- "resistance"
dat <- as.data.frame(scale(dat))

set.seed(1028)
biomass.forest <- randomForest(resistance ~ ., data = dat, importance = TRUE, proximity = TRUE)
# Calculate the overall significance and explanatory power (r-values) of the model
biomass.pvalue <- rf.significance(biomass.forest, dat[, -1], nperm = 1000, ntree = 500)
biomass.pvalue
sink(file = paste0(prefix, "_P value.txt"))
print(biomass.pvalue)
sink()
# Calculate the IncMSE(%) and significance of each variable
biomass.rfP <- rfPermute(resistance ~ .,
  data = dat, importance = TRUE, ntree = 500, nrep = 1000,
  num.cores = 15, proximity = TRUE
)
index.importance.rfp <- data.frame(importance(biomass.rfP, scale = TRUE), check.names = FALSE)
index.importance.rfp$domain <- prefix
save.file(index.importance.rfp, prefix = prefix, filename = "index.importance.rfp")

## Protist ####
prefix.m <- "Protist"

divindex.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\alpha\\diversity index\\alpha index_", prefix.m, ".csv")
divindex.data <- read.csv(divindex.file, row.names = 1, header = T, sep = ",")

result.data.file <- paste0(
  "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\", index.name,
  "\\resistance\\", index.name, "_resistance_", prefix.m, ".csv"
)
result.data <- read.csv(result.data.file, header = T, row.names = 1)

Bac.result.data.file <- paste0(
  "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\", index.name,
  "\\resistance\\", index.name, "_resistance_", "Bacteria", ".csv"
)
Bac.result.data <- read.csv(Bac.result.data.file, header = T, row.names = 1)

prefix.m <- paste0(prefix.m, "_", index.name, "_", "TS")
result.data <- subset(result.data, treat == j &
  plant1 == "TS")
Bac.result.data <- subset(Bac.result.data, treat == j &
  plant1 == "TS")

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"

comm <- t(read.table(com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))
name.row <- rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] <- "YD3L3"
rownames(comm) <- name.row

library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, comm = comm), silent = T)
dim(comm)
comm <- sampc$comm
treat <- sampc$treat
comm <- comm[, colSums(comm) > 0]
dim(comm)

clas.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\18S_tax.txt"
clas <- read.table(clas.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
)
### shell ####
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

shell.result <- shell.weighted.result

### locomotion ####
locomotion.trait <- clas[, c("locomotion"), drop = F]
locomotion.trait <- locomotion.trait[locomotion.trait$locomotion != "", , drop = F]

locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion == c("non_motile"), 0, locomotion.trait$locomotion)
locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion == c("pseudopodia_and_flagella"), 2, locomotion.trait$locomotion)
locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion %in% c(
  "cilia", "flagella", "pseudopodia",
  "pseudopodia_and_flagella_or_pseudopodia_or_flagella"
), 1, locomotion.trait$locomotion)
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

locomotion.result <- locomotion.weighted.result

### bacterivore ####
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

bacterivore.result <- bacterivore.weighted.result
colnames(bacterivore.result)[1] <- "bacterivore"

### body size ####
size.trait <- clas[, c("body.size"), drop = F]
size.trait <- size.trait[!is.na(size.trait$body.size), , drop = F]

sampc <- match.name(cn.list = list(comm = comm), rn.list = list(size.trait = size.trait), silent = T)
size.comm <- sampc$comm
size.trait <- sampc$size.trait

size.otu.table <- t(size.comm)

size.weighted.result <- data.frame(row.names = colnames(size.otu.table))
w.sample.aver <- as.matrix(size.comm) %*% as.matrix(size.trait) / colSums(size.otu.table)
size.weighted.result <- w.sample.aver

size.unweighted.result <- data.frame(row.names = colnames(size.otu.table))
uw.sample.aver <- as.matrix(size.comm > 0) %*% as.matrix(size.trait) / colSums(size.otu.table > 0)
size.unweighted.result <- uw.sample.aver

size.weighted.result <- as.data.frame(size.weighted.result)
size.unweighted.result <- as.data.frame(size.unweighted.result)

size.result <- size.weighted.result
colnames(size.result)[1] <- "size"

trait.result <- cbind(
  shell.result,
  locomotion.result,
  bacterivore.result,
  size.result
)

### others ####
trait.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, trait.result = trait.result), , silent = T)
  used.trait <- sampc$trait.result
  ct.treat <- sampc$ct.treat

  used.trait <- as.data.frame(used.trait)
  used.trait$block <- ct.treat$block
  used.trait$Layer <- ct.treat$Layer
  used.trait$treat <- ct.treat$combined_treat1
  trait.table <- rbind(trait.table, used.trait)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, trait.result = trait.result), , silent = T)
  used.trait <- sampc$trait.result
  ct.treat <- sampc$ct.treat

  used.trait <- as.data.frame(used.trait)
  used.trait$block <- ct.treat$block
  used.trait$Layer <- ct.treat$Layer
  used.trait$treat <- ct.treat$combined_treat1
  used.trait$treat <- paste0(i, ".", used.trait$treat)
  trait.table <- rbind(trait.table, used.trait)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, trait.result = trait.result), silent = T)
  used.trait <- sampc$trait.result
  ct.treat <- sampc$ct.treat

  used.trait <- as.data.frame(used.trait)
  used.trait$block <- ct.treat$block
  used.trait$Layer <- ct.treat$Layer
  used.trait$treat <- ct.treat$combined_treat1
  used.trait$treat <- paste0(i, ".", "W")
  trait.table <- rbind(trait.table, used.trait)
}
trait.table <- trait.table[order(trait.table$treat), ]
trait.table <- trait.table[order(trait.table$block), ]
trait.table <- trait.table[order(trait.table$Layer), ]

env.data <- read.csv(env.data.file, header = T, row.names = 1)

env.data <- env.data[, c(1:4, 11)]
sames <- match.name(rn.list = list(
  env.data = env.data,
  treat = treat
), silent = T)
env.data <- sames$env.data
treat <- sames$treat

divindex.data <- divindex.data[rownames(divindex.data) %in% rownames(treat), ]

index.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), , silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  index.table <- rbind(index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), , silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", used.index$treat)
  index.table <- rbind(index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", "W")
  index.table <- rbind(index.table, used.index)
}
index.table <- index.table[order(index.table$treat), ]
index.table <- index.table[order(index.table$block), ]
index.table <- index.table[order(index.table$Layer), ]

env.table <- data.frame(stringsAsFactors = F)
# 同上，保留对照的环境因子数据
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1

  env.table <- rbind(env.table, mean.env)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1
  mean.env$treat <- paste0(i, ".", mean.env$treat)
  env.table <- rbind(env.table, mean.env)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1
  mean.env$treat <- paste0(i, ".", "W")
  env.table <- rbind(env.table, mean.env)
}
env.table <- env.table[order(env.table$treat), ]
env.table <- env.table[order(env.table$block), ]
env.table <- env.table[order(env.table$Layer), ]

result.data$treat3 <- ifelse(result.data$treat3 %in% c("EP", "RP", "W", "WRP", "WEP"), "C",
  result.data$treat3
)
result.data <- result.data[order(result.data$treat3), ]
result.data <- result.data[order(result.data$block1), ]
result.data <- result.data[order(result.data$layer1), ]

Bac.result.data$treat3 <- ifelse(Bac.result.data$treat3 %in% c("EP", "RP", "W", "WRP", "WEP"), "C",
  Bac.result.data$treat3
)
Bac.result.data <- Bac.result.data[order(Bac.result.data$treat3), ]
Bac.result.data <- Bac.result.data[order(Bac.result.data$block1), ]
Bac.result.data <- Bac.result.data[order(Bac.result.data$layer1), ]

used.index.table <- index.table[, c("richness"), drop = F]
used.env.table <- env.table[, c(
  "root.biomass",
  "NO3", "NH4", "moisture.one.nighbor.aver"
)]
used.trait.table <- trait.table[, c(
  "shell",
  "locomotion",
  "bacterivore",
  "size"
), drop = F]
dat <- cbind(
  result.data$resistance,
  Bac.result.data$resistance,
  used.index.table,
  used.env.table, used.trait.table
)

prefix <- prefix.m

dat[is.na(dat)] <- 0
colnames(dat)[1] <- "resistance"
colnames(dat)[2] <- "Bacteria"
dat <- as.data.frame(scale(dat))

set.seed(1234)
biomass.forest <- randomForest(resistance ~ ., data = dat, importance = TRUE, proximity = TRUE)
# Calculate the overall significance and explanatory power (r-values) of the model
biomass.pvalue <- rf.significance(biomass.forest, dat[, -1], nperm = 1000, ntree = 500)
biomass.pvalue
sink(file = paste0(prefix, "_P value.txt"))
print(biomass.pvalue)
sink()
# Calculate the IncMSE(%) and significance of each variable
biomass.rfP <- rfPermute(resistance ~ .,
  data = dat, importance = TRUE, ntree = 500, nrep = 1000,
  num.cores = 15, proximity = TRUE
)
index.importance.rfp <- data.frame(importance(biomass.rfP, scale = TRUE), check.names = FALSE)
index.importance.rfp$domain <- prefix
save.file(index.importance.rfp, prefix = prefix, filename = "index.importance.rfp")

## Fungi ####
prefix.m <- "Fungi"

divindex.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\alpha\\diversity index\\alpha index_", prefix.m, ".csv")
divindex.data <- read.csv(divindex.file, row.names = 1, header = T, sep = ",")

result.data.file <- paste0(
  "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\", index.name,
  "\\resistance\\", index.name, "_resistance_", prefix.m, ".csv"
)
result.data <- read.csv(result.data.file, header = T, row.names = 1)
prefix.m <- paste0(prefix.m, "_", index.name, "_", "TS")
result.data <- subset(result.data, treat == j &
  plant1 == "TS")

env.data <- read.csv(env.data.file, header = T, row.names = 1)

env.data <- env.data[, c(1:4, 11)]
sames <- match.name(rn.list = list(
  env.data = env.data,
  treat = treat
), silent = T)
env.data <- sames$env.data
treat <- sames$treat

divindex.data <- divindex.data[rownames(divindex.data) %in% rownames(treat), ]

index.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), , silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  index.table <- rbind(index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), , silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", used.index$treat)
  index.table <- rbind(index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(rn.list = list(ct.treat = ct.treat, divindex.data = divindex.data), silent = T)
  used.index <- sampc$divindex.data
  ct.treat <- sampc$ct.treat

  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", "W")
  index.table <- rbind(index.table, used.index)
}
index.table <- index.table[order(index.table$treat), ]
index.table <- index.table[order(index.table$block), ]
index.table <- index.table[order(index.table$Layer), ]

env.table <- data.frame(stringsAsFactors = F)
# 同上，保留对照的环境因子数据
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1

  env.table <- rbind(env.table, mean.env)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1
  mean.env$treat <- paste0(i, ".", mean.env$treat)
  env.table <- rbind(env.table, mean.env)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")

  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]

  sampc <- match.name(
    rn.list = list(
      ct.treat = ct.treat,
      env.data = env.data
    ),
    silent = T
  )
  ct.env.data <- sampc$env.data
  ct.treat <- sampc$ct.treat

  mean.env <- ct.env.data

  mean.env <- as.data.frame(mean.env)
  mean.env$block <- ct.treat$block
  mean.env$Layer <- ct.treat$Layer
  mean.env$treat <- ct.treat$combined_treat1
  mean.env$treat <- paste0(i, ".", "W")
  env.table <- rbind(env.table, mean.env)
}
env.table <- env.table[order(env.table$treat), ]
env.table <- env.table[order(env.table$block), ]
env.table <- env.table[order(env.table$Layer), ]

result.data$treat3 <- ifelse(result.data$treat3 %in% c("EP", "RP", "W", "WRP", "WEP"), "C",
  result.data$treat3
)
result.data <- result.data[order(result.data$treat3), ]
result.data <- result.data[order(result.data$block1), ]
result.data <- result.data[order(result.data$layer1), ]

used.index.table <- index.table[, c("richness"), drop = F]
used.env.table <- env.table[, c(
  "root.biomass",
  "NO3", "NH4", "moisture.one.nighbor.aver"
)]

dat <- cbind(
  result.data$resistance,
  used.index.table,
  used.env.table
)

prefix <- prefix.m

dat[is.na(dat)] <- 0
colnames(dat)[1] <- "resistance"
dat <- as.data.frame(scale(dat))

set.seed(1028)
biomass.forest <- randomForest(resistance ~ ., data = dat, importance = TRUE, proximity = TRUE)
# Calculate the overall significance and explanatory power (r-values) of the model
biomass.pvalue <- rf.significance(biomass.forest, dat[, -1], nperm = 1000, ntree = 500)
biomass.pvalue
sink(file = paste0(prefix, "_P value.txt"))
print(biomass.pvalue)
sink()
# Calculate the IncMSE(%) and significance of each variable
biomass.rfP <- rfPermute(resistance ~ .,
  data = dat, importance = TRUE, ntree = 500, nrep = 1000,
  num.cores = 15, proximity = TRUE
)
index.importance.rfp <- data.frame(importance(biomass.rfP, scale = TRUE), check.names = FALSE)
index.importance.rfp$domain <- prefix
save.file(index.importance.rfp, prefix = prefix, filename = "index.importance.rfp")

# plot ####
library(Rmisc)
i <- 0
plot.data <- data.frame(stringsAsFactors = F)

data.file <- file.choose()
i <- i + 1
data <- read.csv(data.file)

data <- data[order(data$X.IncMSE, decreasing = T), ]

names <- data[, 1][1:6]
used <- data[data$ID %in% names, ]

used <- used[order(used$X.IncMSE, decreasing = T), ]

used <- used[, 1:3]

colnames(used) <- c("ID", "value", "p.value")

all.used <- used
all.used <- rbind(all.used, c("NA", 1, 1))
all.used$ID <- paste0(all.used$ID, ".", i)
all.used$domain <- "Protist"
plot.data <- rbind(plot.data, all.used)

write.csv(plot.data, "plot.data.csv")
# 添加group列
plot.data <- read.csv(file.choose(), row.names = 1)
plot.data$value <- as.numeric(plot.data$value)
plot.data$x <- factor(plot.data$ID, levels = unique(plot.data$ID))
library(ggplot2)
library(ggpattern)

p1 <- ggplot(data = plot.data, aes(x = x, y = value, color = group, fill = group)) +
  geom_bar(colour = "black", stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#CD2626", "#1874CD", "#00A087"), limits = c("Microbial factor", "Edaphic factor", "Botanic factor")) +
  scale_color_manual(values = c("#CD2626", "#1874CD", "#00A087"), limits = c("Microbial factor", "Edaphic factor", "Botanic factor")) +
  xlab("") +
  ylab("IncMSE(%)") +
  theme_bw() +
  coord_flip() +
  theme(
    legend.title = element_text(size = 15, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(
      size = 15,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5,
      color = "black"
    ),
    axis.text = element_text(
      size = 12,
      angle = 0,
      margin = margin(t = 0),
      color = "black"
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(t = 0),
      color = "black"
    )
  )
p1
ggsave("rf_resistance_env_trait_diversity.pdf", width = 7, height = 6, units = "in")
