setwd("C:/Users/True/OneDrive/桌面")
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"
# abundant and rare species ####
com.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"
prefix.m = "Protist"

comm=t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
name.row = rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] = "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] = "YD3L3"
rownames(comm) = name.row

pre.treat <- read.csv(treat.file, header = T, row.names = 1)

pre.treat <- subset(pre.treat, plant.type == "TS")
prefix.m = paste0(prefix.m, "_","TS")

treat = pre.treat
library(ieggr)
sampc=match.name(rn.list = list(treat = treat,comm = comm), silent = T)
dim(comm)
comm = sampc$comm
treat = sampc$treat
comm = comm[,colSums(comm)>0]
dim(comm)

# 统计物种出现频率
library(dplyr)
fre.comm = comm
fre.comm = as.data.frame(fre.comm)
fre.comm$Layer = treat$Layer
fre.comm$Layer = ifelse(fre.comm$Layer %in% c("L1","L2"),
                        "Topsoils","Subsoils")

fre.comm %<>% group_by(Layer) %>%
  summarise_all(sum)
fre.comm = as.data.frame(fre.comm)
rownames(fre.comm) = fre.comm[,1]
fre.comm = fre.comm[,-1]
# 识别交集中的物种和独有物种
inter.otu = names(fre.comm)[which(colSums(fre.comm>0) == nrow(fre.comm))]
# 识别丰富物种和稀有物种，定义源于Dai et al., NC, 2022，但是丰富物种的定义少了出现频率的限制，因为过多限制会导致最后识别的丰富物种太少
abun.otu = which(colSums(comm/rowSums(comm) >= 0.001) > 0.5 & (colSums(comm>0)>(0.5*nrow(comm))))
# abun.otu = which(colSums(comm/rowSums(comm) >= 0.001) > 0.5)
abun.otu = abun.otu[which(names(abun.otu)%in%inter.otu)]

rare.otu = which(colSums(comm/rowSums(comm) < 0.001) == nrow(comm))

abun.comm = comm[,colnames(comm) %in% names(abun.otu)]
rare.comm = comm[,colnames(comm) %in% names(rare.otu)]

# 因为要用物种多样性和丰度数据预测抗性，而抗性是通过处理与对照间的群落差异计算得到的，
# 因此要match对照的物种丰度（应该不需要对照和处理的平均丰度？）
mean.abun.table = data.frame(stringsAsFactors = F)
mean.rare.table = data.frame(stringsAsFactors = F)
for (i in c("RP","EP","W","WEP","WRP")){
  ct.treat = subset(treat, combined_treat1 == "C")
  ct.treat = ct.treat[order(ct.treat$block),]
  ct.treat = ct.treat[order(ct.treat$Layer),]
  
  sampc=match.name(rn.list = list(ct.treat = ct.treat, abun.comm = abun.comm,
                                  rare.comm = rare.comm),, silent = T)
  ct.abun.comm = sampc$abun.comm
  ct.rare.comm = sampc$rare.comm
  ct.treat = sampc$ct.treat
  
  mean.abun.comm = ct.abun.comm
  mean.rare.comm = ct.rare.comm
  
  mean.abun.comm = as.data.frame(mean.abun.comm)
  mean.abun.comm$block = ct.treat$block
  mean.abun.comm$Layer = ct.treat$Layer
  mean.abun.comm$treat = ct.treat$combined_treat1
  
  mean.abun.table = rbind(mean.abun.table, mean.abun.comm)
  
  mean.rare.comm = as.data.frame(mean.rare.comm)
  mean.rare.comm$block = ct.treat$block
  mean.rare.comm$Layer = ct.treat$Layer
  mean.rare.comm$treat = ct.treat$combined_treat1
  mean.rare.table = rbind(mean.rare.table, mean.rare.comm)
}
for (i in c("WEP","WRP")){
  if (i == "WEP"){
    ct.treat = subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat = subset(treat, combined_treat1 == "RP")
  }
  
  ct.treat = ct.treat[order(ct.treat$block),]
  ct.treat = ct.treat[order(ct.treat$Layer),]
  
  sampc=match.name(rn.list = list(ct.treat = ct.treat, abun.comm = abun.comm,
                                  rare.comm = rare.comm),, silent = T)
  ct.abun.comm = sampc$abun.comm
  ct.rare.comm = sampc$rare.comm
  ct.treat = sampc$ct.treat
  
  mean.abun.comm = ct.abun.comm
  mean.rare.comm = ct.rare.comm
  
  mean.abun.comm = as.data.frame(mean.abun.comm)
  mean.abun.comm$block = ct.treat$block
  mean.abun.comm$Layer = ct.treat$Layer
  mean.abun.comm$treat = ct.treat$combined_treat1
  mean.abun.comm$treat = paste0(i,".",mean.abun.comm$treat)
  mean.abun.table = rbind(mean.abun.table, mean.abun.comm)
  
  mean.rare.comm = as.data.frame(mean.rare.comm)
  mean.rare.comm$block = ct.treat$block
  mean.rare.comm$Layer = ct.treat$Layer
  mean.rare.comm$treat = ct.treat$combined_treat1
  mean.rare.comm$treat = paste0(i,".",mean.rare.comm$treat)
  mean.rare.table = rbind(mean.rare.table, mean.rare.comm)
}
for (i in c("WEP","WRP")){
  ct.treat = subset(treat, combined_treat1 == "W")
  ct.treat = ct.treat[order(ct.treat$block),]
  ct.treat = ct.treat[order(ct.treat$Layer),]
  
  sampc=match.name(rn.list = list(ct.treat = ct.treat, abun.comm = abun.comm,
                                  rare.comm = rare.comm), silent = T)
  ct.abun.comm = sampc$abun.comm
  ct.rare.comm = sampc$rare.comm
  ct.treat = sampc$ct.treat
  
  mean.abun.comm = ct.abun.comm
  mean.rare.comm = ct.rare.comm
  
  mean.abun.comm = as.data.frame(mean.abun.comm)
  mean.abun.comm$block = ct.treat$block
  mean.abun.comm$Layer = ct.treat$Layer
  mean.abun.comm$treat = ct.treat$combined_treat1
  mean.abun.comm$treat = paste0(i, ".","W")
  mean.abun.table = rbind(mean.abun.table, mean.abun.comm)
  
  mean.rare.comm = as.data.frame(mean.rare.comm)
  mean.rare.comm$block = ct.treat$block
  mean.rare.comm$Layer = ct.treat$Layer
  mean.rare.comm$treat = ct.treat$combined_treat1
  mean.rare.comm$treat = paste0(i, ".","W")
  mean.rare.table = rbind(mean.rare.table, mean.rare.comm)
}

mean.abun.table = mean.abun.table[order(mean.abun.table$treat),]
mean.abun.table = mean.abun.table[order(mean.abun.table$block),]
mean.abun.table = mean.abun.table[order(mean.abun.table$Layer),]

mean.rare.table = mean.rare.table[order(mean.rare.table$treat),]
mean.rare.table = mean.rare.table[order(mean.rare.table$block),]
mean.rare.table = mean.rare.table[order(mean.rare.table$Layer),]

mean.abun.table = mean.abun.table[,-c((ncol(mean.abun.table)-2):ncol(mean.abun.table))]
mean.rare.table = mean.rare.table[,-c((ncol(mean.rare.table)-2):ncol(mean.rare.table))]

# protistan trait abundance ####
clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\18S_tax.txt"
clas = read.table(clas.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE)

mean.abun.table = t(mean.abun.table)
mean.rare.table <- t(mean.rare.table)
# shell ####
shell.trait = clas[,c("shell"), drop = F]
shell.trait = shell.trait[shell.trait$shell!="", , drop = F]
shell.trait = shell.trait[shell.trait$shell!="naked_or_silica_or_organic", , drop = F]
shell.trait$shell = ifelse(shell.trait$shell == c("naked"), 0, 1)

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
treat <- subset(treat, plant.type == "TS")

spc <- match.name(rn.list = list(mean.abun.table = mean.abun.table, shell.trait = shell.trait), silent = T)
shell.abun.table <- spc$mean.abun.table
abun.shell.trait <- spc$shell.trait

spc <- match.name(rn.list = list(mean.rare.table = mean.rare.table, shell.trait = shell.trait), silent = T)
shell.rare.table <- spc$mean.rare.table
rare.shell.trait <- spc$shell.trait

new.w.shell.abun.table = shell.abun.table
new.w.shell.abun.table = as.data.frame(new.w.shell.abun.table)
for (j in 1:ncol(shell.abun.table)){
  new.w.shell.abun.table[, j] <- shell.abun.table[, j] * abun.shell.trait
}
new.w.shell.rare.table = shell.rare.table
new.w.shell.rare.table = as.data.frame(new.w.shell.rare.table)
for (j in 1:ncol(shell.rare.table)){
  new.w.shell.rare.table[, j] <- shell.rare.table[, j] * rare.shell.trait 
}

new.uw.shell.abun.table = shell.abun.table
new.uw.shell.abun.table = as.data.frame(new.uw.shell.abun.table)
for (j in 1:ncol(shell.abun.table)){
  new.uw.shell.abun.table[, j] <- shell.abun.table[, j] * abun.shell.trait
}
new.uw.shell.rare.table = shell.rare.table
new.uw.shell.rare.table = as.data.frame(new.uw.shell.rare.table)
for (j in 1:ncol(shell.rare.table)){
  new.uw.shell.rare.table[, j] <- shell.rare.table[, j] * rare.shell.trait 
}

abun.weighted <- colSums(new.w.shell.abun.table) / colSums(shell.abun.table)
rare.weighted <- colSums(new.w.shell.rare.table) / colSums(shell.rare.table)
abun.unweighted <- colSums(new.uw.shell.abun.table) / colSums(shell.abun.table > 0) 
rare.unweighted <- colSums(new.uw.shell.rare.table) / colSums(shell.rare.table > 0)

# 去除重复的sample
abun.weighted.unique <- abun.weighted[names(abun.weighted) %in% rownames(treat)]
abun.weighted.unique = as.data.frame(abun.weighted.unique)
colnames(abun.weighted.unique) = "value"
abun.weighted.unique$name = rownames(abun.weighted.unique)
abun.weighted.unique$group = "abundant"
abun.weighted.unique$method = "weighted"

rare.weighted.unique <- rare.weighted[names(rare.weighted) %in% rownames(treat)]
rare.weighted.unique = as.data.frame(rare.weighted.unique)
colnames(rare.weighted.unique) = "value"
rare.weighted.unique$name = rownames(rare.weighted.unique)
rare.weighted.unique$group = "rare"
rare.weighted.unique$method = "weighted"

abun.unweighted.unique <- abun.unweighted[names(abun.unweighted) %in% rownames(treat)]
abun.unweighted.unique = as.data.frame(abun.unweighted.unique)
colnames(abun.unweighted.unique) = "value"
abun.unweighted.unique$name = rownames(abun.unweighted.unique)
abun.unweighted.unique$group = "abundant"
abun.unweighted.unique$method = "unweighted"

rare.unweighted.unique <- rare.unweighted[names(rare.unweighted) %in% rownames(treat)]
rare.unweighted.unique = as.data.frame(rare.unweighted.unique)
colnames(rare.unweighted.unique) = "value"
rare.unweighted.unique$name = rownames(rare.unweighted.unique)
rare.unweighted.unique$group = "rare"
rare.unweighted.unique$method = "unweighted"

result.unique = rbind(abun.weighted.unique, rare.weighted.unique,
                      abun.unweighted.unique, rare.unweighted.unique)
shell.unique = result.unique

# locomotion ####
locomotion.trait = clas[,c("locomotion"), drop = F]
locomotion.trait = locomotion.trait[locomotion.trait$locomotion!="", , drop = F]

locomotion.trait$locomotion = ifelse(locomotion.trait$locomotion == c("non_motile"), 0, locomotion.trait$locomotion)
locomotion.trait$locomotion = ifelse(locomotion.trait$locomotion == c("pseudopodia_and_flagella"), 2, locomotion.trait$locomotion)
locomotion.trait$locomotion = ifelse(locomotion.trait$locomotion %in% c("cilia","flagella","pseudopodia",
                                                                        "pseudopodia_and_flagella_or_pseudopodia_or_flagella"), 1, locomotion.trait$locomotion)
locomotion.trait$locomotion = as.numeric(locomotion.trait$locomotion)

spc <- match.name(rn.list = list(mean.abun.table = mean.abun.table, locomotion.trait = locomotion.trait), silent = T)
locomotion.abun.table <- spc$mean.abun.table
abun.locomotion.trait <- spc$locomotion.trait

spc <- match.name(rn.list = list(mean.rare.table = mean.rare.table, locomotion.trait = locomotion.trait), silent = T)
locomotion.rare.table <- spc$mean.rare.table
rare.locomotion.trait <- spc$locomotion.trait

new.w.locomotion.abun.table = locomotion.abun.table
new.w.locomotion.abun.table = as.data.frame(new.w.locomotion.abun.table)
for (j in 1:ncol(locomotion.abun.table)){
  new.w.locomotion.abun.table[, j] <- locomotion.abun.table[, j] * abun.locomotion.trait
}

new.w.locomotion.rare.table = locomotion.rare.table
new.w.locomotion.rare.table = as.data.frame(new.w.locomotion.rare.table)
for (j in 1:ncol(locomotion.rare.table)){
  new.w.locomotion.rare.table[, j] <- locomotion.rare.table[, j] * rare.locomotion.trait
}

new.uw.locomotion.abun.table = locomotion.abun.table
new.uw.locomotion.abun.table = as.data.frame(new.uw.locomotion.abun.table)
for (j in 1:ncol(locomotion.abun.table)){
  new.uw.locomotion.abun.table[, j] <- locomotion.abun.table[, j] * abun.locomotion.trait
}

new.uw.locomotion.rare.table = locomotion.rare.table
new.uw.locomotion.rare.table = as.data.frame(new.uw.locomotion.rare.table)
for (j in 1:ncol(locomotion.rare.table)){
  new.uw.locomotion.rare.table[, j] <- locomotion.rare.table[, j] * rare.locomotion.trait
}

abun.weighted <- colSums(new.w.locomotion.abun.table) / colSums(locomotion.abun.table)
rare.weighted <- colSums(new.w.locomotion.rare.table) / colSums(locomotion.rare.table)
abun.unweighted <- colSums(new.uw.locomotion.abun.table > 0) / colSums(locomotion.abun.table > 0) 
rare.unweighted <- colSums(new.uw.locomotion.rare.table > 0) / colSums(locomotion.rare.table > 0)

abun.weighted.unique <- abun.weighted[names(abun.weighted) %in% rownames(treat)]
abun.weighted.unique = as.data.frame(abun.weighted.unique)
colnames(abun.weighted.unique) = "value"
abun.weighted.unique$name = rownames(abun.weighted.unique)
abun.weighted.unique$group = "abundant"
abun.weighted.unique$method = "weighted"

rare.weighted.unique <- rare.weighted[names(rare.weighted) %in% rownames(treat)]
rare.weighted.unique = as.data.frame(rare.weighted.unique)
colnames(rare.weighted.unique) = "value"
rare.weighted.unique$name = rownames(rare.weighted.unique)
rare.weighted.unique$group = "rare"
rare.weighted.unique$method = "weighted"

abun.unweighted.unique <- abun.unweighted[names(abun.unweighted) %in% rownames(treat)]
abun.unweighted.unique = as.data.frame(abun.unweighted.unique)
colnames(abun.unweighted.unique) = "value"
abun.unweighted.unique$name = rownames(abun.unweighted.unique)
abun.unweighted.unique$group = "abundant"
abun.unweighted.unique$method = "unweighted"

rare.unweighted.unique <- rare.unweighted[names(rare.unweighted) %in% rownames(treat)]
rare.unweighted.unique = as.data.frame(rare.unweighted.unique)
colnames(rare.unweighted.unique) = "value"
rare.unweighted.unique$name = rownames(rare.unweighted.unique)
rare.unweighted.unique$group = "rare"
rare.unweighted.unique$method = "unweighted"

result.unique = rbind(abun.weighted.unique, rare.weighted.unique,
                      abun.unweighted.unique, rare.unweighted.unique)

locomotion.unique = result.unique

# bacterivore ####
bacterivore.trait = clas[,c("main_functional_class", "associated_organism"), drop = F]
bacterivore.trait = bacterivore.trait[bacterivore.trait$main_functional_class %in% c("predator", "predator (add)"), , drop = F]
bacterivore.trait$associated_organism = ifelse(grepl("bacteria", bacterivore.trait$associated_organism), 1, 0)
bacterivore.trait = bacterivore.trait[,-1,drop = F]

spc <- match.name(rn.list = list(mean.abun.table = mean.abun.table, bacterivore.trait = bacterivore.trait), silent = T)
bacterivore.abun.table <- spc$mean.abun.table
abun.bacterivore.trait <- spc$bacterivore.trait

spc <- match.name(rn.list = list(mean.rare.table = mean.rare.table, bacterivore.trait = bacterivore.trait), silent = T)
bacterivore.rare.table <- spc$mean.rare.table
rare.bacterivore.trait <- spc$bacterivore.trait

new.w.bacterivore.abun.table = bacterivore.abun.table
new.w.bacterivore.abun.table = as.data.frame(new.w.bacterivore.abun.table)
for (j in 1:ncol(bacterivore.abun.table)){
  new.w.bacterivore.abun.table[, j] <- bacterivore.abun.table[, j] * abun.bacterivore.trait
}
new.w.bacterivore.rare.table = bacterivore.rare.table
new.w.bacterivore.rare.table = as.data.frame(new.w.bacterivore.rare.table)
for (j in 1:ncol(bacterivore.rare.table)){
  new.w.bacterivore.rare.table[, j] <- bacterivore.rare.table[, j] * rare.bacterivore.trait
}

new.uw.bacterivore.abun.table = bacterivore.abun.table
new.uw.bacterivore.abun.table = as.data.frame(new.uw.bacterivore.abun.table)
for (j in 1:ncol(bacterivore.abun.table)){
  new.uw.bacterivore.abun.table[, j] <- bacterivore.abun.table[, j] * abun.bacterivore.trait
}
new.uw.bacterivore.rare.table = bacterivore.rare.table
new.uw.bacterivore.rare.table = as.data.frame(new.uw.bacterivore.rare.table)
for (j in 1:ncol(bacterivore.rare.table)){
  new.uw.bacterivore.rare.table[, j] <- bacterivore.rare.table[, j] * rare.bacterivore.trait
}

abun.weighted <- colSums(new.w.bacterivore.abun.table) / colSums(bacterivore.abun.table)
rare.weighted <- colSums(new.w.bacterivore.rare.table) / colSums(bacterivore.rare.table)
abun.unweighted <- colSums(new.uw.bacterivore.abun.table) / colSums(bacterivore.abun.table > 0)
rare.unweighted <- colSums(new.uw.bacterivore.rare.table) / colSums(bacterivore.rare.table > 0)

# 去除重复的sample
abun.weighted.unique <- abun.weighted[names(abun.weighted) %in% rownames(treat)]
abun.weighted.unique = as.data.frame(abun.weighted.unique)
colnames(abun.weighted.unique) = "value"
abun.weighted.unique$name = rownames(abun.weighted.unique)
abun.weighted.unique$group = "abundant"
abun.weighted.unique$method = "weighted"

rare.weighted.unique <- rare.weighted[names(rare.weighted) %in% rownames(treat)]
rare.weighted.unique = as.data.frame(rare.weighted.unique)
colnames(rare.weighted.unique) = "value"
rare.weighted.unique$name = rownames(rare.weighted.unique)
rare.weighted.unique$group = "rare"
rare.weighted.unique$method = "weighted"

abun.unweighted.unique <- abun.unweighted[names(abun.unweighted) %in% rownames(treat)]
abun.unweighted.unique = as.data.frame(abun.unweighted.unique)
colnames(abun.unweighted.unique) = "value"
abun.unweighted.unique$name = rownames(abun.unweighted.unique)
abun.unweighted.unique$group = "abundant"
abun.unweighted.unique$method = "unweighted"

rare.unweighted.unique <- rare.unweighted[names(rare.unweighted) %in% rownames(treat)]
rare.unweighted.unique = as.data.frame(rare.unweighted.unique)
colnames(rare.unweighted.unique) = "value"
rare.unweighted.unique$name = rownames(rare.unweighted.unique)
rare.unweighted.unique$group = "rare"
rare.unweighted.unique$method = "unweighted"

result.unique = rbind(abun.weighted.unique, rare.weighted.unique,
                      abun.unweighted.unique, rare.unweighted.unique)

bacterivore.unique = result.unique

# body size ####
body_size.trait = clas[, c("body.size"), drop = F]
body_size.trait = body_size.trait[!is.na(body_size.trait$body.size), , drop = F]

spc <- match.name(rn.list = list(mean.abun.table = mean.abun.table, body_size.trait = body_size.trait), silent = T)
body_size.abun.table <- spc$mean.abun.table
abun.body_size.trait <- spc$body_size.trait

spc <- match.name(rn.list = list(mean.rare.table = mean.rare.table, body_size.trait = body_size.trait), silent = T)
body_size.rare.table <- spc$mean.rare.table
rare.body_size.trait <- spc$body_size.trait

new.w.body_size.abun.table = body_size.abun.table
new.w.body_size.abun.table = as.data.frame(new.w.body_size.abun.table)
for (j in 1:ncol(body_size.abun.table)){
  new.w.body_size.abun.table[, j] <- body_size.abun.table[, j] * abun.body_size.trait
}

new.w.body_size.rare.table = body_size.rare.table
new.w.body_size.rare.table = as.data.frame(new.w.body_size.rare.table)
for (j in 1:ncol(body_size.rare.table)){
  new.w.body_size.rare.table[, j] <- body_size.rare.table[, j] * rare.body_size.trait
}

new.uw.body_size.abun.table = body_size.abun.table
new.uw.body_size.abun.table = as.data.frame(new.uw.body_size.abun.table)
for (j in 1:ncol(body_size.abun.table)){
  new.uw.body_size.abun.table[, j] <- ifelse(body_size.abun.table[, j] > 0,1,0) * abun.body_size.trait
}

new.uw.body_size.rare.table = body_size.rare.table
new.uw.body_size.rare.table = as.data.frame(new.uw.body_size.rare.table)
for (j in 1:ncol(body_size.rare.table)){
  new.uw.body_size.rare.table[, j] <- ifelse(body_size.rare.table[, j] > 0,1,0) * rare.body_size.trait
}

abun.weighted <- colSums(new.w.body_size.abun.table) / colSums(body_size.abun.table)
rare.weighted <- colSums(new.w.body_size.rare.table) / colSums(body_size.rare.table)
abun.unweighted <- colSums(new.uw.body_size.abun.table) / colSums(body_size.abun.table)
rare.unweighted <- colSums(new.uw.body_size.rare.table) / colSums(body_size.rare.table)

# 去除重复的sample
abun.weighted.unique <- abun.weighted[names(abun.weighted) %in% rownames(treat)]
abun.weighted.unique = as.data.frame(abun.weighted.unique)
colnames(abun.weighted.unique) = "value"
abun.weighted.unique$name = rownames(abun.weighted.unique)
abun.weighted.unique$group = "abundant"
abun.weighted.unique$method = "weighted"

rare.weighted.unique <- rare.weighted[names(rare.weighted) %in% rownames(treat)]
rare.weighted.unique = as.data.frame(rare.weighted.unique)
colnames(rare.weighted.unique) = "value"
rare.weighted.unique$name = rownames(rare.weighted.unique)
rare.weighted.unique$group = "rare"
rare.weighted.unique$method = "weighted"

abun.unweighted.unique <- abun.unweighted[names(abun.unweighted) %in% rownames(treat)]
abun.unweighted.unique = as.data.frame(abun.unweighted.unique)
colnames(abun.unweighted.unique) = "value"
abun.unweighted.unique$name = rownames(abun.unweighted.unique)
abun.unweighted.unique$group = "abundant"
abun.unweighted.unique$method = "unweighted"

rare.unweighted.unique <- rare.unweighted[names(rare.unweighted) %in% rownames(treat)]
rare.unweighted.unique = as.data.frame(rare.unweighted.unique)
colnames(rare.unweighted.unique) = "value"
rare.unweighted.unique$name = rownames(rare.unweighted.unique)
rare.unweighted.unique$group = "rare"
rare.unweighted.unique$method = "unweighted"

result.unique = rbind(abun.weighted.unique, rare.weighted.unique,
                      abun.unweighted.unique, rare.unweighted.unique)
body_size.unique = result.unique

result.unique = data.frame(name = shell.unique$name,
                           group = shell.unique$group,
                           method = shell.unique$method,
                           shell = shell.unique$value,
                           locomotion = locomotion.unique$value,
                           bacterivore = bacterivore.unique$value,
                           body_size = body_size.unique$value)
colnames(result.unique)[4:7] = c("Shell","Locomotion","Bacterivore","Body size")

# plot ####
plot.data = result.unique
plot.data = subset(plot.data, plot.data$method == "unweighted")
library(reshape2)
plot.data = melt(
  plot.data,
  id.vars = c("name","group","method"),
  variable.name = "trait",
  value.name = "value"
)

library(ggplot2)
library(gghalves)
plot.data$group <- factor(plot.data$group, levels = c("abundant", "rare"))
ggplot(
  data = plot.data,
  aes(x = group, y = value, fill = group)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.5) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.25,
                    # position = position_dodge(width = 0.76),
                    linewidth = 0.2) +
  geom_half_point_panel(side = "l", shape = 21, size = 2,
                        alpha = 0.7
                        # position = position_dodge(width = 0.78)
                        , color = "white") +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("abundant", "rare")) +
  labs(y = "Value", x = NULL) +
  facet_wrap(~ trait, scales = "free_y", ncol = 2) +
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

ggsave(paste0("abundant vs rare","_","unweighted","_","protist",".pdf"), width = 4.92, height = 4.30, units = "in")

# statistics test ####
plot.data <- result.unique
plot.data <- subset(plot.data, plot.data$method == "unweighted")
plot.data <- plot.data[, -3]

div.table <- data.frame(stringsAsFactors = F)

divindex <- plot.data[, 3:6]
treat.used <- plot.data[, -c(3:6)]
treat.used$group <- ifelse(treat.used$group == "rare", 1, 0)
treat.used$block <- treat$block[match(treat.used$name, rownames(treat))]
treat.used$treat <- treat$combined_treat1[match(treat.used$name, rownames(treat))]
treat.used$layer <- treat$Layer[match(treat.used$name, rownames(treat))]

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
