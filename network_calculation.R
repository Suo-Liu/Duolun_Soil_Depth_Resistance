setwd("C:\\Users\\True\\OneDrive\\桌面")
# global network calculation ####
library(ieggr)
require(MENA)

# setwd("/home/fangke1/Duolun/RMT_network/Bacteria")
# 
# source("/home/fangke1/Duolun/RMT_network/findcutoff2.r")
# 
# treat.file <- "/home/fangke1/Duolun/RMT_network/data for use/treatment.csv"
# 
# com.file <- "/home/fangke1/Duolun/RMT_network/data for use/protist_zotu_2868.txt"
# prefix.m <- "Protist"
# prefix1 <- "Pro"
# 
# com.file <- "/home/fangke1/Duolun/RMT_network/data for use/fungi_zotu_1174.txt"
# prefix.m <- "Fungi18S"
# prefix1 <- "Fungi18S"
# 
# com.file <- "/home/fangke1/Duolun/RMT_network/data for use/nematoda_zotu_259.txt"
# prefix.m <- "Nematode"
# prefix1 <- "Nem"
# 
# com.file <- "/home/fangke1/Duolun/RMT_network/data for use/18S_zotu_6520.txt"
# prefix.m <- "18S"
# prefix1 <- "X18S"

# source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\findcutoff2.r")
# treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"

# com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
# prefix.m <- "Bacteria"
# prefix1 <- "Bac"

# com.file = "/vd04/home/liusuo/inner_mongolia/bacteria.zotu245_resample_35000.txt"
# prefix.m = "Bacteria"
# prefix1 = "Bac"

# com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\ITS\\all samples\\Unoise\\zotutab220_max_2_resample_10000.txt"
# prefix.m <- "Fungi"
# prefix1 <- "Fungi"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"
prefix.m <- "Protist"
prefix1 <- "Pro"

# com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\nematode_resample_248.txt"
# prefix.m <- "Nematode"
# prefix1 <- "Nem"

comm <- t(read.table(com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))
name.row <- rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] <- "YD3L3"
rownames(comm) <- name.row

cor.meth <- "pearson"
maj.1 <- 0.25

prefixm <- cor.meth

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
## find cutoff ####
for (i in c("MS", "TS", "DS")) {
  for (j in c("L1", "L2", "L3", "L4")) {
    used.comm <- comm
    used.treat <- subset(treat, plant.type == i & Layer == j)
    match.rowname <- match.name(rn.list = list(used.comm = used.comm, used.treat = used.treat))
    used.comm <- match.rowname$used.comm
    used.treat <- match.rowname$used.treat

    # filter abundance by majority
    {
      dim(used.comm)
      used.comm <- used.comm[, which((colSums(used.comm > 0) / dim(used.comm)[1]) >= maj.1)]
      dim(used.comm)
    }
    colnames(used.comm) <- paste0(prefix1, colnames(used.comm))

    sum(colSums(used.comm, na.rm = TRUE) == 0)

    assmx1 <- MENA::assmatrix(used.comm,
      method = cor.meth, majority = maj.1,
      missing.data.fill = "fill.pair.BL", logarithm = TRUE,
      fillzero.value = 0.01, samp.time = NULL, time.lag = 0,
      silent = FALSE, mklthread = 4, output.filled.matrix = FALSE,
      CLR.transform = TRUE
    )

    brmt1 <- MENA::brodyrmt(ass.matrix = assmx1, nthread = 50)

    cutoff1out <- findcutoff2(brmt1, outputxy = TRUE, criteria = c(0.1, 0.5, 0.9))
    cutoff1 <- cutoff1out$cutoff
    smthxy1 <- cutoff1out$smth_xy

    cutoff.use1 <- cutoff1[1, 1]
    assmc1 <- MENA::assmcut(ass.matrix = assmx1, cutoff = cutoff.use1)

    output1 <- list(assmx = assmx1, brmt = brmt1, cutoff = cutoff1, smthxy = smthxy1, cutoff.use = cutoff.use1, assmc = assmc1)
    save(output1, file = paste0(i, ".", j, ".", prefixm, ".", prefix1, ".m", maj.1, ".CLR.FillPairLB.rda"))
  }
}
## after find cutoff ####
for (i in c("MS", "TS", "DS")) {
  for (j in c("L1", "L2", "L3", "L4")) {
    used.comm <- comm
    used.treat <- subset(treat, plant.type == i & Layer == j)
    match.rowname <- match.name(rn.list = list(used.comm = used.comm, used.treat = used.treat))
    used.comm <- match.rowname$used.comm
    used.treat <- match.rowname$used.treat

    # filter abundance by majority
    {
      dim(used.comm)
      used.comm <- used.comm[, which((colSums(used.comm > 0) / dim(used.comm)[1]) >= maj.1)]
      dim(used.comm)
    }
    colnames(used.comm) <- paste0(prefix1, colnames(used.comm))

    assmx1 <- MENA::assmatrix(used.comm,
      method = cor.meth, majority = maj.1,
      missing.data.fill = "fill.pair.BL", logarithm = TRUE,
      fillzero.value = 0.01, samp.time = NULL, time.lag = 0,
      silent = FALSE, mklthread = 4, output.filled.matrix = FALSE,
      CLR.transform = TRUE
    )

    cutoff.use1 <- 0.911

    assmc1 <- MENA::assmcut(ass.matrix = assmx1, cutoff = cutoff.use1)
    write.csv(assmc1, file = paste0(i, ".", j, ".", prefixm, ".", prefix1, ".m", maj.1, ".CLR.FillPairLB.csv"))
  }
}
# small networks: one sample as one network ####
library(ieggr)

source("/home/fangke1/Duolun/RMT_network/newnetindex.R")
treat.file <- "/home/fangke1/Duolun/RMT_network/data for use/treatment.csv"

cor.meth <- "pearson"
maj.1 <- 0.25

prefixm <- cor.meth

result.list <- list()
z <- 1
pre.comm <- t(read.table(com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))
name.row <- rownames(pre.comm)
name.row[which(rownames(pre.comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pre.comm) == "YD52L3")] <- "YD3L3"
rownames(pre.comm) <- name.row
# pre.comm <- pre.comm[!rownames(pre.comm) %in% c("YD37L3", "YD38L3","YD39L3","YD52L3","YD53L3","YD54L3"), ]

pre.treat <- read.csv(treat.file, header = TRUE, sep = ",", row.names = 1)
for (i in c("MS", "TS", "DS")) {
  for (j in c("L1", "L2", "L3", "L4")) {
    comm <- pre.comm
    treat <- subset(pre.treat, plant.type == i & Layer == j)

    match.rowname <- match.name(rn.list = list(comm = comm, treat = treat), silent = T)
    comm <- match.rowname$comm
    treat <- match.rowname$treat

    for (k in rownames(treat)) {
      dat <- read.csv(paste0(i, ".", j, ".", prefixm, ".", 
                             prefix1, ".m", maj.1, ".CLR.FillPairLB.csv"), row.names = 1,
                      header = TRUE, sep = ",")

      used.comm <- comm[rownames(comm) == k, , drop = F]
      used.comm <- used.comm[, colSums(used.comm) > 0, drop = F]
      colnames(used.comm) <- paste0(prefix1, colnames(used.comm))

      save.spec <- colnames(used.comm)

      dat <- dat[rownames(dat) %in% save.spec, colnames(dat) %in% save.spec]
      dat <- as.matrix(dat)
      dat[is.na(dat)] <- 0

      # remove isolated node
      dim(dat)
      diag(dat) <- 0
      dat <- dat[, colSums(abs(dat)) > 0]
      dat <- dat[rowSums(abs(dat)) > 0, ]
      dim(dat)

      result.list[[z]] <- dat
      names(result.list)[z] <- k
      z <- z + 1
    }
  }
}
saveRDS(result.list, paste0(prefix1, ".",cor.meth,".", maj.1, ".RDS"))

# network properties ####
library(tidyr)
require(MENA)
net.att <- list()
module.att <- list()
for (z in 1:length(result.list)) {
  cor.data <- result.list[[z]]
  message("-----Now z=", z, " in ", length(result.list), ". ", date())
  netind <- newnetindex(cor.data, fast = F)
  net.att[[z]] <- data.frame(Group = names(result.list)[z], netind$network, stringsAsFactors = FALSE)
  modules <- module(assmc = cor.data, methods = "greedy")
  modules$module$sample <- names(result.list)[z]
  module.att[[z]] <- data.frame(modules$module, stringsAsFactors = FALSE)
}
net.com.index <- data.frame(stringsAsFactors = F)
module.com.index <- data.frame(stringsAsFactors = F)

for (i in 1:length(net.att)) {
  net.com.index <- rbind(net.com.index, net.att[[i]])
  module.com.index <- rbind(module.com.index, module.att[[i]])
}
spr.net.com.index <- spread(net.com.index, NetworkIndex, Value)
rownames(spr.net.com.index) <- spr.net.com.index$Group
rownames(module.com.index) <- module.com.index$sample
module.com.index <- module.com.index[match(rownames(spr.net.com.index), rownames(module.com.index)), ]
spr.net.com.index <- cbind(spr.net.com.index, module.com.index[, 2:3])
spr.net.com.index <- spr.net.com.index[, -1]
write.csv(spr.net.com.index, paste0(prefix1, ".", cor.meth,".",maj.1, ".csv"))

# small networks: all samples in layer&plant as one network ####
library(ieggr)
wd <- iwd("C:\\Users\\True\\OneDrive\\桌面")
save.wd <- iwd("C:\\Users\\True\\OneDrive\\桌面")

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
prefix.m <- "Bacteria"
prefix1 <- "Bac"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\ITS\\all samples\\Unoise\\zotutab220_max_2_resample_10000.txt"
prefix.m <- "Fungi"
prefix1 <- "Fungi"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\protist.zotu_resample_2961.txt"
clas.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_tax.txt"
prefix.m <- "Protist"
prefix1 <- "Pro"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\nematode_resample_248.txt"
clas.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_tax.txt"
prefix.m <- "Nematode"
prefix1 <- "Nem"

cor.meth <- "spearman"
maj.1 <- 0.75

prefixm <- paste0(cor.meth)

result.list <- list()
z <- 1
pre.comm <- t(read.table(com.file,
                         header = TRUE, sep = "\t", row.names = 1,
                         as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                         check.names = FALSE
))
name.row <- rownames(pre.comm)
name.row[which(rownames(pre.comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pre.comm) == "YD52L3")] <- "YD3L3"
rownames(pre.comm) <- name.row
pre.comm <- pre.comm[rownames(pre.comm) != "YD52L3", ]

pre.treat <- read.csv(treat.file, header = TRUE, sep = ",", row.names = 1)
for (i in c("MS", "TS", "DS")) {
  for (j in c("L1", "L2", "L3", "L4")) {
    load(paste0(i, ".", j, ".", prefixm, ".", prefix1, ".m", maj.1, ".CLR.FillPairLB.rda"))
    comm <- pre.comm
    treat <- subset(pre.treat, plant.type == i & Layer == j)
    
    match.rowname <- match.name(rn.list = list(comm = comm, treat = treat))
    comm <- match.rowname$comm
    treat <- match.rowname$treat
    
    for (k in unique(treat$combined_treat1)) {
      dat <- output1$assmc
      used.treat = subset(treat, combined_treat1 == k)
      match.rowname <- match.name(rn.list = list(comm = comm, used.treat = used.treat))
      used.comm <- match.rowname$comm
      used.treat <- match.rowname$used.treat
      
      used.comm <- used.comm[, colSums(used.comm) > 0, drop = F]
      colnames(used.comm) <- paste0(prefix1, colnames(used.comm))
      
      save.spec <- colnames(used.comm)
      
      dat <- dat[rownames(dat) %in% save.spec, colnames(dat) %in% save.spec]
      dat <- as.matrix(dat)
      dat[is.na(dat)] <- 0
      
      # remove isolated node
      dim(dat)
      diag(dat) <- 0
      dat <- dat[, colSums(abs(dat)) > 0]
      dat <- dat[rowSums(abs(dat)) > 0, ]
      dim(dat)
      
      result.list[[z]] <- dat
      names(result.list)[z] <- paste0(i,".",j,".",k)
      z <- z + 1
    }
  }
}
saveRDS(result.list, paste0("3to1",".",prefix1, ".", maj.1, ".", "BTS", ".RDS"))

source("C:/Users/True/OneDrive/桌面/research/analysis methods/R.code/newnetindex.R")
result.list <- result.list
library(tidyr)
require(MENA)
net.att <- list()
module.att <- list()
for (z in 1:length(result.list)) {
  cor.data <- result.list[[z]]
  message("-----Now z=", z, " in ", length(result.list), ". ", date())
  netind <- newnetindex(cor.data, fast = F)
  net.att[[z]] <- data.frame(Group = names(result.list)[z], netind$network, stringsAsFactors = FALSE)
  modules <- module(assmc = cor.data, methods = "greedy")
  modules$module$sample <- names(result.list)[z]
  module.att[[z]] <- data.frame(modules$module, stringsAsFactors = FALSE)
}
net.com.index <- data.frame(stringsAsFactors = F)
module.com.index <- data.frame(stringsAsFactors = F)

for (i in 1:length(net.att)) {
  net.com.index <- rbind(net.com.index, net.att[[i]])
  module.com.index <- rbind(module.com.index, module.att[[i]])
}
spr.net.com.index <- spread(net.com.index, NetworkIndex, Value)
rownames(spr.net.com.index) <- spr.net.com.index$Group
rownames(module.com.index) <- module.com.index$sample
module.com.index <- module.com.index[match(rownames(spr.net.com.index), rownames(module.com.index)), ]
spr.net.com.index <- cbind(spr.net.com.index, module.com.index[, 2:3])
spr.net.com.index <- spr.net.com.index[, -1]
write.csv(spr.net.com.index, paste0(prefix1, ".", maj.1, ".", "BTS", ".csv"))
# LMM ####
divindex.file <- file.choose()
prefix1 <- "Pro"
cor.meth <- "spearman"
maj.1 <- 0.5

prefix <- paste0(prefix1, ".", maj.1, ".", cor.meth)

divindex <- read.csv(divindex.file, row.names = 1)

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"

treat <- read.csv(treat.file, row.names = 1)
treat = subset(treat, Layer %in% c("L1") & plant.type == "TS")
library(ieggr)
sampc=match.name(rn.list = list(treat = treat,divindex = divindex))
divindex = sampc$divindex
treat = sampc$treat

library(ieggr)
library(lme4)
library(car)
library(tidyr)
library(dplyr)
library(magrittr)
divindex<-scale(divindex)
divs1<-sapply(1:ncol(divindex),function(j){
  message("Now j=",j," in ",ncol(divindex),". ",date())
  if (length(unique(divindex[,j]))<3){
    result<-rep(NA,38)
  } else {
    div<-data.frame(divtest=divindex[,j],treat)
    fm<-lmer(divtest~warm*precip+(1|block),data=div)
    presult<-car::Anova(fm,type=2)
    coefs<-coef(summary(fm))[ , "Estimate"]#four coefs
    names(coefs)<-paste0(names(coefs),".mean")
    SEvalues<-coef(summary(fm))[ , "Std. Error"]#standard errors
    names(SEvalues)<-paste0(names(SEvalues),".se")
    tvalues<-coef(summary(fm))[ , "t value"]#t values
    names(tvalues)<-paste0(names(tvalues),".t")
    chisqP<-c(presult[,1],presult[,3])
    names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
    chisqP<-c(presult[,3])
    names(chisqP)<-c(paste0(row.names(presult),".P"))
    result<-c(coefs,tvalues,SEvalues,chisqP)
    # result<-c(coefs,SEvalues,chisqP)
  }
})
colnames(divs1)<-colnames(divindex)

divs1 = divs1[!grepl("Intercept",rownames(divs1)),]
divs1 = as.data.frame(divs1)
divs1$treat = rownames(divs1)
divs1 = separate(divs1,col = treat,sep = '\\.',into = c('treat','type'))
divs1 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)
write.csv(divs1, paste("p.value.", "LMM.", prefix, ".csv", sep = ""))

# node role ####
source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\newnetindex.R")
library(MENA)
prefix1 <- "Pro"
cor.meth <- "spearman"
maj.1 <- 0.5
prefixm <- paste0(cor.meth,".",prefix1, ".","m",maj.1)
net.att <- list()
module.att <- list()
node.att <- list()
z = 1
for (i in c("MS", "TS", "DS")) {
  for (j in c("L1", "L2", "L3", "L4")) {
network.used <- read.csv(paste0(i, ".", j, ".", prefixm, ".CLR.FillPairLB.csv"), 
                         row.names = 1, header = T, sep = ",")

  cor.data <- network.used
  cor.data <- as.matrix(cor.data)
  cor.data[is.na(cor.data)] <- 0

  # remove isolated node
  diag(cor.data) = 0
  cor.data <- cor.data[rowSums(abs(cor.data)) > 0, ]
  cor.data <- cor.data[, colSums(abs(cor.data)) > 0]

  netind <- newnetindex(cor.data, fast = F)
  net.att[[z]] <- data.frame(Group = paste0(i,".",j), netind$network, stringsAsFactors = FALSE)
  modules <- module(assmc = cor.data, methods = "greedy")
  node.module <- modules$node

  node.att[[z]] <- node.module
  names(node.att)[[z]] <- paste0(i,".",j)
  modules$module$sample <- paste0(i,".",j)
  module.att[[z]] <- data.frame(modules$module, stringsAsFactors = FALSE)
  z = z+1
  }
}

net.com.index <- data.frame(stringsAsFactors = F)
module.com.index <- data.frame(stringsAsFactors = F)

for (i in 1:length(net.att)) {
  net.com.index <- rbind(net.com.index, net.att[[i]])
  module.com.index <- rbind(module.com.index, module.att[[i]])
}
library(tidyr)
spr.net.com.index <- spread(net.com.index, NetworkIndex, Value)
rownames(spr.net.com.index) <- spr.net.com.index$Group
rownames(module.com.index) <- module.com.index$sample

write.csv(spr.net.com.index, paste0(prefixm, ".", "network index", ".", prefix2, ".", maj.2, ".csv"))
write.csv(module.com.index, paste0(prefixm, ".", "module index", ".", prefix2, ".", maj.2, ".csv"))

saveRDS(node.att, paste0(prefixm, ".", "node index", ".RDS"))

# node roles of abundant and rare species ####
node.list = readRDS("C:\\Users\\True\\OneDrive\\桌面\\spearman.Pro.m0.5.node index.Pro.0.5.RDS")
library(dplyr)
merged_df <- bind_rows(lapply(names(node.list), function(name) {
  df <- node.list[[name]]
  df$name_column <- name  # 添加新列，值为dataframe名称
  return(df)
}))
library(tidyr)
merged_df = separate(merged_df,name_column, into = c("plant", "layer"), sep = "\\.")
merged_df = subset(merged_df, layer %in% c("L1","L2"))
table(merged_df$Classification)

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\treatment.csv"
com.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\protist.zotu_resample_2961.txt"
tree.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_rooted_tree.nwk"
clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\18S_tax.txt"
prefix.m = "Protist"

comm=t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                  check.names = FALSE))
name.row = rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] = "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] = "YD3L3"
rownames(comm) = name.row

pre.treat = read.csv(treat.file, header = TRUE, sep = ",", row.names = 1)
pre.treat <- subset(pre.treat, Layer %in% c("L1","L2"))
prefix.m = paste0(prefix.m, ".","L12")

treat = pre.treat
library(ieggr)
sampc=match.name(rn.list = list(treat = treat,comm = comm))
dim(comm)
comm = sampc$comm
treat = sampc$treat
comm = comm[,colSums(comm)>0]
dim(comm)
## abundant and rare species ####
# rare species由于其自身定义原因，很难在network出现，而仅有的几个只是普通的node
# 统计物种出现频率
library(dplyr)
fre.comm = comm
fre.comm = as.data.frame(fre.comm)
fre.comm$plant.type = treat$plant.type
fre.comm %<>% group_by(plant.type) %>%
  summarise_all(sum)
fre.comm = as.data.frame(fre.comm)
rownames(fre.comm) = fre.comm[,1]
fre.comm = fre.comm[,-1]
# 识别交集中的物种和各草原类型独有物种
inter.otu = names(fre.comm)[which(colSums(fre.comm>0) == nrow(fre.comm))]
DS.otu = names(fre.comm)[which(fre.comm[1,] != 0 & fre.comm[2,] == 0 &fre.comm[3,] == 0)]
MS.otu = names(fre.comm)[which(fre.comm[2,] != 0 & fre.comm[1,] == 0 &fre.comm[3,] == 0)]
TS.otu = names(fre.comm)[which(fre.comm[3,] != 0 & fre.comm[1,] == 0 &fre.comm[2,] == 0)]
# 识别丰富物种和稀有物种，定义源于Dai et al., NC, 2022，但是丰富物种的定义少了出现频率的限制，因为过多限制会导致最后识别的丰富物种太少
# abun.otu = which(colSums(comm/rowSums(comm) >= 0.001) > 0.5 & (colSums(comm>0)>(0.5*nrow(comm))))
abun.otu = which(colSums(comm/rowSums(comm) >= 0.001) > 0.5)
abun.otu = abun.otu[which(names(abun.otu)%in%inter.otu)]

rare.otu = which(colSums(comm/rowSums(comm) < 0.001) == nrow(comm))
rare.otu.DS = rare.otu[which(names(rare.otu)%in%DS.otu)]
rare.otu.MS = rare.otu[which(names(rare.otu)%in%MS.otu)]
rare.otu.TS = rare.otu[which(names(rare.otu)%in%TS.otu)]

for (i in 1:nrow(merged_df)){
  comm = comm[-which(rownames(comm) == names(rare.otu.DS)[i]),]
}

abun.node.df = merged_df[merged_df$ID%in%paste0("Pro",names(abun.otu)),]
rare.TS.df = merged_df[merged_df$ID%in%paste0("Pro",names(rare.otu.TS)),]
rare.MS.df = merged_df[merged_df$ID%in%paste0("Pro",names(rare.otu.MS)),]
rare.DS.df = merged_df[merged_df$ID%in%paste0("Pro",names(rare.otu.DS)),]

which(paste0("Pro",names(abun.otu)))

names(rare.otu.DS)
names(rare.otu.MS)
names(rare.otu.TS)
