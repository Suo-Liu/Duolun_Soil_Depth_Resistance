setwd("C:/Users/True/OneDrive/桌面")
rand.remov.once <- function(netRaw, rm.percent, abundance.weighted = F) {
  diag(netRaw) <- 0
  #-随机挑选出一定百分比的OTU
  id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
  net.Raw <- netRaw
  # 这些节点和其他节点连接全部去除
  net.Raw[id.rm, ] <- 0 ## remove all the links to these species
  net.Raw[, id.rm] <- 0
  if (abundance.weighted) {
    # 网络矩阵乘以物种平均丰度，改变相关性值的大小
    net.stength <- net.Raw * sp.ra
  } else {
    net.stength <- net.Raw
  }
  # 每一个节点的平均链接数
  sp.meanInteration <- colMeans(net.stength != 0)

  id.rm2 <- which(sp.meanInteration == 0) ## remove species have no interaction with others
  remain.percent <- (ncol(netRaw) - length(id.rm2)) / ncol(netRaw)
  # for simplicity, I only consider the immediate effects of removing the
  #' id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.

  # you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  remain.percent
}
rmsimu <- function(netRaw.1,netRaw.2,netRaw.3, netRaw.4, rm.p.list,
                   abundance.weighted = F,
                   nperm = 100) {
  t(sapply(rm.p.list, function(x) {
    remains.1 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRaw.1, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remains.2 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRaw.2, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remains.3 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRaw.3, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remains.4 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRaw.4, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    
    a <- t.test(c(remains.1, remains.2), c(remains.3, remains.4))
    t.value <- a[["statistic"]][["t"]]
    p.of.t.value <- a[["p.value"]]

    b <- wilcox.test(c(remains.1, remains.2), c(remains.3, remains.4))
    W.value <- b[["statistic"]][["W"]]
    p.of.wilcox.value <- b[["p.value"]]

    remain.mean.top <- mean(c(remains.1, remains.2))
    remain.mean.sub <- mean(c(remains.3, remains.4))

    remain.sd.top <- sd(c(remains.1, remains.2))
    remain.se.top <- sd(c(remains.1, remains.2)) / (nperm^0.5)

    remain.sd.sub <- sd(c(remains.3, remains.4))
    remain.se.sub <- sd(c(remains.3, remains.4)) / (nperm^0.5)

    result <- c(
      remain.mean.top, remain.mean.sub, remain.sd.top, remain.se.top, remain.sd.sub,
      remain.se.sub, t.value, p.of.t.value, W.value, p.of.wilcox.value
    )
    names(result) <- c(
      "remain.mean.top", "remain.mean.sub", "remain.sd.top", "remain.se.top",
      "remain.sd.sub", "remain.se.sub", "t", "P of t", "W", "P of W"
    )
    result
  }))
}
# top vs sub ####
## each layer ####
library(ieggr)
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\treatment.csv"
prefix1 <- "Pro"
cor.meth <- "pearson"
maj.1 <- 0.25

treat <- read.csv(treat.file, row.names = 1, stringsAsFactors = F)
result <- data.frame(stringsAsFactors = F)

  for (i in unique(treat$plant.type)) {
for (j in unique(treat$Layer)) {
    #     load(paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
    #                 "Protist",
    #     "\\network file\\",
    #                 # "each layer\\",
    #                 # i, "_", j, "_",cor.meth, "_",prefix1,"_",maj.1, "_CLR_FillPairLB.rda"))
    #     i, ".", j, ".",cor.meth, ".",prefix1,".m",maj.1, ".CLR.FillPairLB.rda"))
    # dat = output1$assmc
    
    dat <- read.csv(paste0(
      "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
      "Protist",
      "\\network file\\",
      i, ".", j, ".", cor.meth, ".", prefix1, ".m", maj.1, ".CLR.FillPairLB.csv"
    ), row.names = 1)
    
    dat[is.na(dat)] <- 0
    diag(dat) <- 0
    
    dim(dat)
    dat <- dat[, colSums(dat != 0) > 0]
    dat <- dat[rowSums(dat != 0) > 0, ]
    dim(dat)
    if (j == "L1") {
      dat.L1 <- dat
    } else if (j == "L2") {
      dat.L2 <- dat
    } else if (j == "L3") {
      dat.L3 <- dat
    } else {
      dat.L4 <- dat
    }
  }
  data.simu <- rmsimu(
    netRaw.1 = dat.L1, netRaw.2 = dat.L2, netRaw.3 = dat.L3, netRaw.4 = dat.L4,
    rm.p.list = seq(0.2, 0.85, by = 0.05), abundance.weighted = F, nperm = 100
  )
  data.simu <- data.frame(
    Proportion.removed = seq(0.2, 0.85, by = 0.05),
    data.simu
  )
  data.simu$steppes <- i
  result <- rbind(result, data.simu)
}
write.csv(result, paste0(prefix1, "_", cor.meth, "_", maj.1, "_robustness.csv"))

## combine two soil layers ####
rand.remov.once <- function(netRaw, rm.percent, abundance.weighted = F) {
  diag(netRaw) <- 0
  #-随机挑选出一定百分比的OTU
  id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
  net.Raw <- netRaw
  # 这些节点和其他节点连接全部去除
  net.Raw[id.rm, ] <- 0 ## remove all the links to these species
  net.Raw[, id.rm] <- 0
  if (abundance.weighted) {
    # 网络矩阵乘以物种平均丰度，改变相关性值的大小
    net.stength <- net.Raw * sp.ra
  } else {
    net.stength <- net.Raw
  }
  # 每一个节点的平均链接数
  sp.meanInteration <- colMeans(net.stength != 0)
  
  id.rm2 <- which(sp.meanInteration == 0) ## remove species have no interaction with others
  remain.percent <- (ncol(netRaw) - length(id.rm2)) / ncol(netRaw)
  # for simplicity, I only consider the immediate effects of removing the
  #' id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  # you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  remain.percent
}
rmsimu <- function(netRawTr.1, netRawTr.2,
                   netRawCt.1, netRawCt.2,
                   rm.p.list,
                   abundance.weighted = F,
                   nperm = 100) {
  t(sapply(rm.p.list, function(x) {
    remainsTr.1 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRawTr.1, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remainsTr.2 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRawTr.2, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remainsCt.1 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRawCt.1, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remainsCt.2 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRawCt.2, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    
    a <- t.test(c(remainsTr.1, remainsTr.2), c(remainsCt.1, remainsCt.2))
    t.value <- a[["statistic"]][["t"]]
    p.of.t.value <- a[["p.value"]]
    
    b <- wilcox.test(c(remainsTr.1, remainsTr.2), c(remainsCt.1, remainsCt.2))
    W.value <- b[["statistic"]][["W"]]
    p.of.wilcox.value <- b[["p.value"]]
    
    remain.mean.Tr <- mean(c(remainsTr.1, remainsTr.2))
    remain.mean.Ct <- mean(c(remainsCt.1, remainsCt.2))
    
    remain.sd.Tr <- sd(c(remainsTr.1, remainsTr.2))
    remain.se.Tr <- sd(c(remainsTr.1, remainsTr.2)) / (nperm^0.5)
    
    remain.sd.Ct <- sd(c(remainsCt.1, remainsCt.2))
    remain.se.Ct <- sd(c(remainsCt.1, remainsCt.2)) / (nperm^0.5)
    
    result <- c(
      remain.mean.Tr, remain.mean.Ct, remain.sd.Tr, remain.se.Tr, remain.sd.Ct,
      remain.se.Ct, t.value, p.of.t.value, W.value, p.of.wilcox.value
    )
    names(result) <- c(
      "remain.mean.treat", "remain.mean.ctrl", "remain.sd.treat", "remain.se.treat",
      "remain.sd.ctrl", "remain.se.ctrl", "t", "P of t", "W", "P of W"
    )
    result
  }))
}

library(ieggr)
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\treatment.csv"
prefix1 <- "Fungi"
cor.meth <- "pearson"
maj.1 <- 0.25

treat <- read.csv(treat.file, row.names = 1, stringsAsFactors = F)
result <- data.frame(stringsAsFactors = F)

for (i in unique(treat$plant.type)) {
  load(paste0(
    "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
    "Fungi",
    "\\network file\\",
    # "each layer\\",
    i, "_", "L3", "_", cor.meth, "_", prefix1, "_", maj.1, "_CLR_FillPairLB.rda"
  ))
  dat.1 <- output1$assmc
  
  load(paste0(
    "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
    "Fungi",
    "\\network file\\",
    # "each layer\\",
    i, "_", "L4", "_", cor.meth, "_", prefix1, "_", maj.1, "_CLR_FillPairLB.rda"
  ))
  dat.2 <- output1$assmc
  
  # dat.1 = read.csv(paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
  #                       "Protist",
  #                       "\\network file\\",
  #                       i, ".", "L3", ".",cor.meth, ".",prefix1,".m",maj.1, ".CLR.FillPairLB.csv"), row.names = 1)
  #
  # dat.2 = read.csv(paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
  #                       "Protist",
  #                       "\\network file\\",
  #                       i, ".", "L4", ".",cor.meth, ".",prefix1,".m",maj.1, ".CLR.FillPairLB.csv"), row.names = 1)
  
  dat.1[is.na(dat.1)] <- 0
  diag(dat.1) <- 0
  
  dat.2[is.na(dat.2)] <- 0
  diag(dat.2) <- 0
  
  dim(dat.1)
  dat.1 <- dat.1[, colSums(dat.1 != 0) > 0]
  dat.1 <- dat.1[rowSums(dat.1 != 0) > 0, ]
  dim(dat.1)
  
  dim(dat.2)
  dat.2 <- dat.2[, colSums(dat.2 != 0) > 0]
  dat.2 <- dat.2[rowSums(dat.2 != 0) > 0, ]
  dim(dat.2)
  
  if (i == "MS") {
    dat.1.MS <- dat.1
    dat.2.MS <- dat.2
  } else if (i == "TS") {
    dat.1.TS <- dat.1
    dat.2.TS <- dat.2
  } else if (i == "DS") {
    dat.1.DS <- dat.1
    dat.2.DS <- dat.2
  }
}
data.simu.MT <- rmsimu(
  netRawTr.1 = dat.1.MS, netRawTr.2 = dat.2.MS,
  netRawCt.1 = dat.1.TS, netRawCt.2 = dat.2.TS,
  rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
)
data.simu.MD <- rmsimu(
  netRawTr.1 = dat.1.MS, netRawTr.2 = dat.2.MS,
  netRawCt.1 = dat.1.DS, netRawCt.2 = dat.2.DS,
  rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
)
data.simu.TD <- rmsimu(
  netRawTr.1 = dat.1.TS, netRawTr.2 = dat.2.TS,
  netRawCt.1 = dat.1.DS, netRawCt.2 = dat.2.DS,
  rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
)

data.simu.MT <- data.frame(
  Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("M.vs.T", 18),
  data.simu.MT
)
data.simu.MD <- data.frame(
  Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("M.vs.D", 18),
  data.simu.MD
)
data.simu.TD <- data.frame(
  Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("T.vs.D", 18),
  data.simu.TD
)

data.simu <- rbind(data.simu.MT, data.simu.MD, data.simu.TD)
data.simu$layer <- "L34"
result <- rbind(result, data.simu)

write.csv(result, paste0(prefix1, "_", cor.meth, "_", maj.1, "_L34_robustness.csv"))


# steppes ####
## each layer ####
library(ieggr)
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\treatment.csv"
prefix1 <- "Pro"
cor.meth <- "pearson"
maj.1 <- 0.25

treat <- read.csv(treat.file, row.names = 1, stringsAsFactors = F)
result <- data.frame(stringsAsFactors = F)

for (j in unique(treat$Layer)) {
  for (i in unique(treat$plant.type)) {
    #     load(paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
    #                 "Protist",
    #     "\\network file\\",
    #                 # "each layer\\",
    #                 # i, "_", j, "_",cor.meth, "_",prefix1,"_",maj.1, "_CLR_FillPairLB.rda"))
    #     i, ".", j, ".",cor.meth, ".",prefix1,".m",maj.1, ".CLR.FillPairLB.rda"))
    # dat = output1$assmc

    dat <- read.csv(paste0(
      "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
      "Protist",
      "\\network file\\",
      i, ".", j, ".", cor.meth, ".", prefix1, ".m", maj.1, ".CLR.FillPairLB.csv"
    ), row.names = 1)

    dat[is.na(dat)] <- 0
    diag(dat) <- 0

    dim(dat)
    dat <- dat[, colSums(dat != 0) > 0]
    dat <- dat[rowSums(dat != 0) > 0, ]
    dim(dat)
    if (i == "MS") {
      dat.MS <- dat
    } else if (i == "TS") {
      dat.TS <- dat
    } else if (i == "DS") {
      dat.DS <- dat
    }
  }
  data.simu.MT <- rmsimu(
    netRawTr = dat.MS, netRawCt = dat.TS,
    rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
  )
  data.simu.MD <- rmsimu(
    netRawTr = dat.MS, netRawCt = dat.DS,
    rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
  )
  data.simu.TD <- rmsimu(
    netRawTr = dat.TS, netRawCt = dat.DS,
    rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
  )

  data.simu.MT <- data.frame(
    Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("M.vs.T", 18),
    data.simu.MT
  )
  data.simu.MD <- data.frame(
    Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("M.vs.D", 18),
    data.simu.MD
  )
  data.simu.TD <- data.frame(
    Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("T.vs.D", 18),
    data.simu.TD
  )

  data.simu <- rbind(data.simu.MT, data.simu.MD, data.simu.TD)
  data.simu$layer <- j
  result <- rbind(result, data.simu)
}
write.csv(result, paste0(prefix1, "_", cor.meth, "_", maj.1, "_robustness.csv"))

## combine two soil layers ####
rand.remov.once <- function(netRaw, rm.percent, abundance.weighted = F) {
  diag(netRaw) <- 0
  #-随机挑选出一定百分比的OTU
  id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
  net.Raw <- netRaw
  # 这些节点和其他节点连接全部去除
  net.Raw[id.rm, ] <- 0 ## remove all the links to these species
  net.Raw[, id.rm] <- 0
  if (abundance.weighted) {
    # 网络矩阵乘以物种平均丰度，改变相关性值的大小
    net.stength <- net.Raw * sp.ra
  } else {
    net.stength <- net.Raw
  }
  # 每一个节点的平均链接数
  sp.meanInteration <- colMeans(net.stength != 0)

  id.rm2 <- which(sp.meanInteration == 0) ## remove species have no interaction with others
  remain.percent <- (ncol(netRaw) - length(id.rm2)) / ncol(netRaw)
  # for simplicity, I only consider the immediate effects of removing the
  #' id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.

  # you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  remain.percent
}
rmsimu <- function(netRawTr.1, netRawTr.2,
                   netRawCt.1, netRawCt.2,
                   rm.p.list,
                   abundance.weighted = F,
                   nperm = 100) {
  t(sapply(rm.p.list, function(x) {
    remainsTr.1 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRawTr.1, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remainsTr.2 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRawTr.2, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remainsCt.1 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRawCt.1, rm.percent = x, abundance.weighted = abundance.weighted)
    })
    remainsCt.2 <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRawCt.2, rm.percent = x, abundance.weighted = abundance.weighted)
    })

    a <- t.test(c(remainsTr.1, remainsTr.2), c(remainsCt.1, remainsCt.2))
    t.value <- a[["statistic"]][["t"]]
    p.of.t.value <- a[["p.value"]]

    b <- wilcox.test(c(remainsTr.1, remainsTr.2), c(remainsCt.1, remainsCt.2))
    W.value <- b[["statistic"]][["W"]]
    p.of.wilcox.value <- b[["p.value"]]

    remain.mean.Tr <- mean(c(remainsTr.1, remainsTr.2))
    remain.mean.Ct <- mean(c(remainsCt.1, remainsCt.2))

    remain.sd.Tr <- sd(c(remainsTr.1, remainsTr.2))
    remain.se.Tr <- sd(c(remainsTr.1, remainsTr.2)) / (nperm^0.5)

    remain.sd.Ct <- sd(c(remainsCt.1, remainsCt.2))
    remain.se.Ct <- sd(c(remainsCt.1, remainsCt.2)) / (nperm^0.5)

    result <- c(
      remain.mean.Tr, remain.mean.Ct, remain.sd.Tr, remain.se.Tr, remain.sd.Ct,
      remain.se.Ct, t.value, p.of.t.value, W.value, p.of.wilcox.value
    )
    names(result) <- c(
      "remain.mean.treat", "remain.mean.ctrl", "remain.sd.treat", "remain.se.treat",
      "remain.sd.ctrl", "remain.se.ctrl", "t", "P of t", "W", "P of W"
    )
    result
  }))
}

library(ieggr)
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\treatment.csv"
prefix1 <- "Fungi"
cor.meth <- "pearson"
maj.1 <- 0.25

treat <- read.csv(treat.file, row.names = 1, stringsAsFactors = F)
result <- data.frame(stringsAsFactors = F)

for (i in unique(treat$plant.type)) {
  load(paste0(
    "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
    "Fungi",
    "\\network file\\",
    # "each layer\\",
    i, "_", "L3", "_", cor.meth, "_", prefix1, "_", maj.1, "_CLR_FillPairLB.rda"
  ))
  dat.1 <- output1$assmc

  load(paste0(
    "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
    "Fungi",
    "\\network file\\",
    # "each layer\\",
    i, "_", "L4", "_", cor.meth, "_", prefix1, "_", maj.1, "_CLR_FillPairLB.rda"
  ))
  dat.2 <- output1$assmc

  # dat.1 = read.csv(paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
  #                       "Protist",
  #                       "\\network file\\",
  #                       i, ".", "L3", ".",cor.meth, ".",prefix1,".m",maj.1, ".CLR.FillPairLB.csv"), row.names = 1)
  #
  # dat.2 = read.csv(paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\network\\",
  #                       "Protist",
  #                       "\\network file\\",
  #                       i, ".", "L4", ".",cor.meth, ".",prefix1,".m",maj.1, ".CLR.FillPairLB.csv"), row.names = 1)

  dat.1[is.na(dat.1)] <- 0
  diag(dat.1) <- 0

  dat.2[is.na(dat.2)] <- 0
  diag(dat.2) <- 0

  dim(dat.1)
  dat.1 <- dat.1[, colSums(dat.1 != 0) > 0]
  dat.1 <- dat.1[rowSums(dat.1 != 0) > 0, ]
  dim(dat.1)

  dim(dat.2)
  dat.2 <- dat.2[, colSums(dat.2 != 0) > 0]
  dat.2 <- dat.2[rowSums(dat.2 != 0) > 0, ]
  dim(dat.2)

  if (i == "MS") {
    dat.1.MS <- dat.1
    dat.2.MS <- dat.2
  } else if (i == "TS") {
    dat.1.TS <- dat.1
    dat.2.TS <- dat.2
  } else if (i == "DS") {
    dat.1.DS <- dat.1
    dat.2.DS <- dat.2
  }
}
data.simu.MT <- rmsimu(
  netRawTr.1 = dat.1.MS, netRawTr.2 = dat.2.MS,
  netRawCt.1 = dat.1.TS, netRawCt.2 = dat.2.TS,
  rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
)
data.simu.MD <- rmsimu(
  netRawTr.1 = dat.1.MS, netRawTr.2 = dat.2.MS,
  netRawCt.1 = dat.1.DS, netRawCt.2 = dat.2.DS,
  rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
)
data.simu.TD <- rmsimu(
  netRawTr.1 = dat.1.TS, netRawTr.2 = dat.2.TS,
  netRawCt.1 = dat.1.DS, netRawCt.2 = dat.2.DS,
  rm.p.list = seq(0.1, 0.95, by = 0.05), abundance.weighted = F, nperm = 100
)

data.simu.MT <- data.frame(
  Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("M.vs.T", 18),
  data.simu.MT
)
data.simu.MD <- data.frame(
  Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("M.vs.D", 18),
  data.simu.MD
)
data.simu.TD <- data.frame(
  Proportion.removed = seq(0.1, 0.95, by = 0.05), group = rep("T.vs.D", 18),
  data.simu.TD
)

data.simu <- rbind(data.simu.MT, data.simu.MD, data.simu.TD)
data.simu$layer <- "L34"
result <- rbind(result, data.simu)

write.csv(result, paste0(prefix1, "_", cor.meth, "_", maj.1, "_L34_robustness.csv"))
