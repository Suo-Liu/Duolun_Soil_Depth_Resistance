setwd("C:\\Users\\True\\OneDrive\\桌面")
# global network calculation ####
require(MENA)

com.file <- "protist_zotu.txt"
prefix.m <- "Protist"
prefix1 <- "Pro"

comm <- t(read.table(com.file,
                     header = TRUE, sep = "\t", row.names = 1,
                     as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE
))

cor.meth <- "pearson"
maj.1 <- 0.25

prefixm <- cor.meth

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
## find cutoff ####
for (j in c("L1", "L2", "L3", "L4")) {
  used.comm <- comm
  used.treat <- subset(treat, Layer == j)
  
  used.comm = used.comm[match(rownames(used.treat), rownames(used.comm)),]
  
  dim(used.comm)
  used.comm <- used.comm[, which((colSums(used.comm > 0) / dim(used.comm)[1]) >= maj.1)]
  dim(used.comm)
  
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
  save(output1, file = paste0(j, ".", prefixm, ".", prefix1, ".m", maj.1, ".CLR.FillPairLB.rda"))
}
## after find cutoff ####
for (j in c("L1", "L2", "L3", "L4")) {
  used.comm <- comm
  used.treat <- subset(treat, Layer == j)
  used.comm = used.comm[match(rownames(used.treat), rownames(used.comm)),]
  
  dim(used.comm)
  used.comm <- used.comm[, which((colSums(used.comm > 0) / dim(used.comm)[1]) >= maj.1)]
  dim(used.comm)
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
  write.csv(assmc1, file = paste0(j, ".", prefixm, ".", prefix1, ".m", maj.1, ".CLR.FillPairLB.csv"))
}

# small networks: one sample as one network ####
library(ieggr)

source("newnetindex.R")
treat.file <- "treatment.csv"

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

pre.treat <- read.csv(treat.file, header = TRUE, sep = ",", row.names = 1)
for (j in c("L1", "L2", "L3", "L4")) {
  comm <- pre.comm
  treat <- subset(pre.treat, Layer == j)
  
  comm = comm[match(rownames(treat), rownames(comm)),]
  
  for (k in rownames(treat)) {
    dat <- read.csv(paste0(j, ".", prefixm, ".", 
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
