# Functions ####
# install.packages("OmicsPLS")
# library("OmicsPLS")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ropls")
library(ropls)
R2ff<-function(xm,ym,o2p)
{
  c(R2x=mean(sapply(1:ncol(xm),function(i){1-sum((o2p$X_hat[,i]-xm[,i])^2)/sum((xm[,i]-mean(xm[,i]))^2)})),
    R2y=mean(sapply(1:ncol(ym),function(i){1-sum((o2p$Y_hat[,i]-ym[,i])^2)/sum((ym[,i]-mean(ym[,i]))^2)})))
}
R2sep<-function(xm,pls){(((getVipVn(pls))^2)/ncol(xm))*getSummaryDF(pls)[,'R2Y(cum)']}
plstest<-function(xm,ym,rand=100)
{
  nx=ncol(xm)
  combs=list()
  for(i in 1:nx)
  {
    message("i=",i," ",date())
    cbni=combn(nx,i)
    combs=c(combs,lapply(1:ncol(cbni),function(i){cbni[,i]}))
  }
  message("Total of ",length(combs)," combinations. ",date())
  #trac=seq(from=1,to=length(combs),by=50)
  id=t(rep(0,nx))
  colnames(id)=colnames(xm)
  
  dfna=t(rep(NA,8))
  colnames(dfna)=c('R2X(cum)','R2Y(cum)','Q2(cum)','RMSEE','pre','ort','pR2Y','pQ2')
  
  res=lapply(1:length(combs),
             function(i)
             {
               message("----- PLS i=",i," ",date())
               xmi=xm[,combs[[i]],drop=FALSE]
               pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
               if(class(pls)=="try-error")
               {
                 pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
               }
               out=dfna
               if(class(pls)!="try-error"){dfi=getSummaryDF(pls);out[1,match(colnames(dfi),colnames(dfna))]=as.vector(as.matrix(dfi[1,]))}
               idni=id
               idni[combs[[i]]]=1
               cbind(idni,out)
             })
  res
}
plsfw<-function(xm,ym,r2buf=0.98,Q2ck=FALSE,SEEck=FALSE,rand=100)
{
  nx=ncol(xm)
  fs<-list(integer(0))
  dfna=t(rep(NA,8))
  colnames(dfna)=c('R2X(cum)','R2Y(cum)','Q2(cum)','RMSEE','pre','ort','pR2Y','pQ2')
  id=t(rep(0,nx))
  colnames(id)=colnames(xm)
  R2Yn=list()
  fijn=list()
  outrc=list()
  if(Q2ck){Q2n=list()}
  if(SEEck){SEEn=list()}
  k=1
  kn=1
  for(fn in 1:nx)
  {
    for(i in 1:length(fs))
    {
      fsi=fs[[i]]
      fai=which(!((1:nx) %in% fsi))
      for(j in 1:length(fai))
      {
        fij=c(fsi,fai[j])
        
        message("----- PLS fn=",fn," in ",nx,", i=",i," in ",length(fs),", j=",j," in ",length(fai),". ",date())
        xmi=xm[,fij,drop=FALSE]
        pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
        if(class(pls)=="try-error")
        {
          pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=rand,fig.pdfC='none',info.txtC='none'))
        }
        out=dfna
        if(class(pls)!="try-error"){dfi=getSummaryDF(pls);out[1,match(colnames(dfi),colnames(dfna))]=as.vector(as.matrix(dfi[1,]))}
        idni=id
        idni[fij]=1
        outrc[[k]]=cbind(idni,out)
        k=k+1
        R2Yn[[kn]]=out[,"R2Y(cum)"][[1]]
        if(Q2ck){Q2n[[kn]]=out[,"Q2(cum)"][[1]]}
        if(SEEck){SEEn[[kn]]=out[,"RMSEE"][[1]]}
        fijn[[kn]]=fij
        kn=kn+1
      }
    }
    maxR2n=max(unlist(R2Yn),na.rm = TRUE)
    kns=which(unlist(R2Yn)>=(maxR2n*r2buf))
    if(Q2ck)
    {
      if(length(kns)>1)
      {
        Q2nk=Q2n[kns]
        maxQ2nk=max(unlist(Q2nk),na.rm = TRUE)
        if(maxQ2nk>=0){maxQ2nkb=maxQ2nk*r2buf}else{maxQ2nkb=maxQ2nk*(1+(1-r2buf))}
        knsk1=which(unlist(Q2nk)>=maxQ2nkb)
        kns=kns[knsk1]
      }
    }
    if(SEEck)
    {
      if(length(kns)>1)
      {
        SEEnk=SEEn[kns]
        minSEE=min(unlist(SEEnk),na.rm = TRUE)
        knsk2=which(unlist(SEEnk)<=(minSEE*(1+(1-r2buf))))
        kns=kns[knsk2]
      }
    }
    EPS <- (.Machine$double.eps)
    if((max(kns)<=length(fs)) | sum(is.na(kns))>0 | maxR2n >= (1-EPS)){break}else{
      fs=fijn[kns[which(kns>length(fs))]]
      R2Yn=R2Yn[kns]
      fijn=fijn[kns]
      kn=length(R2Yn)+1
    }
  }
  outrcm=Reduce(rbind,outrc)
}
rp.pls<-function(xm,ym,rand=100)
{
  R2sepf<-function(xm,ym,predI)
  {
    pls=try(opls(x=xm,y=ym,predI=predI,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))
    if(class(pls)=="try-error"){pls=try(opls(x=xm,y=ym,predI=1,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))}
    if(class(pls)=="try-error"){out1=rep(NA,1+ncol(xm))}else{
      out1=c(R2Y=getSummaryDF(pls)[,"R2Y(cum)"],(((getVipVn(pls))^2)/ncol(xm))*getSummaryDF(pls)[,'R2Y(cum)'])
    }
    out1
  }
  pls=try(opls(x=xm,y=ym,predI=NA,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))
  predI=NA
  if(class(pls)=="try-error"){predI=1;pls=try(opls(x=xm,y=ym,predI=1,orthoI=0,permI=100,info.txtC="none",fig.pdfC="none"))}
  if(class(pls)=="try-error"){out=rep(NA,2*(1+ncol(xm)))}else{
    R2obs=R2sepf(xm,ym,predI)
    perm=permute::shuffleSet(nrow(xm),nset = rand)
    tracs=seq(from=1,to=nrow(perm),by=20)
    R2rm=sapply(1:nrow(perm),
                function(i)
                {
                  if(i %in% tracs){message("i=",i,". ",date())}
                  ymri=ym[perm[i,],,drop=FALSE]
                  rownames(ymri)=rownames(ym)
                  R2sepf(xm,ymri,predI)
                })
    EPS <- (.Machine$double.eps)
    dR2=((R2rm-R2obs)>=(-EPS))
    out=c(R2obs,rowSums(dR2,na.rm = TRUE)/rand)
  }
  names(out)=c(paste0("R2.",c('Y',colnames(xm))),paste0("P.",c('R2Y',colnames(xm))))
  out
}
r2adj<-function(r2,n,p)
{
  idx=which((n-p-1)<0)
  out=1-((1-r2)*((n-1)/(n-p-1)))
  out[idx]=NA
  out
}
# PLS construction ####
setwd("C:\\Users\\True\\OneDrive\\桌面")
library("ieggr")
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"
env.data.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\env_used.csv"

treat <- read.csv(treat.file, header = T, row.names = 1)
env.data <- read.csv(env.data.file, header = T, row.names = 1)
treat <- subset(treat, plant.type == "TS")

# divertsity index ####
Bac.divindex.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\alpha\\diversity index\\alpha index_", "Bacteria", ".csv")
Bac.divindex <- read.csv(Bac.divindex.file, row.names = 1, header = T, sep = ",")

Bac.divindex <- Bac.divindex[rownames(Bac.divindex) %in% rownames(treat), ]

Bac.index.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Bac.divindex = Bac.divindex),silent = T)
  used.index <- sampc$Bac.divindex
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  Bac.index.table <- rbind(Bac.index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }
  
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Bac.divindex = Bac.divindex), silent = T)
  used.index <- sampc$Bac.divindex
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", used.index$treat)
  Bac.index.table <- rbind(Bac.index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Bac.divindex = Bac.divindex), silent = T)
  used.index <- sampc$Bac.divindex
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", "W")
  Bac.index.table <- rbind(Bac.index.table, used.index)
}
Bac.index.table <- Bac.index.table[order(Bac.index.table$treat), ]
Bac.index.table <- Bac.index.table[order(Bac.index.table$block), ]
Bac.index.table <- Bac.index.table[order(Bac.index.table$Layer), ]

Pro.divindex.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\alpha\\diversity index\\alpha index_", "Protist", ".csv")
Pro.divindex <- read.csv(Pro.divindex.file, row.names = 1, header = T, sep = ",")

Pro.divindex <- Pro.divindex[rownames(Pro.divindex) %in% rownames(treat), ]

Pro.index.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Pro.divindex = Pro.divindex), , silent = T)
  used.index <- sampc$Pro.divindex
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  Pro.index.table <- rbind(Pro.index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }
  
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Pro.divindex = Pro.divindex),silent = T)
  used.index <- sampc$Pro.divindex
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", used.index$treat)
  Pro.index.table <- rbind(Pro.index.table, used.index)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Pro.divindex = Pro.divindex), silent = T)
  used.index <- sampc$Pro.divindex
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", "W")
  Pro.index.table <- rbind(Pro.index.table, used.index)
}
Pro.index.table <- Pro.index.table[order(Pro.index.table$treat), ]
Pro.index.table <- Pro.index.table[order(Pro.index.table$block), ]
Pro.index.table <- Pro.index.table[order(Pro.index.table$Layer), ]

# resistance ####
Bac.result.data.file <- paste0(
  "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\", "beta",
  "\\resistance\\", "beta", "_resistance_", "Bacteria", ".csv"
)
Bac.result.data <- read.csv(Bac.result.data.file, header = T, row.names = 1)
Bac.result.data <- subset(Bac.result.data, treat == "Sorenson" &
                            plant1 == "TS")

Pro.result.data.file <- paste0(
  "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\", "beta",
  "\\resistance\\", "beta", "_resistance_", "Protist", ".csv"
)
Pro.result.data <- read.csv(Pro.result.data.file, header = T, row.names = 1)
Pro.result.data <- subset(Pro.result.data, treat == "Sorenson" &
                            plant1 == "TS")
# Bacterial trait ####
Bac.com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"

Bac.com <- t(read.table(Bac.com.file,
                        header = TRUE, sep = "\t", row.names = 1,
                        as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                        check.names = FALSE
))
name.row <- rownames(Bac.com)
name.row[which(rownames(Bac.com) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(Bac.com) == "YD52L3")] <- "YD3L3"
rownames(Bac.com) <- name.row

library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, Bac.com = Bac.com), silent = T)
dim(Bac.com)
Bac.com <- sampc$Bac.com
treat <- sampc$treat
Bac.com <- Bac.com[, colSums(Bac.com) > 0]
dim(Bac.com)

Bac.trait.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\predicted_traits.csv"

Bac.trait <- read.csv(Bac.trait.file, header = T, row.names = 1)
Bac.trait$cell_diameter <- 10^Bac.trait$cell_diameter
Bac.trait$cell_length <- 10^Bac.trait$cell_length
Bac.trait$doubling_h <- 10^Bac.trait$doubling_h
Bac.trait <- Bac.trait[, 1:11]

Bac.otu.table <- t(Bac.com)

spc <- match.name(rn.list = list(Bac.otu.table = Bac.otu.table, Bac.trait = Bac.trait), silent = T)
Bac.otu.table <- spc$Bac.otu.table
Bac.trait <- spc$Bac.trait

data.list <- list()
for (i in 1:ncol(Bac.trait)) {
  for (j in 1:ncol(Bac.otu.table)) {
    Bac.otu.table[, j] <- Bac.otu.table[, j] * Bac.trait[, i]
  }
  data.list[[i]] <- Bac.otu.table
  names(data.list)[i] <- colnames(Bac.trait)[i]
}
trait_names <- names(Bac.trait)
used.adjusted_mats <- list()
for (trait in trait_names) {
  adjusted_mat <- Bac.otu.table / Bac.trait[[trait]]
  adjusted_mat[is.infinite(adjusted_mat)] <- NA
  
  used.adjusted_mats[[trait]] <- t(adjusted_mat)
}
Bac.weighted.result <- data.frame(row.names = colnames(Bac.otu.table))
for (trait in trait_names) {
  w.sample.aver <- colSums(Bac.otu.table) / colSums(Bac.otu.table / Bac.trait[[trait]], na.rm = TRUE)
  Bac.weighted.result[[trait]] <- w.sample.aver
}
Bac.unweighted.result <- data.frame(row.names = colnames(Bac.otu.table))
for (trait in trait_names) {
  uw.sample.aver <- colSums(Bac.otu.table > 0) / colSums((Bac.otu.table > 0) / Bac.trait[[trait]], na.rm = TRUE)
  Bac.unweighted.result[[trait]] <- uw.sample.aver
}

Bac.weighted.result$name <- rownames(Bac.weighted.result)
Bac.weighted.result$method <- "weighted"

Bac.weighted.result$name <- rownames(Bac.unweighted.result)
Bac.unweighted.result$method <- "unweighted"

Bac.trait.table <- data.frame(stringsAsFactors = F)
for (i in c("RP", "EP", "W", "WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "C")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Bac.weighted.result = Bac.weighted.result), , silent = T)
  used.index <- sampc$Bac.weighted.result
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  Bac.trait.table <- rbind(Bac.trait.table, used.index)
}
for (i in c("WEP", "WRP")) {
  if (i == "WEP") {
    ct.treat <- subset(treat, combined_treat1 == "EP")
  } else {
    ct.treat <- subset(treat, combined_treat1 == "RP")
  }
  
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Bac.weighted.result = Bac.weighted.result), , silent = T)
  used.index <- sampc$Bac.weighted.result
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", used.index$treat)
  Bac.trait.table <- rbind(Bac.trait.table, used.index)
}
for (i in c("WEP", "WRP")) {
  ct.treat <- subset(treat, combined_treat1 == "W")
  ct.treat <- ct.treat[order(ct.treat$block), ]
  ct.treat <- ct.treat[order(ct.treat$Layer), ]
  
  sampc <- match.name(rn.list = list(ct.treat = ct.treat, Bac.weighted.result = Bac.weighted.result), silent = T)
  used.index <- sampc$Bac.weighted.result
  ct.treat <- sampc$ct.treat
  
  used.index <- as.data.frame(used.index)
  used.index$block <- ct.treat$block
  used.index$Layer <- ct.treat$Layer
  used.index$treat <- ct.treat$combined_treat1
  used.index$treat <- paste0(i, ".", "W")
  Bac.trait.table <- rbind(Bac.trait.table, used.index)
}
Bac.trait.table <- Bac.trait.table[order(Bac.trait.table$treat), ]
Bac.trait.table <- Bac.trait.table[order(Bac.trait.table$block), ]
Bac.trait.table <- Bac.trait.table[order(Bac.trait.table$Layer), ]

# Protistan trait ####
Pro.com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"

Pro.com <- t(read.table(Pro.com.file,
                        header = TRUE, sep = "\t", row.names = 1,
                        as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                        check.names = FALSE
))
name.row <- rownames(Pro.com)
name.row[which(rownames(Pro.com) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(Pro.com) == "YD52L3")] <- "YD3L3"
rownames(Pro.com) <- name.row

library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, Pro.com = Pro.com), silent = T)
dim(Pro.com)
Pro.com <- sampc$Pro.com
treat <- sampc$treat
Pro.com <- Pro.com[, colSums(Pro.com) > 0]
dim(Pro.com)

clas.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\18S_tax.txt"
clas <- read.table(clas.file,
                   header = TRUE, sep = "\t", row.names = 1,
                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                   check.names = FALSE
)
## shell ####
shell.trait <- clas[, c("shell"), drop = F]
shell.trait <- shell.trait[shell.trait$shell != "", , drop = F]
shell.trait <- shell.trait[shell.trait$shell != "naked_or_silica_or_organic", , drop = F]
shell.trait$shell <- ifelse(shell.trait$shell == c("naked"), 0, 1)

sampc <- match.name(cn.list = list(Pro.com = Pro.com), rn.list = list(shell.trait = shell.trait), silent = T)
shell.comm <- sampc$Pro.com
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

## locomotion ####
locomotion.trait <- clas[, c("locomotion"), drop = F]
locomotion.trait <- locomotion.trait[locomotion.trait$locomotion != "", , drop = F]

locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion == c("non_motile"), 0, locomotion.trait$locomotion)
locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion == c("pseudopodia_and_flagella"), 2, locomotion.trait$locomotion)
locomotion.trait$locomotion <- ifelse(locomotion.trait$locomotion %in% c(
  "cilia", "flagella", "pseudopodia",
  "pseudopodia_and_flagella_or_pseudopodia_or_flagella"
), 1, locomotion.trait$locomotion)
locomotion.trait$locomotion <- as.numeric(locomotion.trait$locomotion)
sampc <- match.name(cn.list = list(Pro.com = Pro.com), rn.list = list(locomotion.trait = locomotion.trait), silent = T)
locomotion.comm <- sampc$Pro.com
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

## bacterivore ####
bacterivore.trait <- clas[, c("main_functional_class", "associated_organism"), drop = F]
bacterivore.trait <- bacterivore.trait[bacterivore.trait$main_functional_class %in% c("predator", "predator (add)"), , drop = F]
bacterivore.trait$associated_organism <- ifelse(grepl("bacteria", bacterivore.trait$associated_organism), 1, 0)
bacterivore.trait <- bacterivore.trait[, -1, drop = F]

sampc <- match.name(cn.list = list(Pro.com = Pro.com), rn.list = list(bacterivore.trait = bacterivore.trait), silent = T)
bacterivore.comm <- sampc$Pro.com
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

## body size ####
size.trait <- clas[, c("body.size"), drop = F]
size.trait <- size.trait[!is.na(size.trait$body.size), , drop = F]

sampc <- match.name(cn.list = list(Pro.com = Pro.com), rn.list = list(size.trait = size.trait), silent = T)
size.comm <- sampc$Pro.com
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

Pro.trait.table <- data.frame(stringsAsFactors = F)
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
  Pro.trait.table <- rbind(Pro.trait.table, used.trait)
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
  Pro.trait.table <- rbind(Pro.trait.table, used.trait)
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
  Pro.trait.table <- rbind(Pro.trait.table, used.trait)
}
Pro.trait.table <- Pro.trait.table[order(Pro.trait.table$treat), ]
Pro.trait.table <- Pro.trait.table[order(Pro.trait.table$block), ]
Pro.trait.table <- Pro.trait.table[order(Pro.trait.table$Layer), ]

# env ####
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

sem.data = cbind(env.table[,c(1:4,11,19,23,26,27,28,33)],Bac.index.table[,1],
                 Bac.result.data[,1],Bac.trait.table[,1:11],Pro.index.table[,1],
                 Pro.result.data[,1], Pro.trait.table[,1:4])
colnames(sem.data)[c(12:13,25:26)] = c("Bac.richness","Bac.resistance","Pro.richness","Pro.resistance")
colnames(sem.data)
library(dplyr)
sem.data = scale(sem.data)
sem.data = as.data.frame(sem.data)
sem.data = sem.data[,c(1:13, 14:26)]

xms1 = sem.data
# 定义所需要解释的变量
# 1 ####
be=xms1[,'HR',drop=FALSE]
ym=as.matrix(be)
Yname=paste(colnames(ym),collapse = "_")

colnames(xms1)

# 选择一些数据用于构建PLS
xms2=as.matrix(xms1[,c("NO3", 
                       "NH4",
                       "moisture.one.nighbor.aver",
                       "root.biomass",
                       "BNPP",
                       "Bac.resistance",
                       "Pro.resistance",
                       "Bac.richness",
                       "Pro.richness"
), drop=FALSE])

head(xms2)

plstm=plsfw(xms2,ym,Q2ck=T,SEEck=T)
ieggr::save.file(plstm,prefix="PLS",filename = paste0(Yname,".PLSFWtest"))

# Optimize result: R2Y larger than 98% of maximum R2Y, and P values for R2Y and Q2 should be <0.05; 
# the least factor number; if multiple hits, use minimum RMSEE; if still multiple hits, choose the one with better biological sense.
sel.id = c(2,3,5,6,7,8,9)
xmi=xms2[,sel.id,drop=FALSE]
head(xmi)

pls=try(opls(x=xmi,y=ym,predI=NA,orthoI=0,permI=1000,fig.pdfC="none"))
if(class(pls)=="try-error"){
  pls=try(opls(x=xmi,y=ym,predI=1,orthoI=0,permI=1000,fig.pdfC="none"))
}

vip=getVipVn(pls)
names(vip)=paste0("VIP.",names(vip))
rpi=R2sep(xmi,pls)
names(rpi)=paste0("R2.",names(rpi))

R2sig<-Psig<-rep(NA,ncol(xmi))
names(R2sig)<-paste0("R2.Single.",colnames(xmi))
names(Psig)<-paste0("P.Single.",colnames(xmi))

for(j in 1:ncol(xmi))
{
  message("-----Single j=",j,". ",date())
  plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=NA,orthoI=0,
                permI=1000,fig.pdfC="none"))
  # if(class(plsj)=="try-error"){plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=1,
  #                                            orthoI=0,permI=1000,fig.pdfC="none"))}
  if(class(plsj)=="try-error"|sum(dim(plsj@modelDF)) ==0){plsj=try(opls(x=xmi[,j,drop=FALSE],y=ym,predI=1,
                                                                        orthoI=0,permI=1000,fig.pdfC="none"))}
  sdfj=getSummaryDF(plsj)
  R2sig[j]=sdfj[,'R2Y(cum)'][[1]]
  if('pR2Y' %in% colnames(sdfj)){Psig[j]=sdfj[,'pR2Y'][[1]]}
  if(is.na(Psig[j])){rpj=rp.pls(xm=xmi[,j,drop=FALSE],ym=ym,rand = 100);Psig[j]=rpj['P.R2Y']}
}

save.file(t(cbind(t(c(Y=Yname,R2sig,Psig,rpi,vip)),(getSummaryDF(pls)))),
          prefix="PLS",filename = paste0(Yname,".modelselectedFW"))
