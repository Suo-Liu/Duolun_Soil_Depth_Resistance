setwd("C:\\Users\\True\\OneDrive\\桌面")

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"
library(ieggr)

bac.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\alpha\\diversity index\\alpha index_Bacteria.csv"
fungi.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\alpha\\diversity index\\alpha index_Fungi.csv"
pro.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\alpha\\diversity index\\alpha index_Protist.csv"

bac.divindex.data <- read.csv(bac.divindex.file, row.names = 1, header = T, sep = ",")
fungi.divindex.data <- read.csv(fungi.divindex.file, row.names = 1, header = T, sep = ",")
pro.divindex.data <- read.csv(pro.divindex.file, row.names = 1, header = T, sep = ",")

bac.divindex.data <- bac.divindex.data[, 1, drop = F]
fungi.divindex.data <- fungi.divindex.data[, 1, drop = F]
pro.divindex.data <- pro.divindex.data[, 1, drop = F]

divindex.data <- cbind(bac.divindex.data, fungi.divindex.data, pro.divindex.data)
colnames(divindex.data) <- c("Bac_richness", "Fungi_richness", "Protist_richness")

treat <- read.csv(treat.file, header = T, row.names = 1)
treat <- subset(treat, plant.type == "TS")

library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, divindex.data = divindex.data))
divindex.data <- sampc$divindex.data
treat <- sampc$treat

## plot ####
source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\rsquaredglmm.r")
tdcm_index<-function(betai,treat,prefixi = NULL,scale.num = F,rand=1000)
{
  library(lme4)
  library(car)
  betai = as.data.frame(betai)
  treat = as.data.frame(treat)
  
  tdc.lmm<-function(betai,treat,save.output=FALSE,scale.num = F,prefixi=NULL)
  {
    Layer=treat$Layer
    Layer.lev=unique(Layer)
    out=list()
    for(j in 1:length(Layer.lev))
    {
      idj=which(Layer==Layer.lev[j])
      betaij <- betai[, 1][idj]
      xij <- betai[, 2][idj]
      blockij <- treat$block[idj]
      treatij <- treat$combined_treat1[idj]
      if (scale.num){
        betaij = scale(betaij)
        xij = scale(xij)
      }
      lmij <- lmer(betaij ~ xij + (1 | blockij) + (1 | treatij))
      lmijsm=summary(lmij)
      AIC1=AIC(lmij)
      r2ij=rsquared.glmm(lmij)
      lmijCS=Anova(lmij,type = "II")
      if(save.output)
      {
        dataij=data.frame(layer = Layer.lev[j],
                          richness_1 = betaij,
                          richness_2 = xij,
                          block = blockij,
                          treat = treatij)
        save.file(dataij,prefix = prefixi,filename = paste0("index.Data.",Layer.lev[j]))
        sink(file=paste0(prefixi,".index.LMM.",Layer.lev[j],".txt"))
        print("------ LMM model result -------")
        print(lmij)
        print("------ LMM model R squared -------")
        print(r2ij)
        print("------ LMM model summary -------")
        print(lmijsm)
        print("------ LMM model Wald Type II chisquare test -------")
        print(lmijCS)
        sink()
      }
      out[[j]]=c(slope.fix=lmijsm$coefficients[2,1],slope.se=lmijsm$coefficients[2, 2],
                 R2M=r2ij$Marginal,R2C=r2ij$Conditional,AIC1=AIC1,AIC2=r2ij$AIC,
                 P.typeII=lmijCS[[3]],Chisq=lmijCS[[1]])
    }
    outs=Reduce(cbind,out)
    colnames(outs)=Layer.lev
    outs
  }
  
  tdci=tdc.lmm(betai = betai,treat = treat,
               scale.num = scale.num,
               save.output = FALSE,
               prefixi = prefixi)
  r2.obs=as.vector(tdci[3:4,])
  aic.obs=as.vector(tdci[5:6,])
  ds.obs=(-tdci[1,1])-(-tdci[1,2])
  
  # randomize time points and other the same as observed.
  layer.lev = unique(treat$Layer)
  layer.perm = vegan:::getPermuteMatrix(rand,nrow(betai))
  trace.seq = seq(from = 1,to = rand,by = 100)
  ind.rand = lapply(1:nrow(layer.perm),
                    function(k)
                    {
                      if(k %in% trace.seq) message("-------Now randomizing k=",k,". ",date())
                      out=list()
                      idi=layer.perm[k,]
                      perm.treat=treat
                      perm.treat[,"Layer"][idi] = c(rep(layer.lev[1], nrow(betai)/2),
                                                    rep(layer.lev[2], nrow(betai)/2))
                      tdcr=tdc.lmm(betai = betai, treat = perm.treat,scale.num = scale.num)
                      out$r2=as.vector(tdcr[3:4,])
                      out$aic=as.vector(tdcr[5:6,])
                      out$ds=(-tdcr[1,1])-(-tdcr[1,2])
                      out
                    })
  r2.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$r2})
  aic.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$aic})
  ds.ran=sapply(1:length(ind.rand),function(k){ind.rand[[k]]$ds})
  EPS <- sqrt(.Machine$double.eps)
  p.r2=(rowSums(r2.ran>=(matrix(r2.obs,nr=nrow(r2.ran),nc=ncol(r2.ran))-EPS))+1)/(ncol(r2.ran)+1)
  p.aic=(rowSums(aic.ran<=(matrix(aic.obs,nr=nrow(aic.ran),nc=ncol(aic.ran))+EPS))+1)/(ncol(aic.ran)+1)
  
  r2.m.diff = r2.obs[1]-r2.obs[3]
  r2.c.diff = r2.obs[2]-r2.obs[4]
  
  p.m.r2.diff = sum((r2.ran[1,] - r2.ran[3,]) <= r2.m.diff)/ncol(r2.ran)
  p.c.r2.diff = sum((r2.ran[2,] - r2.ran[4,]) <= r2.c.diff)/ncol(r2.ran)
  
  if(ds.obs>0)
  {
    p.perm=(sum(ds.ran>=(ds.obs-EPS))+1)/(length(ds.ran)+1)
  }else{
    p.perm=(sum(ds.ran<=(ds.obs+EPS))+1)/(length(ds.ran)+1)
  }
  p.values=rbind(matrix(p.r2,2,2),matrix(p.aic,2,2),c(p.perm,NA),
                 c(p.m.r2.diff,NA),c(p.c.r2.diff,NA))
  rownames(p.values)=c("P.R2M","P.R2C","P.AIC1","P.AIC2","P.dS.perm",
                       "P.R2M.diff","P.R2C.diff")
  output=rbind(tdci,p.values)
  output = as.data.frame(output)
  output
}
betai = divindex.data[,c(1,3)]
treat$Layer = ifelse(treat$Layer %in% c("L1", "L2"), "Topsoil", "Subsoil")
results.lmm <- tdcm_index(betai = betai, treat = treat, 
                          scale.num = F, rand=1000)
write.csv(results.lmm, file = paste0("richness_LMM_perm test", ".csv"))
