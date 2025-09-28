setwd("C:\\Users\\True\\OneDrive\\桌面")

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"
library(ieggr)

bac.com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
fungi.com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\ITS\\all samples\\Unoise\\zotutab220_max_2_resample_10000.txt"
pro.com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"

bac.comm <- t(read.table(bac.com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))
name.row <- rownames(bac.comm)
name.row[which(rownames(bac.comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(bac.comm) == "YD52L3")] <- "YD3L3"
rownames(bac.comm) <- name.row

fungi.comm <- t(read.table(fungi.com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))
name.row <- rownames(fungi.comm)
name.row[which(rownames(fungi.comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(fungi.comm) == "YD52L3")] <- "YD3L3"
rownames(fungi.comm) <- name.row

pro.comm <- t(read.table(pro.com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))
name.row <- rownames(pro.comm)
name.row[which(rownames(pro.comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pro.comm) == "YD52L3")] <- "YD3L3"
rownames(pro.comm) <- name.row

treat <- read.csv(treat.file, header = T, row.names = 1)
treat <- subset(treat, plant.type == "TS")

library(ieggr)
sampc <- match.name(
  rn.list = list(
    treat = treat, bac.comm = bac.comm, fungi.comm = fungi.comm,
    pro.comm = pro.comm
  ),
  silent = TRUE
)
bac.comm <- sampc$bac.comm
fungi.comm <- sampc$fungi.comm
pro.comm <- sampc$pro.comm
treat <- sampc$treat

bac.comm <- bac.comm[, colSums(bac.comm) > 0]
fungi.comm <- fungi.comm[, colSums(fungi.comm) > 0]
pro.comm <- pro.comm[, colSums(pro.comm) > 0]

library(vegan)
dist.used <- vegdist(bac.comm, method = "bray", binary = T)
Bac.table <- dist.3col(dist.used)

dist.used <- vegdist(fungi.comm, method = "bray", binary = T)
Fungi.table <- dist.3col(dist.used)

dist.used <- vegdist(pro.comm, method = "bray", binary = T)
Pro.table <- dist.3col(dist.used)

Bac.table$plot1 <- treat$plot[match(Bac.table$name1, rownames(treat))]
Bac.table$plot2 <- treat$plot[match(Bac.table$name2, rownames(treat))]
Bac.table$layer1 <- treat$Layer[match(Bac.table$name1, rownames(treat))]
Bac.table$layer2 <- treat$Layer[match(Bac.table$name2, rownames(treat))]
Bac.table$treat1 <- treat$combined_treat1[match(Bac.table$name1, rownames(treat))]
Bac.table$treat2 <- treat$combined_treat1[match(Bac.table$name2, rownames(treat))]
Bac.table$block1 <- treat$block[match(Bac.table$name1, rownames(treat))]
Bac.table$block2 <- treat$block[match(Bac.table$name2, rownames(treat))]

Fungi.table$plot1 <- treat$plot[match(Fungi.table$name1, rownames(treat))]
Fungi.table$plot2 <- treat$plot[match(Fungi.table$name2, rownames(treat))]
Fungi.table$layer1 <- treat$Layer[match(Fungi.table$name1, rownames(treat))]
Fungi.table$layer2 <- treat$Layer[match(Fungi.table$name2, rownames(treat))]
Fungi.table$treat1 <- treat$combined_treat1[match(Fungi.table$name1, rownames(treat))]
Fungi.table$treat2 <- treat$combined_treat1[match(Fungi.table$name2, rownames(treat))]
Fungi.table$block1 <- treat$block[match(Fungi.table$name1, rownames(treat))]
Fungi.table$block2 <- treat$block[match(Fungi.table$name2, rownames(treat))]

Pro.table$plot1 <- treat$plot[match(Pro.table$name1, rownames(treat))]
Pro.table$plot2 <- treat$plot[match(Pro.table$name2, rownames(treat))]
Pro.table$layer1 <- treat$Layer[match(Pro.table$name1, rownames(treat))]
Pro.table$layer2 <- treat$Layer[match(Pro.table$name2, rownames(treat))]
Pro.table$treat1 <- treat$combined_treat1[match(Pro.table$name1, rownames(treat))]
Pro.table$treat2 <- treat$combined_treat1[match(Pro.table$name2, rownames(treat))]
Pro.table$block1 <- treat$block[match(Pro.table$name1, rownames(treat))]
Pro.table$block2 <- treat$block[match(Pro.table$name2, rownames(treat))]

Bac.table <- Bac.table[Bac.table$treat1 == Bac.table$treat2 &
  Bac.table$layer1 == Bac.table$layer2, ]
Fungi.table <- Fungi.table[Fungi.table$treat1 == Fungi.table$treat2 &
  Fungi.table$layer1 == Fungi.table$layer2, ]
Pro.table <- Pro.table[Pro.table$treat1 == Pro.table$treat2 &
  Pro.table$layer1 == Pro.table$layer2, ]

Dist.table <- cbind(
  Bac.table$dis, Fungi.table$dis, Pro.table$dis,
  Bac.table$plot1, Bac.table$layer1, Bac.table$treat1,
  Bac.table$block1, Bac.table$block2
)
Dist.table <- as.data.frame(Dist.table)
colnames(Dist.table) <- c("Bac", "Fungi", "Pro", "plot", "Layer", "treat", "block1", "block2")
Dist.table$Bac <- as.numeric(Dist.table$Bac)
Dist.table$Fungi <- as.numeric(Dist.table$Fungi)
Dist.table$Pro <- as.numeric(Dist.table$Pro)
Dist.table$block <- paste0(Dist.table$block1, "_", Dist.table$block2)

## plot ####
source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\rsquaredglmm.r")
tdcm_index <- function(betai, treat, prefixi = NULL, scale.num = F, rand = 1000) {
  library(lme4)
  library(car)
  betai <- as.data.frame(betai)
  treat <- as.data.frame(treat)

  tdc.lmm <- function(betai, treat, save.output = FALSE, scale.num = F, prefixi = NULL) {
    Layer <- treat$Layer
    Layer.lev <- unique(Layer)
    out <- list()
    for (j in 1:length(Layer.lev))
    {
      idj <- which(Layer == Layer.lev[j])
      betaij <- betai[, 1][idj]
      xij <- betai[, 2][idj]
      blockij <- treat$block[idj]
      treatij <- treat$treat[idj]
      if (scale.num) {
        betaij <- scale(betaij)
        xij <- scale(xij)
      }
      lmij <- lmer(betaij ~ xij + (1 | blockij) + (1 | treatij))
      lmijsm <- summary(lmij)
      AIC1 <- AIC(lmij)
      r2ij <- rsquared.glmm(lmij)
      lmijCS <- Anova(lmij, type = "II")
      if (save.output) {
        dataij <- data.frame(
          layer = Layer.lev[j],
          richness_1 = betaij,
          richness_2 = xij,
          block = blockij,
          treat = treatij
        )
        save.file(dataij, prefix = prefixi, filename = paste0("index.Data.", Layer.lev[j]))
        sink(file = paste0(prefixi, ".index.LMM.", Layer.lev[j], ".txt"))
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
      out[[j]] <- c(
        slope.fix = lmijsm$coefficients[2, 1], slope.se = lmijsm$coefficients[2, 2],
        R2M = r2ij$Marginal, R2C = r2ij$Conditional, AIC1 = AIC1, AIC2 = r2ij$AIC,
        P.typeII = lmijCS[[3]], Chisq = lmijCS[[1]]
      )
    }
    outs <- Reduce(cbind, out)
    colnames(outs) <- Layer.lev
    outs
  }

  tdci <- tdc.lmm(
    betai = betai, treat = treat,
    scale.num = scale.num,
    save.output = FALSE,
    prefixi = prefixi
  )
  r2.obs <- as.vector(tdci[3:4, ])
  aic.obs <- as.vector(tdci[5:6, ])
  ds.obs <- (-tdci[1, 1]) - (-tdci[1, 2])

  # randomize time points and other the same as observed.
  layer.lev <- unique(treat$Layer)
  layer.perm <- vegan:::getPermuteMatrix(rand, nrow(betai))
  trace.seq <- seq(from = 1, to = rand, by = 100)
  ind.rand <- lapply(
    1:nrow(layer.perm),
    function(k) {
      if (k %in% trace.seq) message("-------Now randomizing k=", k, ". ", date())
      out <- list()
      idi <- layer.perm[k, ]
      perm.treat <- treat
      perm.treat[, "Layer"][idi] <- c(
        rep(layer.lev[1], nrow(betai) / 2),
        rep(layer.lev[2], nrow(betai) / 2)
      )
      tdcr <- tdc.lmm(betai = betai, treat = perm.treat, scale.num = scale.num)
      out$r2 <- as.vector(tdcr[3:4, ])
      out$aic <- as.vector(tdcr[5:6, ])
      out$ds <- (-tdcr[1, 1]) - (-tdcr[1, 2])
      out
    }
  )
  r2.ran <- sapply(1:length(ind.rand), function(k) {
    ind.rand[[k]]$r2
  })
  aic.ran <- sapply(1:length(ind.rand), function(k) {
    ind.rand[[k]]$aic
  })
  ds.ran <- sapply(1:length(ind.rand), function(k) {
    ind.rand[[k]]$ds
  })
  EPS <- sqrt(.Machine$double.eps)
  p.r2 <- (rowSums(r2.ran >= (matrix(r2.obs, nr = nrow(r2.ran), nc = ncol(r2.ran)) - EPS)) + 1) / (ncol(r2.ran) + 1)
  p.aic <- (rowSums(aic.ran <= (matrix(aic.obs, nr = nrow(aic.ran), nc = ncol(aic.ran)) + EPS)) + 1) / (ncol(aic.ran) + 1)

  r2.m.diff <- r2.obs[1] - r2.obs[3]
  r2.c.diff <- r2.obs[2] - r2.obs[4]

  p.m.r2.diff <- sum((r2.ran[1, ] - r2.ran[3, ]) <= r2.m.diff) / ncol(r2.ran)
  p.c.r2.diff <- sum((r2.ran[2, ] - r2.ran[4, ]) <= r2.c.diff) / ncol(r2.ran)

  if (ds.obs > 0) {
    p.perm <- (sum(ds.ran >= (ds.obs - EPS)) + 1) / (length(ds.ran) + 1)
  } else {
    p.perm <- (sum(ds.ran <= (ds.obs + EPS)) + 1) / (length(ds.ran) + 1)
  }
  p.values <- rbind(
    matrix(p.r2, 2, 2), matrix(p.aic, 2, 2), c(p.perm, NA),
    c(p.m.r2.diff, NA), c(p.c.r2.diff, NA)
  )
  rownames(p.values) <- c(
    "P.R2M", "P.R2C", "P.AIC1", "P.AIC2", "P.dS.perm",
    "P.R2M.diff", "P.R2C.diff"
  )
  output <- rbind(tdci, p.values)
  output <- as.data.frame(output)
  output
}
betai <- Dist.table[, c(1, 3)]
treat <- Dist.table
treat$Layer <- ifelse(treat$Layer %in% c("L1", "L2"), "Topsoil", "Subsoil")
results.lmm <- tdcm_index(
  betai = betai, treat = treat,
  scale.num = F, rand = 1000
)
write.csv(results.lmm, file = paste0("beta_LMM_perm test", ".csv"))
