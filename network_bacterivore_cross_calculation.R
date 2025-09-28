# Big networks ####
## bacteria&protists ####
new_assmatrix = function (comm, method = c("spearman", "pearson", "kendall", 
                                           "clr_pearson", "dcor", "bray", "lss", "mi", "mic"), majority = 0.5, 
                          missing.data.fill = c("fill.pair", "fill.all", "blank", 
                                                "fill.pair.BL", "fill.all.BL"), fillzero.value = 0.01, 
                          logarithm = TRUE, samp.time = NULL, time.lag = 1, silent = FALSE, 
                          mklthread = 4, output.filled.matrix = FALSE, CLR.transform = FALSE) 
{
  if (!missing.data.fill %in% c("fill.pair", "fill.all", "blank", 
                                "fill.pair.BL", "fill.all.BL")) {
    stop("missing.data.fill is wrong.")
  }
  zero.as.na = TRUE
  if (sum(comm < 0, na.rm = TRUE) > 0) 
    stop("Negative value is not allowed.")
  if (zero.as.na) {
    comm[which(comm == 0, arr.ind = TRUE)] = NA
  }
  else {
    if (sum(comm[!is.na(comm)] == 0) > 0 & sum(is.na(comm)) > 
        0) 
      warning("Both zero and NA appeared. NA was treated as missing data, zero was not.")
  }
  if (!is.null(samp.time)) {
    sampint = intersect(rownames(comm), rownames(samp.time))
    if (length(sampint) == 0) 
      stop("samp.time should have the same rownames as comm.")
    if (length(sampint) < nrow(comm)) 
      warning("some rownames of comm are not found in samp.time input.")
    comm = comm[match(sampint, rownames(comm)), , drop = FALSE]
    samp.time = samp.time[match(sampint, rownames(samp.time)), 
                          , drop = FALSE]
  }
  if (majority > 1) {
    min.samp = majority
  }
  else {
    min.samp = nrow(comm) * majority
  }
  if (min.samp < 3) 
    stop("the mininal number should not be less than three.")
  missing.data.fill = missing.data.fill[1]
  sampn = nrow(comm)
  if (majority <= 1) {
    keep.id = which(colSums(comm > 0, na.rm = TRUE) >= (majority * 
                                                          sampn))
  }
  else {
    keep.id = which(colSums(comm > 0, na.rm = TRUE) >= majority)
  }
  if (length(keep.id) == 0) 
    stop("No taxa passed majority rule.")
  comm.old = comm
  comm = comm.old[, keep.id, drop = FALSE]
  if (!silent) 
    message("Taxa number passing majority: ", ncol(comm), 
            ". ", date())
  if (!is.null(samp.time)) {
    obj.lev = unique(samp.time[, 1])
    idlag1 <- idlag2 <- list()
    for (i in 1:length(obj.lev)) {
      idi = which(samp.time[, 1] == obj.lev[i])
      sampi = rownames(samp.time)[idi]
      timei = samp.time[idi, 2]
      ti.lev = sort(unique(timei))
      if ((length(ti.lev) - time.lag) < 3) 
        stop("overlap lower than three. not applicable.")
      t1 = list()
      t2 = list()
      for (j in 1:(length(ti.lev) - time.lag)) {
        t1j = which(timei == ti.lev[j])
        t2j = which(timei == ti.lev[j + time.lag])
        if (length(t1j) != 1 | length(t2j) != 1) 
          warning("time point not unique, the original rank matters.")
        if (length(t1j) != length(t2j)) 
          warning("some time point number not match after lagged.")
        if (length(t1j) > length(t2j)) {
          t1[[j]] = t1j
          t2[[j]] = t1j
          t2[[j]][] = t2j
        }
        else {
          t2[[j]] = t2j
          t1[[j]] = t2j
          t1[[j]][] = t1j
        }
      }
      idlag1[[i]] = idi[unlist(t1)]
      idlag2[[i]] = idi[unlist(t2)]
    }
    idlagA = unlist(idlag1)
    idlagB = unlist(idlag2)
  }
  meth = method[1]
  minab = min(comm[comm > 0], na.rm = TRUE)
  xm = comm/minab
  if (missing.data.fill %in% c("fill.pair.BL", "fill.all.BL")) {
    require(zCompositions)
    xm0 = xm
    xm0[is.na(xm)] = 0
    xmn0 = zCompositions::cmultRepl(xm0, label = 0, method = "BL", 
                                    delta = 0.65, z.delete = FALSE,
                                    suppress.print = FALSE,, output = "p-counts")
  }
  else if (missing.data.fill %in% c("fill.pair", "fill.all")) {
    xmn0 = xm
    xmn0[is.na(xm)] = fillzero.value
  }
  else {
    xmn0 = xm
  }
  if (logarithm & (meth != "clr_pearson") & (!CLR.transform)) {
    xm = log(xm)
    if (missing.data.fill != "blank") {
      xmn0 = log(xmn0)
    }
    else {
      xmn0 = xm
    }
  }
  else if (meth == "clr_pearson" | CLR.transform) {
    x_log = log(xm)
    x_lnam = rowMeans(x_log, na.rm = TRUE)
    x_clr = log(xmn0) - x_lnam
    xmn0 = x_clr
    xm[which(!is.na(xm), arr.ind = TRUE)] = x_clr[which(!is.na(xm), 
                                                        arr.ind = TRUE)]
  }
  nr = nrow(xm)
  nc = ncol(xm)
  xp = expand.grid(1:nc, 1:nc)
  xp = xp[xp$Var1 > xp$Var2, ]
  xcor = vector("double")
  fillmiss <- function(idx, idy, idrx = NULL, idry = NULL, 
                       ...) {
    if (is.null(idrx)) {
      idrx = 1:nrow(xm)
    }
    if (is.null(idry)) {
      idry = 1:nrow(xm)
    }
    x = xm[idrx, idx]
    y = xm[idry, idy]
    if (missing.data.fill %in% c("fill.pair", "fill.pair.BL")) {
      x[which(is.na(x) & (!is.na(y)))] = xmn0[idrx[which(is.na(x) & 
                                                           (!is.na(y)))], idx]
      y[which(is.na(y) & (!is.na(x)))] = xmn0[idry[which(is.na(y) & 
                                                           (!is.na(x)))], idy]
    }
    else if (missing.data.fill == "fill.all") {
      x = xmn0[idrx, idx]
      y = xmn0[idry, idy]
    }
    if (sum(!is.na(x) & !is.na(y)) < 3) {
      x[] = NA
      y[] = NA
    }
    else {
      keep.id = which(!is.na(x) & !is.na(y))
      x = x[keep.id]
      y = y[keep.id]
    }
    list(x = x, y = y)
  }
  xpnum = nrow(xp)
  trac = seq(from = 1, to = xpnum, by = 2000)
  if (is.null(samp.time)) {
    if (!silent) 
      message("Correlation calculation without time lagging begin.", 
              date())
    if (meth == "pearson" || meth == "spearman" || meth == 
        "kendall") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        stats::cor(xy$x, xy$y, method = meth)
      })
    }
    else if (meth == "clr_pearson") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        stats::cor(xy$x, xy$y, method = "pearson")
      })
    }
    else if (meth == "dcor") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        if (length(xy$x) > 0) {
          out = energy::dcor(xy$x, xy$y, 1)
        }
        else {
          out = NA
        }
        out
      })
    }
    else if (meth == "bray") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        if (length(xy$x) > 0) {
          vec = rbind(xy$x, xy$y)
          out = 1 - (vegan::vegdist(vec, method = "bray")[[1]])
        }
        else {
          out = NA
        }
        out
      })
    }
    else if (meth == "lss") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        for (i in 1:nrow(xp)) {
          xy = fillmiss(idx = xp[i, 1], idy = xp[i, 
                                                 2])
          if (length(xy$x) > 0) {
            tr1 = normalTransform(xy$x)
            tr2 = normalTransform(xy$y)
            lss_vec = LocalSimilarity3(tr1, tr2, 0, 
                                       length(tr1), F)
            score = lss_vec[1]
            sign = lss_vec[5]
            if (sign != 1) {
              sign = -1
            }
            out = score * sign
          }
          else {
            out = NA
          }
          out
        }
      })
    }
    else if (meth == "mi") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        if (length(xy$x) > 0) {
          xymm = cbind(xy$x, xy$y)
          dvec = infotheo::discretize(xymm, nbins = ceiling(1 + 
                                                              log(nr, 2)))
          out = infotheo::mutinformation(dvec)[1, 2]
        }
        else {
          out = NA
        }
        out
      })
    }
    else if (meth == "mic") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        if (length(xy$x) > 0) {
          vec = cbind(xy$x, xy$y)
          out = minerva::mine(vec, alpha = 0.65)$MIC[1, 
                                                     2]
        }
        else {
          out = NA
        }
        out
      })
    }
    else {
      print("type of method not suppoted yet!\n")
    }
  }
  if (!is.null(samp.time)) {
    if (!silent) 
      message("Correlation calculation with time lagging begin.", 
              date())
    if (meth == "pearson" || meth == "spearman" || meth == 
        "kendall") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        cori = vector("double")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        cori[1] = stats::cor(xy$x, xy$y, method = meth)
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagA, idry = idlagB)
        cori[2] = stats::cor(xy$x, xy$y, method = meth)
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagB, idry = idlagA)
        cori[3] = stats::cor(xy$x, xy$y, method = meth)
        if (sum(!is.na(cori)) == 0) {
          out = NA
        }
        else {
          out = max(cori, na.rm = TRUE)
        }
        out
      })
    }
    else if (meth == "clr_pearson") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        cori = vector("double")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        cori[1] = stats::cor(xy$x, xy$y, method = "pearson")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagA, idry = idlagB)
        cori[2] = stats::cor(xy$x, xy$y, method = "pearson")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagB, idry = idlagA)
        cori[3] = stats::cor(xy$x, xy$y, method = "pearson")
        if (sum(!is.na(cori)) == 0) {
          out = NA
        }
        else {
          out = max(cori, na.rm = TRUE)
        }
        out
      })
    }
    else if (meth == "dcor") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        cori = vector("double")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        if (length(xy$x) > 0) {
          cori[1] = energy::dcor(xy$x, xy$y, 1)
        }
        else {
          cori[1] = NA
        }
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagA, idry = idlagB)
        if (length(xy$x) > 0) {
          cori[2] = energy::dcor(xy$x, xy$y, 1)
        }
        else {
          cori[2] = NA
        }
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagB, idry = idlagA)
        if (length(xy$x) > 0) {
          cori[3] = energy::dcor(xy$x, xy$y, 1)
        }
        else {
          cori[3] = NA
        }
        if (sum(!is.na(cori)) == 0) {
          out = NA
        }
        else {
          out = max(cori, na.rm = TRUE)
        }
        out
      })
    }
    else if (meth == "bray") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        cori = vector("double")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        if (length(xy$x) > 0) {
          vec = rbind(xy$x, xy$y)
          cori[1] = 1 - (vegan::vegdist(vec, method = "bray")[[1]])
        }
        else {
          cori[1] = NA
        }
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagA, idry = idlagB)
        if (length(xy$x) > 0) {
          vec = rbind(xy$x, xy$y)
          cori[2] = 1 - (vegan::vegdist(vec, method = "bray")[[1]])
        }
        else {
          cori[2] = NA
        }
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagB, idry = idlagA)
        if (length(xy$x) > 0) {
          vec = rbind(xy$x, xy$y)
          cori[3] = 1 - (vegan::vegdist(vec, method = "bray")[[1]])
        }
        else {
          cori[3] = NA
        }
        if (sum(!is.na(cori)) == 0) {
          out = NA
        }
        else {
          out = max(cori, na.rm = TRUE)
        }
        out
      })
    }
    else if (meth == "lss") {
      xylss <- function(xy, ...) {
        if (length(xy$x) > 0) {
          tr1 = normalTransform(xy$x)
          tr2 = normalTransform(xy$y)
          lss_vec = LocalSimilarity3(tr1, tr2, 0, nr, 
                                     F)
          score = lss_vec[1]
          sign = lss_vec[5]
          if (sign != 1) {
            sign = -1
          }
          out = score * sign
        }
        else {
          out = NA
        }
        out
      }
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        cori = vector("double")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        cori[1] = xylss(xy)
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagA, idry = idlagB)
        cori[2] = xylss(xy)
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagB, idry = idlagA)
        cori[3] = xylss(xy)
        if (sum(!is.na(cori)) == 0) {
          out = NA
        }
        else {
          out = max(cori, na.rm = TRUE)
        }
        out
      })
    }
    else if (meth == "mi") {
      xymi <- function(xy, ...) {
        if (length(xy$x) > 0) {
          xymm = cbind(xy$x, xy$y)
          dvec = infotheo::discretize(xymm, nbins = ceiling(1 + 
                                                              log(nr, 2)))
          out = infotheo::mutinformation(dvec)[1, 2]
        }
        else {
          out = NA
        }
        out
      }
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        cori = vector("double")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        cori[1] = xymi(xy)
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagA, idry = idlagB)
        cori[2] = xymi(xy)
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagB, idry = idlagA)
        cori[3] = xymi(xy)
        if (sum(!is.na(cori)) == 0) {
          out = NA
        }
        else {
          out = max(cori, na.rm = TRUE)
        }
        out
      })
    }
    else if (meth == "mic") {
      xcor = sapply(1:nrow(xp), function(i) {
        if ((!silent) & (i %in% trac)) 
          message("Now i=", i, " in ", xpnum, ". ", 
                  date())
        cori = vector("double")
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2])
        if (length(xy$x) > 0) {
          vec = cbind(xy$x, xy$y)
          cori[1] = minerva::mine(vec, alpha = 0.65)$MIC[1, 
                                                         2]
        }
        else {
          cori[1] = NA
        }
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagA, idry = idlagB)
        if (length(xy$x) > 0) {
          vec = cbind(xy$x, xy$y)
          cori[2] = minerva::mine(vec, alpha = 0.65)$MIC[1, 
                                                         2]
        }
        else {
          cori[2] = NA
        }
        xy = fillmiss(idx = xp[i, 1], idy = xp[i, 2], 
                      idrx = idlagB, idry = idlagA)
        if (length(xy$x) > 0) {
          vec = cbind(xy$x, xy$y)
          cori[3] = minerva::mine(vec, alpha = 0.65)$MIC[1, 
                                                         2]
        }
        else {
          cori[3] = NA
        }
        if (sum(!is.na(cori)) == 0) {
          out = NA
        }
        else {
          out = max(cori, na.rm = TRUE)
        }
        out
      })
    }
    else {
      print("type of method not suppoted yet!\n")
    }
  }
  mcor = matrix(rep(0, nc * nc), nrow = nc)
  mcor[lower.tri(mcor, diag = FALSE)] = xcor
  mcor = mcor + t(mcor)
  if (meth == "lss" || meth == "mi") {
    mcor_zero_diag = mcor
    diag(mcor_zero_diag) = 0
    mcor = mcor/max(abs(mcor_zero_diag))
  }
  diag(mcor) = 1
  if (sum(is.na(mcor)) > 0) 
    warning("Some correlations returned NA, set as zero.")
  mcor[is.na(mcor)] = 0
  rownames(mcor) <- colnames(mcor) <- colnames(comm)
  if (output.filled.matrix) {
    output = list(ass.matrix = mcor, filled.matrix = xmn0)
  }
  else {
    output = mcor
  }
  output
}

# setwd("C:\\Users\\True\\OneDrive\\桌面")
# com1.f = "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
# com2.f <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"
# clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\18S_tax.txt"
# treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"
# source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\findcutoff2.r")
# source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\brodyrmt2.r")
# source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\assmcut2.r")

setwd("/home/fangke1/Duolun/RMT_network/BP")
com1.f = "/home/fangke1/Duolun/RMT_network/data for use/bacteria.zotu245_resample_35000.txt"
com2.f <- "/home/fangke1/Duolun/RMT_network/data for use/protist_zotu_2868.txt"
clas.file = "/home/fangke1/Duolun/RMT_network/data for use/18S_tax.txt"
treat.file <- "/home/fangke1/Duolun/RMT_network/data for use/treatment.csv"
source("/home/fangke1/Duolun/RMT_network/findcutoff2.r")
source("/home/fangke1/Duolun/RMT_network/brodyrmt2.r")
source("/home/fangke1/Duolun/RMT_network/assmcut2.r")

prefix1 <- "Bac"
prefix2 <- "Pro"
cor.meth <- "pearson"
maj.1 <- 0.75
maj.2 <- 0.25

### find cutoff2 ####
pre.treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
pre.treat <- subset(pre.treat, plant.type == "TS")

pre.comm1 <- t(read.table(com1.f,
                          header = TRUE, sep = "\t", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE
))
name.row <- rownames(pre.comm1)
name.row[which(rownames(pre.comm1) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pre.comm1) == "YD52L3")] <- "YD3L3"
rownames(pre.comm1) <- name.row

pre.comm2 <- t(read.table(com2.f,
                          header = TRUE, sep = "\t", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE
))
name.row <- rownames(pre.comm2)
name.row[which(rownames(pre.comm2) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pre.comm2) == "YD52L3")] <- "YD3L3"
rownames(pre.comm2) <- name.row

pre.clas <- read.table(clas.file,
                       header = TRUE, sep = "\t", row.names = 1,
                       as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                       check.names = FALSE
)

pre.clas$main_functional_class = ifelse(pre.clas$main_functional_class == "predator (add)", "predator", pre.clas$main_functional_class)
pre.clas = pre.clas[pre.clas$main_functional_class == "predator" & grepl("bacteria", pre.clas$associated_organism),]

library(ieggr)
sampc <- match.name(cn.list = list(pre.comm2 = pre.comm2), rn.list = list(pre.clas = pre.clas), silent = T)
pre.comm2 <- sampc$pre.comm2
pre.clas <- sampc$pre.clas

sampc <- match.name(rn.list = list(pre.treat = pre.treat, pre.comm1 = pre.comm1, pre.comm2 = pre.comm2), silent = T)
dim(pre.comm1)
dim(pre.comm2)
pre.comm1 <- sampc$pre.comm1
pre.comm2 <- sampc$pre.comm2
pre.treat <- sampc$pre.treat
pre.comm1 <- pre.comm1[, colSums(pre.comm1) > 0]
pre.comm2 <- pre.comm2[, colSums(pre.comm2) > 0]
dim(pre.comm1)
dim(pre.comm2)

# for (j in c("L1", "L2", "L3", "L4")) {
for (j in c("L2", "L3", "L4")) {
  prefixm <- paste0(j, ".", cor.meth)
  library(ieggr)
  treat = pre.treat
  comm1 = pre.comm1
  comm2 = pre.comm2
  
  treat <- subset(treat, Layer == j)
  
  match.rowname <- match.name(rn.list = list(comm1 = comm1, treat = treat), silent = T)
  comm1 <- match.rowname$comm1
  treat <- match.rowname$treat
  
  match.rowname <- match.name(rn.list = list(comm2 = comm2, treat = treat), silent = T)
  comm2 <- match.rowname$comm2
  treat <- match.rowname$treat
  
  colnames(comm1) <- paste0(prefix1, colnames(comm1))
  colnames(comm2) <- paste0(prefix2, colnames(comm2))
  # filter abundance by majority
  # if not filter here, there will be "Error in cor.test.default(com1cn0[idij,i], com2cn0[idij,j], method = cor.method):
  # not enough finite observations" in 3.1
  {
    dim(comm1)
    comm1 <- comm1[, which((colSums(comm1 > 0) / dim(comm1)[1]) >= maj.1)]
    dim(comm1)
    
    dim(comm2)
    comm2 <- comm2[, which((colSums(comm2 > 0) / dim(comm2)[1]) >= maj.2)]
    dim(comm2)
  }
  
  require(MENA)
  # 2.1 comm1 within kingdom
  assmx1 <- new_assmatrix(comm1,
                          method = cor.meth, majority = maj.1,
                          missing.data.fill = "fill.pair.BL", logarithm = TRUE,
                          fillzero.value = 0.01, samp.time = NULL, time.lag = 0,
                          silent = FALSE, mklthread = 4, output.filled.matrix = FALSE,
                          CLR.transform = TRUE
  )
  
  cutoff.use1 <- 0.747
  assmc1 <- MENA::assmcut(ass.matrix = assmx1, cutoff = cutoff.use1)
  
  # 2.2 comm2 within kingdom
  assmx2 <- new_assmatrix(comm2,
                          method = cor.meth, majority = maj.2,
                          missing.data.fill = "fill.pair.BL", logarithm = TRUE,
                          fillzero.value = 0.01, samp.time = NULL, time.lag = 0,
                          silent = FALSE, mklthread = 4, output.filled.matrix = FALSE,
                          CLR.transform = TRUE
  )
  cutoff.use2 <- 0.86
  assmc2 <- MENA::assmcut(ass.matrix = assmx2, cutoff = cutoff.use2)
  
  # 3.1 association matrix between two kingdoms: correlation matrix
  new.CLRbipcor <- function(comm1, comm2,
                            cor.method = "spearman",
                            prefixi = prefixi, save.wd = getwd()) {
    CLRtrsf <- function(xm) {
      xm[xm == 0] <- NA
      xm0 <- xm
      xm0[is.na(xm0)] <- 0
      xmn0 <- zCompositions::cmultRepl(xm0, label = 0, method = "BL", 
                                       delta = 0.65, z.delete = FALSE,
                                       suppress.print = FALSE,, output = "p-counts")
      x_log <- log(xm)
      x_lnam <- rowMeans(x_log, na.rm = TRUE)
      x_clr <- log(xmn0) - x_lnam
      xm[which(!is.na(xm), arr.ind = TRUE)] <- x_clr[which(!is.na(xm), arr.ind = TRUE)]
      list(xm = xm, xmn0 = x_clr)
    }
    
    com1.clr <- CLRtrsf(comm1)
    com1cna <- com1.clr$xm
    com1cn0 <- com1.clr$xmn0
    com2.clr <- CLRtrsf(comm2)
    com2cna <- com2.clr$xm
    com2cn0 <- com2.clr$xmn0
    
    trck <- seq(from = 1, to = ncol(comm1), by = 100)
    corm12 <- sapply(1:ncol(comm1), function(i) {
      if (i %in% trck) {
        message("i=", i, " in ", ncol(comm1), ". ", date())
      }
      sapply(1:ncol(comm2), function(j) {
        idij <- which(!(is.na(com1cna[, i]) & is.na(com2cna[, j])))
        ctij <- stats::cor.test(com1cn0[idij, i], com2cn0[idij, j], method = cor.method)
        c(r = ctij$estimate[[1]], p = ctij$p.value)
      })
    }, simplify = "array")
    dim(corm12)
    dimnames(corm12)[[2]] <- colnames(comm2)
    dimnames(corm12)[[3]] <- colnames(comm1)
    corm12.r <- corm12[1, , ]
    corm12.p <- corm12[2, , ]
    corm12.padj <- matrix(p.adjust(as.vector(corm12.p), method = "fdr"), nrow = nrow(corm12.p), ncol = ncol(corm12.p))
    rownames(corm12.padj) <- rownames(corm12.p)
    colnames(corm12.padj) <- colnames(corm12.p)
    
    dim(corm12.r)
    dim(corm12.p)
    dim(corm12.padj)
    corm12.r[1:5, 1:5]
    corm12.p[1:5, 1:5]
    corm12.padj[1:5, 1:5]
    
    # save.file(corm12.r, prefix = prefixi, filename = "CorCoef", rowname.title = "TaxID", folder = save.wd)
    # save.file(corm12.p, prefix = prefixi, filename = "Pvalue", rowname.title = "TaxID", folder = save.wd)
    # save.file(corm12.padj, prefix = prefixi, filename = "Padjust", rowname.title = "TaxID", folder = save.wd)
    
    list(
      rho = corm12.r, p = corm12.p, p.adj = corm12.padj,
      stats = c(
        ra0.3 = sum(abs(corm12.r) > 0.3),
        ra0.5 = sum(abs(corm12.r) > 0.5),
        pu0.05 = sum(corm12.p < 0.05),
        padju0.05 = sum(corm12.padj < 0.05)
      )
    )
  }
  prefixi <- paste0(prefixm, ".All.", prefix1, ".", maj.1, ".", prefix2, ".", maj.2)
  cbcor.all <- new.CLRbipcor(comm1, comm2,
                             cor.method = cor.meth,
                             prefixi = prefixi, save.wd = getwd()
  )
  cbcor.all$stats
  assmbt <- cbcor.all$rho
  # 3.2 combined association matrix
  assm12 <- matrix(NA, nrow = (ncol(comm1) + ncol(comm2)), ncol = (ncol(comm1) + ncol(comm2)))
  rownames(assm12) <- colnames(assm12) <- c(colnames(comm1), colnames(comm2))
  assm12[match(rownames(assmc1), rownames(assm12)), match(colnames(assmc1), colnames(assm12))] <- assmc1
  assm12[match(rownames(assmc2), rownames(assm12)), match(colnames(assmc2), colnames(assm12))] <- assmc2
  assm12[match(rownames(assmbt), rownames(assm12)), match(colnames(assmbt), colnames(assm12))] <- assmbt
  assm12[match(colnames(assmbt), rownames(assm12)), match(rownames(assmbt), colnames(assm12))] <- t(assmbt)
  assm12[is.na(assm12)] <- 0
  
  id.const <- rbind(
    as.matrix(expand.grid(list(match(colnames(comm1), rownames(assm12)), match(colnames(comm1), colnames(assm12))))),
    as.matrix(expand.grid(list(match(colnames(comm2), rownames(assm12)), match(colnames(comm2), colnames(assm12)))))
  )
  
  brmt12 <- brodyrmt2(ass.matrix = assm12, id.const = id.const, nthread = 100)
  # 3.3 search the cutoff based on RMT
  cutoff12out <- findcutoff2(brmt12, outputxy = TRUE, criteria = c(0.1, 0.5, 0.9))
  cutoff12 <- cutoff12out$cutoff
  smthxy12 <- cutoff12out$smth_xy
  cutoff.use12 <- cutoff12[1, 1]
  assmc12 <- assmcut2(ass.matrix = assm12, cutoff = cutoff.use12, id.const = id.const)
  
  output3 <- list(assmx = assm12, brmt = brmt12, cutoff = cutoff12, smthxy = smthxy12, cutoff.use = cutoff.use12, assmc = assmc12)
  save(output3, file = paste0(prefixm, ".", prefix1, ".", maj.1, ".", prefix2, ".", maj.2, ".CLR.FillPairLB.rda"))
}

### after finding cutoff ####
pre.treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
pre.treat <- subset(pre.treat, plant.type == "TS")

pre.comm1 <- t(read.table(com1.f,
                          header = TRUE, sep = "\t", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE
))
name.row <- rownames(pre.comm1)
name.row[which(rownames(pre.comm1) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pre.comm1) == "YD52L3")] <- "YD3L3"
rownames(pre.comm1) <- name.row

pre.comm2 <- t(read.table(com2.f,
                          header = TRUE, sep = "\t", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE
))
name.row <- rownames(pre.comm2)
name.row[which(rownames(pre.comm2) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pre.comm2) == "YD52L3")] <- "YD3L3"
rownames(pre.comm2) <- name.row

pre.clas <- read.table(clas.file,
                       header = TRUE, sep = "\t", row.names = 1,
                       as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                       check.names = FALSE
)

pre.clas$main_functional_class = ifelse(pre.clas$main_functional_class == "predator (add)", "predator", pre.clas$main_functional_class)
pre.clas = pre.clas[pre.clas$main_functional_class == "predator" & grepl("bacteria", pre.clas$associated_organism),]

sampc <- match.name(cn.list = list(pre.comm2 = pre.comm2), rn.list = list(pre.clas = pre.clas), silent = T)
pre.comm2 <- sampc$pre.comm2
pre.clas <- sampc$pre.clas

sampc <- match.name(rn.list = list(pre.treat = pre.treat, pre.comm1 = pre.comm1, pre.comm2 = pre.comm2), silent = T)
dim(pre.comm1)
dim(pre.comm2)
pre.comm1 <- sampc$pre.comm1
pre.comm2 <- sampc$pre.comm2
pre.treat <- sampc$pre.treat
pre.comm1 <- pre.comm1[, colSums(pre.comm1) > 0]
pre.comm2 <- pre.comm2[, colSums(pre.comm2) > 0]
dim(pre.comm1)
dim(pre.comm2)

for (j in c("L1", "L2", "L3", "L4")) {
  prefixm <- paste0(j, ".", cor.meth)
  library(ieggr)
  treat = pre.treat
  comm1 = pre.comm1
  comm2 = pre.comm2
  
  treat <- subset(treat, Layer == j)
  
  match.rowname <- match.name(rn.list = list(comm1 = comm1, treat = treat), silent = T)
  comm1 <- match.rowname$comm1
  treat <- match.rowname$treat
  
  match.rowname <- match.name(rn.list = list(comm2 = comm2, treat = treat), silent = T)
  comm2 <- match.rowname$comm2
  treat <- match.rowname$treat
  
  colnames(comm1) <- paste0(prefix1, colnames(comm1))
  colnames(comm2) <- paste0(prefix2, colnames(comm2))
  # filter abundance by majority
  # if not filter here, there will be "Error in cor.test.default(com1cn0[idij,i], com2cn0[idij,j], method = cor.method):
  # not enough finite observations" in 3.1
  {
    dim(comm1)
    comm1 <- comm1[, which((colSums(comm1 > 0) / dim(comm1)[1]) >= maj.1)]
    dim(comm1)
    
    dim(comm2)
    comm2 <- comm2[, which((colSums(comm2 > 0) / dim(comm2)[1]) >= maj.2)]
    dim(comm2)
  }
  
  require(MENA)
  # 2.1 comm1 within kingdom
  assmx1 <- new_assmatrix(comm1,
                          method = cor.meth, majority = maj.1,
                          missing.data.fill = "fill.pair.BL", logarithm = TRUE,
                          fillzero.value = 0.01, samp.time = NULL, time.lag = 0,
                          silent = FALSE, mklthread = 4, output.filled.matrix = FALSE,
                          CLR.transform = TRUE
  )
  
  cutoff.use1 <- 0.747
  assmc1 <- MENA::assmcut(ass.matrix = assmx1, cutoff = cutoff.use1)
  
  # 2.2 comm2 within kingdom
  assmx2 <- new_assmatrix(comm2,
                          method = cor.meth, majority = maj.2,
                          missing.data.fill = "fill.pair.BL", logarithm = TRUE,
                          fillzero.value = 0.01, samp.time = NULL, time.lag = 0,
                          silent = FALSE, mklthread = 4, output.filled.matrix = FALSE,
                          CLR.transform = TRUE
  )
  cutoff.use2 <- 0.86
  assmc2 <- MENA::assmcut(ass.matrix = assmx2, cutoff = cutoff.use2)
  
  # 3.1 association matrix between two kingdoms: correlation matrix
  new.CLRbipcor <- function(comm1, comm2,
                            cor.method = "spearman",
                            prefixi = prefixi, save.wd = getwd()) {
    CLRtrsf <- function(xm) {
      xm[xm == 0] <- NA
      xm0 <- xm
      xm0[is.na(xm0)] <- 0
      xmn0 <- zCompositions::cmultRepl(xm0,
                                       label = 0, method = "BL",
                                       z.delete = FALSE,
                                       delta = 0.65, suppress.print = FALSE, output = "p-counts"
      )
      x_log <- log(xm)
      x_lnam <- rowMeans(x_log, na.rm = TRUE)
      x_clr <- log(xmn0) - x_lnam
      xm[which(!is.na(xm), arr.ind = TRUE)] <- x_clr[which(!is.na(xm), arr.ind = TRUE)]
      list(xm = xm, xmn0 = x_clr)
    }
    
    com1.clr <- CLRtrsf(comm1)
    com1cna <- com1.clr$xm
    com1cn0 <- com1.clr$xmn0
    com2.clr <- CLRtrsf(comm2)
    com2cna <- com2.clr$xm
    com2cn0 <- com2.clr$xmn0
    
    trck <- seq(from = 1, to = ncol(comm1), by = 100)
    corm12 <- sapply(1:ncol(comm1), function(i) {
      if (i %in% trck) {
        message("i=", i, " in ", ncol(comm1), ". ", date())
      }
      sapply(1:ncol(comm2), function(j) {
        idij <- which(!(is.na(com1cna[, i]) & is.na(com2cna[, j])))
        ctij <- stats::cor.test(com1cn0[idij, i], com2cn0[idij, j], method = cor.method)
        c(r = ctij$estimate[[1]], p = ctij$p.value)
      })
    }, simplify = "array")
    dim(corm12)
    dimnames(corm12)[[2]] <- colnames(comm2)
    dimnames(corm12)[[3]] <- colnames(comm1)
    corm12.r <- corm12[1, , ]
    corm12.p <- corm12[2, , ]
    corm12.padj <- matrix(p.adjust(as.vector(corm12.p), method = "fdr"), nrow = nrow(corm12.p), ncol = ncol(corm12.p))
    rownames(corm12.padj) <- rownames(corm12.p)
    colnames(corm12.padj) <- colnames(corm12.p)
    
    dim(corm12.r)
    dim(corm12.p)
    dim(corm12.padj)
    corm12.r[1:5, 1:5]
    corm12.p[1:5, 1:5]
    corm12.padj[1:5, 1:5]
    
    # save.file(corm12.r, prefix = prefixi, filename = "CorCoef", rowname.title = "TaxID", folder = save.wd)
    # save.file(corm12.p, prefix = prefixi, filename = "Pvalue", rowname.title = "TaxID", folder = save.wd)
    # save.file(corm12.padj, prefix = prefixi, filename = "Padjust", rowname.title = "TaxID", folder = save.wd)
    
    list(
      rho = corm12.r, p = corm12.p, p.adj = corm12.padj,
      stats = c(
        ra0.3 = sum(abs(corm12.r) > 0.3),
        ra0.5 = sum(abs(corm12.r) > 0.5),
        pu0.05 = sum(corm12.p < 0.05),
        padju0.05 = sum(corm12.padj < 0.05)
      )
    )
  }
  prefixi <- paste0(prefixm, ".All.", prefix1, ".", maj.1, ".", prefix2, ".", maj.2)
  cbcor.all <- new.CLRbipcor(comm1, comm2,
                             cor.method = cor.meth,
                             prefixi = prefixi, save.wd = getwd()
  )
  cbcor.all$stats
  assmbt <- cbcor.all$rho
  # 3.2 combined association matrix
  assm12 <- matrix(NA, nrow = (ncol(comm1) + ncol(comm2)), ncol = (ncol(comm1) + ncol(comm2)))
  rownames(assm12) <- colnames(assm12) <- c(colnames(comm1), colnames(comm2))
  assm12[match(rownames(assmc1), rownames(assm12)), match(colnames(assmc1), colnames(assm12))] <- assmc1
  assm12[match(rownames(assmc2), rownames(assm12)), match(colnames(assmc2), colnames(assm12))] <- assmc2
  assm12[match(rownames(assmbt), rownames(assm12)), match(colnames(assmbt), colnames(assm12))] <- assmbt
  assm12[match(colnames(assmbt), rownames(assm12)), match(rownames(assmbt), colnames(assm12))] <- t(assmbt)
  assm12[is.na(assm12)] <- 0
  
  id.const <- rbind(
    as.matrix(expand.grid(list(match(colnames(comm1), rownames(assm12)), match(colnames(comm1), colnames(assm12))))),
    as.matrix(expand.grid(list(match(colnames(comm2), rownames(assm12)), match(colnames(comm2), colnames(assm12)))))
  )
  
  # 3.3 search the cutoff based on RMT
  cutoff.use12 <- 0.789
  assmc12 <- assmcut2(ass.matrix = assm12, cutoff = cutoff.use12, id.const = id.const)
  
  output3 <- list(assmx = assm12, cutoff.use = cutoff.use12, assmc = assmc12)
  save(output3, file = paste0(prefixm, ".", prefix1, ".", maj.1, ".", prefix2, ".", maj.2, ".CLR.FillPairLB.rda"))
}

# small networks: one sample as one network ####
library(ieggr)
# com1.f = "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
# com2.f <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"
# clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\18S_tax.txt"
# treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

com1.f = "/home/fangke1/Duolun/RMT_network/data for use/bacteria.zotu245_resample_35000.txt"
com2.f <- "/home/fangke1/Duolun/RMT_network/data for use/protist_zotu_2868.txt"
clas.file = "/home/fangke1/Duolun/RMT_network/data for use/18S_tax.txt"
treat.file <- "/home/fangke1/Duolun/RMT_network/data for use/treatment.csv"

pre.comm1 <- t(read.table(com1.f,
                          header = TRUE, sep = "\t", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE
))
name.row <- rownames(pre.comm1)
name.row[which(rownames(pre.comm1) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pre.comm1) == "YD52L3")] <- "YD3L3"
rownames(pre.comm1) <- name.row

pre.comm2 <- t(read.table(com2.f,
                          header = TRUE, sep = "\t", row.names = 1,
                          as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                          check.names = FALSE
))
name.row <- rownames(pre.comm2)
name.row[which(rownames(pre.comm2) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(pre.comm2) == "YD52L3")] <- "YD3L3"
rownames(pre.comm2) <- name.row

pre.clas <- read.table(clas.file,
                       header = TRUE, sep = "\t", row.names = 1,
                       as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                       check.names = FALSE
)
pre.clas$main_functional_class = ifelse(pre.clas$main_functional_class == "predator (add)", "predator", pre.clas$main_functional_class)
pre.clas = pre.clas[pre.clas$main_functional_class == "predator" & grepl("bacteria", pre.clas$associated_organism),]

pre.treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
pre.treat <- subset(pre.treat, plant.type == "TS")

prefix1 <- "Bac"
prefix2 <- "Pro"
cor.meth <- "pearson"
maj.1 <- 0.75
maj.2 <- 0.25

z = 1
result.list <- list()
for (j in c("L1","L2","L3","L4")) {
  prefixm <- paste0(j, ".", cor.meth)
  library(ieggr)
  
  treat = pre.treat
  comm1 = pre.comm1
  comm2 = pre.comm2
  clas = pre.clas
  
  load(paste0(prefixm, ".", prefix1, ".", maj.1, ".", prefix2, ".", maj.2, ".CLR.FillPairLB.rda"))
  
  treat <- subset(treat, Layer == j)
  
  match.rowname <- match.name(rn.list = list(clas = clas), cn.list = list(comm2 = comm2), silent = T)
  comm2 <- match.rowname$comm2
  clas <- match.rowname$clas
  
  match.rowname <- match.name(rn.list = list(comm1 = comm1, treat = treat), silent = T)
  comm1 <- match.rowname$comm1
  treat <- match.rowname$treat
  
  match.rowname <- match.name(rn.list = list(comm2 = comm2, treat = treat), silent = T)
  comm2 <- match.rowname$comm2
  treat <- match.rowname$treat
  for (k in rownames(treat)) {
    dat <- output3$assmc
    
    used.comm1 <- comm1[rownames(comm1) == k, , drop = F]
    used.comm1 <- used.comm1[, colSums(used.comm1) > 0, drop = F]
    colnames(used.comm1) <- paste0(prefix1, colnames(used.comm1))
    
    used.comm2 <- comm2[rownames(comm2) == k, , drop = F]
    used.comm2 <- used.comm2[, colSums(used.comm2) > 0, drop = F]
    colnames(used.comm2) <- paste0(prefix2, colnames(used.comm2))
    
    save.spec <- c(colnames(used.comm1), colnames(used.comm2))
    
    dat <- dat[rownames(dat) %in% save.spec, colnames(dat) %in% save.spec]
    dat <- as.matrix(dat)
    dat[is.na(dat)] <- 0
    
    # remove isolated node
    dim(dat)
    dat <- dat[, colSums(abs(dat)) > 0]
    dat <- dat[rowSums(abs(dat)) > 0, ]
    dim(dat)
    
    result.list[[z]] <- dat
    names(result.list)[z] <- k
    z <- z + 1
  }
}
saveRDS(result.list, paste0(prefix1, "_",maj.1, "_", prefix2, "_", maj.2,".RDS"))

library(bipartite)
library(magrittr)
# source("C:\\Users\\True\\OneDrive\\桌面\\research\\analysis methods\\R.code\\new bipartite index.R")
source("/home/fangke1/Duolun/RMT_network/new bipartite index.R")

result.data <- data.frame(stringsAsFactors = FALSE)
for (z in 1:length(result.list)) {
  cor.data <- result.list[[z]]
  # as the correlation matrix includes the correlations within kingdom, we need to remove them
  cor.data <- cor.data[grepl(prefix2, rownames(cor.data)), ]
  cor.data <- cor.data[, grepl(prefix1, colnames(cor.data))]
  
  cor.data <- as.matrix(cor.data)
  cor.data[is.na(cor.data)] <- 0
  
  # remove isolated node
  cor.data <- cor.data[rowSums(abs(cor.data)) > 0, ]
  cor.data <- cor.data[, colSums(abs(cor.data)) > 0]
  
  # add sample column
  
  result.vector = c(names(result.list)[z],
                    # "Total links"
                    sum(cor.data != 0),
                    # positive links
                    sum(cor.data > 0),
                    # negative links
                    sum(cor.data < 0),
                    # links ratio of Pos to Neg
                    sum(cor.data > 0) / sum(cor.data < 0),
                    # "Total nodes"
                    sum(dim(cor.data)),
                    # high nodes
                    nrow(cor.data),
                    # low nodes
                    ncol(cor.data),
                    # nodes ratio of high to low
                    nrow(cor.data) / ncol(cor.data))
  
  # change the correlation matrix to binary matrix to calculate the average degree
  cor.data <- ifelse(cor.data != 0, 1, 0)
  result.vector = c(result.vector,
                    # average degree all
                    cor.data %>%
                      sum() %>%
                      {
                        # because one link is shared by two nodes, so the links should be multiplied by 2
                        2 * .
                      } %>%
                      {
                        . / sum(dim(cor.data))
                      },
                    # average degree low level
                    cor.data %>%
                      sum() %>%
                      {
                        . / ncol(cor.data)
                      },
                    # average degree high level
                    cor.data %>%
                      sum() %>%
                      {
                        . / nrow(cor.data)
                      })
  result.vector = c(result.vector,
                    new.bipartite.index(cor.data, 
                                        complexity = T,
                                        group.number = T, 
                                        modular = F))
  result.data <- rbind(result.data, result.vector)
  for (i in 1:ncol(result.data)) {
    result.data[, i] <- as.character(result.data[, i])
  }
}
colnames(result.data) <- colnames(result.data) <- c(
  "Sample", "Total links", "Positive links", "Negative links", "Links ratio of Pos to Neg", "Total nodes",
  "High nodes", "Low nodes", "Nodes ratio of high to low", "Average degree all", "Average degree low",
  "Average degree high",
  "Density","Connectance", 
  # "web.asymmetry",
  "linkage.density", "weighted.connectance",
  "modularity", "nestedness", "NODF","H2",
  "number.of.species.HL", "number.of.species.LL", 
  # "mean.number.of.shared.partners.HL", "mean.number.of.shared.partners.LL", 
  "cluster.coefficient.HL", "cluster.coefficient.LL", 
  "weighted.cluster.coefficient.HL", "weighted.cluster.coefficient.LL", 
  "niche.overlap.HL", "niche.overlap.LL", 
  # "togetherness.HL", 
  # "togetherness.LL", "C.score.HL", "C.score.LL", "V.ratio.HL", 
  # "V.ratio.LL", "discrepancy.HL", "discrepancy.LL", 
  # "extinction.slope.HL", "extinction.slope.LL", "robustness.HL", 
  # "robustness.LL", 
  "functional.complementarity.HL", 
  "functional.complementarity.LL", 
  "partner.diversity.HL", "partner.diversity.LL"
  # "generality.HL", "vulnerability.LL"
)

write.csv(result.data, paste0(prefix1, "_", maj.1, "_", prefix2, "_", maj.2, "_",cor.meth,".csv"))

