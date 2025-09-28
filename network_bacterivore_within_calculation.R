setwd("C:\\Users\\True\\OneDrive\\桌面")
# global network calculation ####
library(ieggr)
require(MENA)
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
com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"
prefix.m <- "Protist"
prefix1 <- "Pro"

clas.file = "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\18S\\all samples\\Unoise\\18S_tax.txt"

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

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
clas$main_functional_class = ifelse(clas$main_functional_class == "predator (add)", "predator", clas$main_functional_class)
clas = clas[clas$main_functional_class == "predator" & grepl("bacteria", clas$associated_organism),]

treat <- read.csv(treat.file, header = T, row.names = 1)
treat <- subset(treat, plant.type == "TS")
prefix.m <- paste0(prefix.m, "_", "TS")

sampc <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), silent = T)
comm <- sampc$comm
clas <- sampc$clas

sampc <- match.name(rn.list = list(treat = treat, comm = comm), silent = T)
dim(comm)
comm <- sampc$comm
treat <- sampc$treat
comm <- comm[, colSums(comm) > 0]
dim(comm)

cor.meth <- "pearson"
maj.1 <- 0.25

prefixm <- cor.meth

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
## find cutoff ####
  for (j in c("L1", "L2", "L3", "L4")) {
    used.comm <- comm
    used.treat <- subset(treat, Layer == j)
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
    
    assmx1 <- new_assmatrix(used.comm,
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

