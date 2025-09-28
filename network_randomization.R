setwd("/home/fangke1/Duolun/RMT_network/Fungi")
# The function "net.rand" can be used to generate random networks and 
# compare them with empirical networks
net.pro = function (assmc) 
{
  require(igraph)
  amg = abs(as.matrix(assmc))
  amg[is.na(amg)] = 0
  diag(amg) = 0
  amg = amg[,colSums(amg!=0)>0]
  amg = amg[rowSums(amg!=0)>0,]
  ag = igraph::graph_from_adjacency_matrix(adjmatrix = amg, 
                                           mode = "undirected", weighted = TRUE, diag = FALSE)
  nod.name = V(ag)$name
  centr.deg = igraph::centr_degree(ag, mode = "all", loops = FALSE, 
                                   normalized = TRUE)
  nod.deg = centr.deg$res
  net.deg = mean(nod.deg)
  nod.wdeg = igraph::strength(ag, mode = "all", loops = FALSE)
  net.wdeg = mean(nod.wdeg)
  require(sna)
  nod.transit = igraph::transitivity(ag, type = "local")
  net.noden = nrow(assmc)
  net.edgen = sum(abs(amg[upper.tri(amg)]) > 0)
  net.meandis = igraph::mean_distance(ag, directed = FALSE)
  net.density = igraph::edge_density(ag, loops = FALSE)
  net.cc = mean(nod.transit, na.rm = TRUE)
  net.transit = igraph::transitivity(ag, type = "global")
  net.connect = sna::connectedness(amg)
  fitpowerlaw = function(graph) {
    d = igraph::degree(graph, mode = "all")
    dd = igraph::degree_distribution(graph, mode = "all", 
                                     loops = FALSE, cumulative = FALSE)
    degree = 1:max(d)
    probability = dd[-1]
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
    lmijCS=Anova(reg,type = "II")
    p = lmijCS$`Pr(>F)`[[1]]
    list(alpha = alpha, R.square = R.square,p = p)
  }
  net.plr2 = fitpowerlaw(ag)$R.square
  net.alpha = fitpowerlaw(ag)$alpha
  net.p = fitpowerlaw(ag)$p
  net.att = data.frame(NetworkIndex = c("Total nodes", "Total links", "Average degree (avgK)", "Average weighted degree (avgKw)", 
                                        "Average clustering coefficient (avgCC)", "Average path distance (GD)", 
                                        "Density (D)","Transitivity (Trans)", "Krackhardt Connectedness (Con)",
                                        "alpha","r2","P"), 
                       Value = c(net.noden, net.edgen,net.deg, net.wdeg,  
                                 net.cc, net.meandis, net.density,net.transit, net.connect,
                                 net.alpha,net.plr2,net.p), 
                       stringsAsFactors = FALSE)
  net.att
}
net.rand = function (assmc)
{
  require(igraph)
  amg = abs(as.matrix(assmc))
  amg[is.na(amg)] = 0
  diag(amg) = 0
  amg = amg[,colSums(amg!=0)>0]
  amg = amg[rowSums(amg!=0)>0,]
  ag = igraph::graph_from_adjacency_matrix(adjmatrix = amg, 
                                           mode = "undirected", weighted = TRUE, diag = FALSE)
  neta.obs = net.pro(assmc)
  modul.obs = module(assmc = assmc, absolute = TRUE, methods = c("greedy"))
  mod.obs = modul.obs$module
  mod.node = modul.obs$node
  
  netp.obs = rbind(neta.obs, data.frame(NetworkIndex = c(paste0("Module.number.", mod.obs$Method), paste0("Modularity.", mod.obs$Method)), 
                                        Value = c(mod.obs$Module.number, mod.obs$Modularity), 
                                        stringsAsFactors = FALSE))
  rand.method = c("swap")
  randone <- function(ag, rand.method, swap.iter = 100, 
                      endpoint.p = 0.5) {
    require(igraph)
    if (rand.method == "swap") {
      randm = igraph::rewire(ag, with = igraph::keeping_degseq(loops = FALSE, 
                                                               niter = igraph::vcount(ag) * 10))
      randm = igraph::set_edge_attr(graph = randm, 
                                    name = "weight", value = E(ag)$weight)
    }
    else if (rand.method == "endpoints") {
      randm = igraph::rewire(ag, with = igraph::each_edge(prob = endpoint.p, 
                                                          loops = FALSE, multiple = FALSE))
    }
    else {
      stop("rand.method is not valid.")
    }
    igraph::as_adj(randm, sparse = FALSE, names = TRUE, 
                   attr = "weight")
  }
  netp.rand <- sapply(1:100, function(i) {
    message("Now i=", i, " in ", 100, ". ", date())
    ramgi = randone(ag = ag, rand.method = rand.method, 
                    swap.iter = swap.iter, endpoint.p = endpoint.p)
    netai = net.pro(ramgi)
    moduli = module(assmc = ramgi, absolute = TRUE,methods = c("greedy"))
    modi = moduli$module
    as.numeric(c(netai[,2], modi$Module.number, 
                 modi$Modularity))
  })
  z.test = function(a, mu, two.tail = FALSE) {
    zeta = (mean(a, na.rm = TRUE) - mu)/sd(a, na.rm = TRUE)
    if (two.tail) {
      p = 2 * pnorm(-abs(zeta))
    }
    else {
      p = pnorm(-abs(zeta))
    }
    list(z = zeta, p = p)
  }
  random.mean = rowMeans(netp.rand, na.rm = TRUE)
  random.sd = apply(netp.rand, 1, sd, na.rm = TRUE)
  z.value = (random.mean - as.numeric(netp.obs$Value))/random.sd
  p.ztest = pnorm(-abs(z.value))
  EPS <- sqrt(.Machine$double.eps)
  p.lower = (rowSums(netp.rand >= (matrix(as.numeric(netp.obs$Value), 
                                          nrow = nrow(netp.rand), ncol = ncol(netp.rand)) - EPS), 
                     na.rm = TRUE) + 1)/(1 + rowSums(!is.na(netp.rand)))
  p.higher = (rowSums(netp.rand <= (matrix(as.numeric(netp.obs$Value), 
                                           nrow = nrow(netp.rand), ncol = ncol(netp.rand)) + EPS), 
                      na.rm = TRUE) + 1)/(1 + rowSums(!is.na(netp.rand)))
  p.count = apply(cbind(p.lower, p.higher), 1, min)
  out = data.frame(netp.obs, random.mean = random.mean, random.sd = random.sd, 
                   z.value = z.value, p.ztest = p.ztest, p.count = p.count,stringsAsFactors = FALSE)
  colnames(out) = c("NetworkIndex", "Empirical.Results", "Random.Mean", 
                    "Random.Stdev", "Z.value", "P.Ztest", "P.count")
  out
}
library(car)
library(MENA)

prefix1 <- "Fungi"
cor.meth <- "pearson"
maj.1 <- 0.25
prefixm <- cor.meth

random.table = data.frame(stringsAsFactors = FALSE)
for (i in c("MS", "TS", "DS")) {
  for (j in c("L1", "L2", "L3", "L4")) {
    # dat <- read.csv(paste0(i, ".", j, ".", prefixm, ".", 
    #                        prefix1, ".m", maj.1, ".CLR.FillPairLB.csv"), row.names = 1,
    #                 header = TRUE, sep = ",")
    load(paste0(i, "_", j, "_", prefixm, "_", prefix1,"_", maj.1, "_CLR_FillPairLB.rda"))
    dat = output1[["assmc"]]
  result = net.rand(dat)
  result$steppes = i
  result$layer = j
  random.table = rbind(random.table, result)
  }
}
write.csv(random.table,paste0(prefix1, ".",cor.meth,".", maj.1,".","network random compare.csv"))
