setwd("C:\\Users\\True\\OneDrive\\桌面")
prefix1 <- "Bac"
prefix2 <- "Pro"
col_g <- "#C1C1C1"
cols <- c(
  "#DEB99B", "#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D", "#6E4821", "#A4B423",
  "#C094DF", "#DC95D8", "#326530", "#50C0C9", "#67C021", "#DC69AF", "#8C384F", "#30455C", "#F96C72", "#5ED2BF"
)
for (Layer in c("L1", "L2", "L3", "L4")) {
  load(paste0(
    Layer, "_network_cross domain.rda"
  ))
  occor.r <- output3$assmc
  cor.data.1 <- occor.r
  cor.data.2 <- occor.r

  cor.data.1 <- as.matrix(cor.data.1)
  cor.data.1[is.na(cor.data.1)] <- 0

  cor.data.2 <- as.matrix(cor.data.2)
  cor.data.2[is.na(cor.data.2)] <- 0

  # as the correlation matrix includes the correlations within kingdom, we need to remove them
  cor.data.1 <- cor.data.1[grepl(prefix2, rownames(cor.data.1)), ]
  cor.data.1 <- cor.data.1[, grepl(prefix1, colnames(cor.data.1))]

  # remove isolated node
  cor.data.1 <- cor.data.1[rowSums(abs(cor.data.1)) > 0, ]
  cor.data.1 <- cor.data.1[, colSums(abs(cor.data.1)) > 0]

  cor.data.2[grepl(prefix2, rownames(cor.data.2)), grepl(prefix2, rownames(cor.data.2))] <- 0
  cor.data.2[grepl(prefix1, colnames(cor.data.2)), grepl(prefix1, colnames(cor.data.2))] <- 0

  cor.data.2 <- cor.data.2[rownames(cor.data.2) %in% c(rownames(cor.data.1), colnames(cor.data.1)), ]
  cor.data.2 <- cor.data.2[, colnames(cor.data.2) %in% c(rownames(cor.data.1), colnames(cor.data.1))]

  occor.r <- cor.data.2
  library(igraph)
  g <- graph.adjacency(occor.r, weighted = TRUE, mode = "undirected")
  g <- simplify(g)
  g <- delete.vertices(g, which(degree(g) == 0))

  g1 <- g

  E(g1)$correlation <- E(g1)$weight
  E(g1)$weight <- abs(E(g1)$weight)

  set.seed(007)
  V(g1)$modularity <- membership(cluster_louvain(g1))
  V(g1)$label <- V(g1)$name
  V(g1)$label <- NA

  modu_sort <- V(g1)$modularity %>%
    table() %>%
    sort(decreasing = T)
  top_num <- 16
  modu_name <- names(modu_sort[1:16])

  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name

  V(g1)$color <- V(g1)$modularity
  V(g1)$color[!(V(g1)$color %in% modu_name)] <- col_g
  V(g1)$color[(V(g1)$color %in% modu_name)] <- modu_cols[match(V(g1)$color[(V(g1)$color %in% modu_name)], modu_name)]
  V(g1)$frame.color <- V(g1)$color

  E(g1)$color <- col_g
  for (i in modu_name) {
    col_edge <- cols[which(modu_name == i)]
    otu_same_modu <- V(g1)$name[which(V(g1)$modularity == i)]
    E(g1)$color[(data.frame(as_edgelist(g1))$X1 %in% otu_same_modu) & (data.frame(as_edgelist(g1))$X2 %in% otu_same_modu)] <- col_edge
  }
  V(g1)$color[grepl(prefix1, V(g1)$name)] <- NA

  sub_net_layout <- layout_with_fr(g1, niter = 999, grid = "auto")

  if (Layer == "L1") {
    ga <- g1
    layout1 <- sub_net_layout
  } else if (Layer == "L2") {
    gb <- g1
    layout2 <- sub_net_layout
  } else if (Layer == "L3") {
    gc <- g1
    layout3 <- sub_net_layout
  } else {
    gd <- g1
    layout4 <- sub_net_layout
  }
}

pdf(paste0("network_cross BP", "_", "L1", ".pdf"), width = 5.5, height = 5.5)
plot(ga,
  layout = layout1, edge.color = E(ga)$color,
  vertex.size = 3,
  edge.width = 3
)
title(main = paste0(
  "Nodes=", length(V(ga)$name),
  ", ", "Edges=", nrow(data.frame(as_edgelist(ga)))
))
dev.off()

pdf(paste0("network_cross BP", "_", "L2", ".pdf"), width = 5.5, height = 5.5)
plot(gb,
  layout = layout2, edge.color = E(gb)$color,
  vertex.size = 3,
  edge.width = 3
)
title(main = paste0(
  "Nodes=", length(V(gb)$name),
  ", ", "Edges=", nrow(data.frame(as_edgelist(gb)))
))
dev.off()

pdf(paste0("network_cross BP", "_", "L3", ".pdf"), width = 5.5, height = 5.5)
plot(gc,
  layout = layout3, edge.color = E(gc)$color,
  vertex.size = 3,
  edge.width = 3
)
title(main = paste0(
  "Nodes=", length(V(gc)$name),
  ", ", "Edges=", nrow(data.frame(as_edgelist(gc)))
))
dev.off()

pdf(paste0("network_cross BP", "_", "L4", ".pdf"), width = 5.5, height = 5.5)
plot(gd,
  layout = layout4, edge.color = E(gd)$color,
  vertex.size = 3,
  edge.width = 3
)
title(main = paste0(
  "Nodes=", length(V(gd)$name),
  ", ", "Edges=", nrow(data.frame(as_edgelist(gd)))
))
dev.off()
