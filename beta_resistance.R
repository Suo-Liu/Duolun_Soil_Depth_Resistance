setwd("C:/Users/True/OneDrive/桌面")

# calculation ####
library(ieggr)
library(dplyr)
library(vegan)
library(Rmisc)
library(reshape2)

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

## single community for resistance ####
com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"
prefix.m <- "Bacteria"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\ITS\\all samples\\Unoise\\zotutab220_max_2_resample_10000.txt"
prefix.m <- "Fungi"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\18S\\all samples\\Unoise\\protist_zotu_2868.txt"
prefix.m <- "Protist"

com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\18S\\all samples\\Unoise\\nematoda_zotu_259.txt"
prefix.m <- "Nematode"

comm <- t(read.table(com.file,
                     header = TRUE, sep = "\t", row.names = 1,
                     as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE
))
name.row <- rownames(comm)
name.row[which(rownames(comm) == "YD3L3")] <- "YD52L3"
name.row[which(rownames(comm) == "YD52L3")] <- "YD3L3"
rownames(comm) <- name.row

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")

treat$plot <- as.character(treat$plot)
treat$Layer <- as.character(treat$Layer)
treat$block <- as.character(treat$block)
treat$plant.type <- as.character(treat$plant.type)
treat$combined_treat1 <- as.character(treat$combined_treat1)

sampc <- match.name(rn.list = list(treat = treat, comm = comm))
dim(comm)
comm <- sampc$comm
treat <- sampc$treat
comm <- comm[, colSums(comm) > 0]
dim(comm)

## multiple community for resistance ####
prefix.m <- "mixed"

bac.com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\16S\\all samples\\Unoise\\bacteria.zotu245_resample_35000.txt"

fungi.com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\ITS\\all samples\\Unoise\\zotutab220_max_2_resample_10000.txt"

pro.com.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\work\\inner mongolia\\data use for analysis\\18S\\all samples\\Unoise\\protist.zotu_resample_2961.txt"


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

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")

sampc <- match.name(rn.list = list(treat = treat, bac.comm = bac.comm, fungi.comm = fungi.comm, pro.comm = pro.comm))
bac.comm <- sampc$bac.comm
fungi.comm <- sampc$fungi.comm
pro.comm <- sampc$pro.comm
treat <- sampc$treat

bac.comm <- bac.comm[, colSums(bac.comm) > 0]
fungi.comm <- fungi.comm[, colSums(fungi.comm) > 0]
pro.comm <- pro.comm[, colSums(pro.comm) > 0]

comm <- cbind(bac.comm, fungi.comm, pro.comm)
dim(comm)

dist.used <- vegdist(comm)
table.used1 <- dist.3col(dist.used)

table.used1$plot1 <- treat$plot[match(table.used1$name1, rownames(treat))]
table.used1$plot2 <- treat$plot[match(table.used1$name2, rownames(treat))]
table.used1$layer1 <- treat$Layer[match(table.used1$name1, rownames(treat))]
table.used1$layer2 <- treat$Layer[match(table.used1$name2, rownames(treat))]
table.used1$treat1 <- treat$combined_treat1[match(table.used1$name1, rownames(treat))]
table.used1$treat2 <- treat$combined_treat1[match(table.used1$name2, rownames(treat))]
table.used1$block1 <- treat$block[match(table.used1$name1, rownames(treat))]
table.used1$block2 <- treat$block[match(table.used1$name2, rownames(treat))]
table.used1$plant1 <- treat$plant.type[match(table.used1$name1, rownames(treat))]
table.used1$plant2 <- treat$plant.type[match(table.used1$name2, rownames(treat))]

table.used1 <- table.used1[!(table.used1$treat1 == table.used1$treat2), ]
table.used1 <- table.used1[table.used1$block1 == table.used1$block2 &
  table.used1$layer1 == table.used1$layer2 &
  table.used1$plant1 == table.used1$plant2, ]

table.used1 <- table.used1[
  (table.used1$treat1 %in% c("EP", "C") & table.used1$treat2 %in% c("EP", "C")) |
    (table.used1$treat1 %in% c("WEP", "C") & table.used1$treat2 %in% c("WEP", "C")) |
    (table.used1$treat1 %in% c("WRP", "C") & table.used1$treat2 %in% c("WRP", "C")) |
    (table.used1$treat1 %in% c("W", "C") & table.used1$treat2 %in% c("W", "C")) |
    (table.used1$treat1 %in% c("RP", "C") & table.used1$treat2 %in% c("RP", "C")) |
    (table.used1$treat1 %in% c("WEP", "W") & table.used1$treat2 %in% c("WEP", "W")) |
    (table.used1$treat1 %in% c("WRP", "W") & table.used1$treat2 %in% c("WRP", "W")) |
    (table.used1$treat1 %in% c("WRP", "RP") & table.used1$treat2 %in% c("WRP", "RP")) |
    (table.used1$treat1 %in% c("WEP", "EP") & table.used1$treat2 %in% c("WEP", "EP")),
]

table.used1$treat3 <- ifelse(table.used1$treat1 == "C", table.used1$treat2, table.used1$treat1)
table.used1$treat3 <- ifelse(table.used1$treat1 %in% c("W", "WRP") & table.used1$treat2 %in% c("W", "WRP"),
  "WRP.W", ifelse(table.used1$treat1 %in% c("W", "WEP") & table.used1$treat2 %in% c("W", "WEP"), "WEP.W",
    ifelse(table.used1$treat1 %in% c("WRP", "RP") & table.used1$treat2 %in% c("WRP", "RP"), "WRP.RP",
      ifelse(table.used1$treat1 %in% c("WEP", "EP") & table.used1$treat2 %in% c("WEP", "EP"), "WEP.EP",
        table.used1$treat3
      )
    )
  )
)

table.used1$treat <- "Bray"

dist.used <- vegdist(comm, binary = T)
table.used2 <- dist.3col(dist.used)

table.used2$plot1 <- treat$plot[match(table.used2$name1, rownames(treat))]
table.used2$plot2 <- treat$plot[match(table.used2$name2, rownames(treat))]
table.used2$layer1 <- treat$Layer[match(table.used2$name1, rownames(treat))]
table.used2$layer2 <- treat$Layer[match(table.used2$name2, rownames(treat))]
table.used2$treat1 <- treat$combined_treat1[match(table.used2$name1, rownames(treat))]
table.used2$treat2 <- treat$combined_treat1[match(table.used2$name2, rownames(treat))]
table.used2$block1 <- treat$block[match(table.used2$name1, rownames(treat))]
table.used2$block2 <- treat$block[match(table.used2$name2, rownames(treat))]
table.used2$plant1 <- treat$plant.type[match(table.used2$name1, rownames(treat))]
table.used2$plant2 <- treat$plant.type[match(table.used2$name2, rownames(treat))]

table.used2 <- table.used2[!(table.used2$treat1 == table.used2$treat2), ]
table.used2 <- table.used2[table.used2$block1 == table.used2$block2 &
  table.used2$layer1 == table.used2$layer2 &
  table.used2$plant1 == table.used2$plant2, ]

table.used2 <- table.used2[
  (table.used2$treat1 %in% c("EP", "C") & table.used2$treat2 %in% c("EP", "C")) |
    (table.used2$treat1 %in% c("WEP", "C") & table.used2$treat2 %in% c("WEP", "C")) |
    (table.used2$treat1 %in% c("WRP", "C") & table.used2$treat2 %in% c("WRP", "C")) |
    (table.used2$treat1 %in% c("W", "C") & table.used2$treat2 %in% c("W", "C")) |
    (table.used2$treat1 %in% c("RP", "C") & table.used2$treat2 %in% c("RP", "C")) |
    (table.used2$treat1 %in% c("WEP", "W") & table.used2$treat2 %in% c("WEP", "W")) |
    (table.used2$treat1 %in% c("WRP", "W") & table.used2$treat2 %in% c("WRP", "W")) |
    (table.used2$treat1 %in% c("WRP", "RP") & table.used2$treat2 %in% c("WRP", "RP")) |
    (table.used2$treat1 %in% c("WEP", "EP") & table.used2$treat2 %in% c("WEP", "EP")),
]

table.used2$treat3 <- ifelse(table.used2$treat1 == "C", table.used2$treat2, table.used2$treat1)
table.used2$treat3 <- ifelse(table.used2$treat1 %in% c("W", "WRP") & table.used2$treat2 %in% c("W", "WRP"),
  "WRP.W", ifelse(table.used2$treat1 %in% c("W", "WEP") & table.used2$treat2 %in% c("W", "WEP"), "WEP.W",
    ifelse(table.used2$treat1 %in% c("WRP", "RP") & table.used2$treat2 %in% c("WRP", "RP"), "WRP.RP",
      ifelse(table.used2$treat1 %in% c("WEP", "EP") & table.used2$treat2 %in% c("WEP", "EP"), "WEP.EP",
        table.used2$treat3
      )
    )
  )
)

table.used2$treat <- "Sorenson"

table.used <- rbind(table.used1, table.used2)

library(tidyr)

table.used <- table.used[, c("dis", "block1", "treat3", "treat", "plant1", "layer1")]
table.used$dis = 1-table.used$dis
colnames(table.used)[1] <- "resistance"
table.used$group <- paste0(
  table.used$treat3, "_", table.used$treat
)
write.csv(table.used, paste0("beta","_","resistance","_",prefix.m, ".csv"))

# top vs sub ####
library(ggplot2)
library(tidyverse)
library(gghalves)

prefix.m <- "Protist"
data.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\beta\\resistance\\", "beta_", "resistance", "_", prefix.m, ".csv")
result <- read.csv(data.file, row.names = 1, header = T, sep = ",")
result <- subset(result, treat == "Sorenson" & plant1 == "TS")
# result <- subset(result, treat == "Sorenson" & plant1 == "TS" & treat3 %in% c("RP", "EP", "W", "WRP", "WEP"))

prefix.m.1 <- paste0(prefix.m, "_","Sorenson")

result$layer1 <- ifelse(result$layer1 %in% c("L1", "L2"), "Topsoil", "Subsoil")
result$layer1 <- factor(result$layer1, levels = c("Topsoil", "Subsoil"))

## plot all ####
plot.data = result

ggplot(
  data = plot.data,
  aes(x = layer1, y = resistance, fill = layer1)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 2, color = "white") +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("Topsoil", "Subsoil")) +
  # scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  # scale_x_discrete(labels = c('Artio','Carn','Prim','Rod','Others')) +
  labs(y = "Resistance", x = NULL) +
  # facet_wrap(~treat, scales = "free_y", ncol = 3)+
  theme_bw() +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(
      size = 15,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      colour = "#000000"
    )
  )

# ggsave(paste0("all","_",prefix.m.1, "_", "all treatments.paired treatment.pdf"), width = 2.97, height = 2.81, units = "in")
ggsave(paste0("all", "_", prefix.m.1, "_", "unique treatments.paired treatment.pdf"), width = 2.97, height = 2.81, units = "in")

## statistics test ####
div.table <- data.frame(stringsAsFactors = F)

divindex <- plot.data[, 1, drop = F]
treat.used <- plot.data[, -1]
treat.used$layer <- ifelse(treat.used$layer1 == "Topsoil", 0, 1)
divindex <- scale(divindex)
library(car)
library(lme4)
divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ layer + (1 | block1) + (1 | treat3), data = div)
  presult <- car::Anova(fm, type = 2)
  coefs <- coef(summary(fm))[, "Estimate"] # four coefs
  names(coefs) <- paste0(names(coefs), ".mean")
  SEvalues <- coef(summary(fm))[, "Std. Error"] # standard errors
  names(SEvalues) <- paste0(names(SEvalues), ".se")
  tvalues <- coef(summary(fm))[, "t value"] # t values
  names(tvalues) <- paste0(names(tvalues), ".t")
  chisqP <- c(presult[, 1], presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
  chisqP <- c(presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".P"))
  result <- c(coefs, tvalues, SEvalues, chisqP)
})
colnames(divs1) <- colnames(divindex)
divs1 <- divs1[!grepl("Intercept", rownames(divs1)), ]
divs1 <- as.data.frame(divs1)
divs1$treat <- rownames(divs1)
divs1 <- separate(divs1, col = treat, sep = "\\.", into = c("treat", "type"))
divs1 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)
divs1$method <- "with treat"
div.table <- rbind(div.table, divs1)

divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ layer + (1 | block1), data = div)
  presult <- car::Anova(fm, type = 2)
  coefs <- coef(summary(fm))[, "Estimate"] # four coefs
  names(coefs) <- paste0(names(coefs), ".mean")
  SEvalues <- coef(summary(fm))[, "Std. Error"] # standard errors
  names(SEvalues) <- paste0(names(SEvalues), ".se")
  tvalues <- coef(summary(fm))[, "t value"] # t values
  names(tvalues) <- paste0(names(tvalues), ".t")
  chisqP <- c(presult[, 1], presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
  chisqP <- c(presult[, 3])
  names(chisqP) <- c(paste0(row.names(presult), ".P"))
  result <- c(coefs, tvalues, SEvalues, chisqP)
})
colnames(divs1) <- colnames(divindex)
divs1 <- divs1[!grepl("Intercept", rownames(divs1)), ]
divs1 <- as.data.frame(divs1)
divs1$treat <- rownames(divs1)
divs1 <- separate(divs1, col = treat, sep = "\\.", into = c("treat", "type"))
divs1 %<>%
  select(treat, type, everything()) %>%
  arrange(treat, type)
divs1$method <- "without treat"
div.table <- rbind(div.table, divs1)
write.csv(div.table, paste0(prefix.m, "_", "all treatments", "_", "P value.csv"))

## per treatments ####
plot.data = result

ggplot(
  data = plot.data,
  aes(x = layer1, y = resistance, fill = layer1)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 2, color = "white") +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("Topsoil", "Subsoil")) +
  labs(y = "Resistance", x = NULL) +
  facet_wrap(~treat3, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 15),
    axis.title = element_text(
      size = 15,
      face = "bold",
      vjust = 0.5,
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      colour = "#000000"
    )
  )

ggsave(paste0("per treatments", "_", prefix.m.1, ".pdf"), width = 6.26, height = 5, units = "in")

## statistics test ####
table.used <- plot.data
result.table <- data.frame(stringsAsFactors = F)

for (i in unique(table.used$treat3)) {
  temp.data <- subset(table.used, treat3 == i)
  temp.data <- temp.data[order(temp.data$block1), ]
  temp.data <- temp.data[order(temp.data$treat3), ]
  temp.data <- temp.data[order(temp.data$layer1), ]
  
  top.temp.data <- subset(temp.data, layer1 == "Topsoil")
  sub.temp.data <- subset(temp.data, layer1 == "Subsoil")
  
  top.ks <- ks.test(scale(jitter(top.temp.data$resistance, 0.00001)), "pnorm")
  sub.ks <- ks.test(scale(jitter(sub.temp.data$resistance, 0.00001)), "pnorm")
  
  a.mean <- (mean(top.temp.data$resistance) - mean(sub.temp.data$resistance)) / mean(sub.temp.data$resistance)
  
  prefix.pair <- "unpaired"
  a.t <- t.test(
    top.temp.data$resistance,
    sub.temp.data$resistance
  )
  a.wilcox <- wilcox.test(
    top.temp.data$resistance,
    sub.temp.data$resistance
  )
  re.vector <- c(
    top.ks[["p.value"]], sub.ks[["p.value"]],
    a.mean, a.t[["statistic"]][["t"]], a.t[["p.value"]],
    a.wilcox[["statistic"]][["W"]], a.wilcox[["p.value"]],
    i
  )
  result.table <- rbind(result.table, re.vector)
}
colnames(result.table) <- c(
  "top.ks.p.value", "sub.ks.p.value", "richness.mean", "richness.t", "richness.t.P", "richness.W", "richness.W.P",
  "treat"
)
write.csv(result.table, paste0(prefix.m.1, "_", "per treatment", "_", "t&wilcox", "_", "P value.csv"))
