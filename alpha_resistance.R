setwd("C:/Users/True/OneDrive/桌面")
# alpha diversity resistance under treatments among different steppes
# resistance calculation ####
# adopt the resistance index
# reference:  Orwin, K. H. & Wardle, D. A. New indices for quantifying the resistance and resilience of soil biota to exogenous disturbances. Soil Biol. Biochem. 36, 1907–1912 (2004)
# Peng, Z., van der Heijden, M.G.A., Liu, Y. et al. Agricultural subsoil microbiomes and functions exhibit lower resistance to global change than topsoils in Chinese agroecosystems. Nat FoodIF 23.6SCIEJCR (2025)
# https://mp.weixin.qq.com/s/XZ3PDK6vetVO2y5CmNhhAQ

multi <- F
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\data for use\\treatment.csv"

if (!multi) {
  prefix.m <- "Bacteria"
  divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\alpha\\diversity index\\alpha index_Bacteria.csv"

  prefix.m <- "Fungi"
  divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\alpha\\diversity index\\alpha index_Fungi.csv"

  prefix.m <- "Protist"
  divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\alpha\\diversity index\\alpha index_Protist.csv"

  prefix.m <- "Nematode"
  divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\alpha\\diversity index\\alpha index_Nematode.csv"

  divindex.data <- read.csv(divindex.file, row.names = 1, header = T, sep = ",")
  divindex.data <- divindex.data[, 1, drop = F]
} else {
  prefix.m <- "multidiversity"

  bac.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\alpha\\diversity index\\alpha index_Bacteria.csv"

  fungi.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\alpha\\diversity index\\alpha index_Fungi.csv"

  pro.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\alpha\\diversity index\\alpha index_Protist.csv"

  nem.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 transplant Duolun\\alpha\\diversity index\\alpha index_Nematode.csv"

  bac.divindex.data <- read.csv(bac.divindex.file, row.names = 1, header = T, sep = ",")
  fungi.divindex.data <- read.csv(fungi.divindex.file, row.names = 1, header = T, sep = ",")
  pro.divindex.data <- read.csv(pro.divindex.file, row.names = 1, header = T, sep = ",")
  nem.divindex.data <- read.csv(nem.divindex.file, row.names = 1, header = T, sep = ",")

  bac.divindex.data <- bac.divindex.data[, 1, drop = F]
  fungi.divindex.data <- fungi.divindex.data[, 1, drop = F]
  pro.divindex.data <- pro.divindex.data[, 1, drop = F]
  nem.divindex.data <- nem.divindex.data[, 1, drop = F]

  bac.divindex.data <- apply(bac.divindex.data, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  fungi.divindex.data <- apply(fungi.divindex.data, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  pro.divindex.data <- apply(pro.divindex.data, 2, function(x) (x - min(x)) / (max(x) - min(x)))
  nem.divindex.data <- apply(nem.divindex.data, 2, function(x) (x - min(x)) / (max(x) - min(x)))

  divindex.data <- bac.divindex.data + fungi.divindex.data + pro.divindex.data + nem.divindex.data
}

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")

library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, divindex.data = divindex.data))
divindex.data <- sampc$divindex.data
treat <- sampc$treat

divindex.data <- divindex.data[, "richness", drop = F]
result <- data.frame(Divisor = character(), Dividend = character(), Difference = numeric(), stringsAsFactors = FALSE)

for (i in 1:(nrow(divindex.data) - 1)) {
  for (j in (i + 1):nrow(divindex.data)) {
    result <- rbind(result, data.frame(
      Divisor = rownames(divindex.data)[i],
      Dividend = rownames(divindex.data)[j],
      Dif_richness = divindex.data[i, 1] - divindex.data[j, 1]
    ))
  }
}

result$layer1 <- treat$Layer[match(result$Divisor, rownames(treat))]
result$layer2 <- treat$Layer[match(result$Dividend, rownames(treat))]

result$plant1 <- treat$plant.type[match(result$Divisor, rownames(treat))]
result$plant2 <- treat$plant.type[match(result$Dividend, rownames(treat))]

result$block1 <- treat$block[match(result$Divisor, rownames(treat))]
result$block2 <- treat$block[match(result$Dividend, rownames(treat))]

result$treat1 <- treat$combined_treat1[match(result$Divisor, rownames(treat))]
result$treat2 <- treat$combined_treat1[match(result$Dividend, rownames(treat))]

result <- result[result$block1 == result$block2 &
  result$layer1 == result$layer2 &
  result$plant1 == result$plant2, ]

result <- result[
  (result$treat1 %in% c("EP", "C") & result$treat2 %in% c("EP", "C")) |
    (result$treat1 %in% c("C", "W") & result$treat2 %in% c("C", "W")) |
    (result$treat1 %in% c("RP", "C") & result$treat2 %in% c("RP", "C")) |
    (result$treat1 %in% c("WEP", "C") & result$treat2 %in% c("WEP", "C")) |
    (result$treat1 %in% c("WRP", "C") & result$treat2 %in% c("WRP", "C")) |
    (result$treat1 %in% c("WEP", "W") & result$treat2 %in% c("WEP", "W")) |
    (result$treat1 %in% c("WRP", "W") & result$treat2 %in% c("WRP", "W")) |
    (result$treat1 %in% c("WRP", "RP") & result$treat2 %in% c("WRP", "RP")) |
    (result$treat1 %in% c("WEP", "EP") & result$treat2 %in% c("WEP", "EP")),
]

result$treat3 <- ifelse(result$treat1 == "C", result$treat2, result$treat1)
result$treat3 <- ifelse(result$treat1 %in% c("W", "WRP") & result$treat2 %in% c("W", "WRP"),
  "WRP.W", ifelse(result$treat1 %in% c("W", "WEP") & result$treat2 %in% c("W", "WEP"), "WEP.W",
    ifelse(result$treat1 %in% c("WRP", "RP") & result$treat2 %in% c("WRP", "RP"), "WRP.RP",
      ifelse(result$treat1 %in% c("WEP", "EP") & result$treat2 %in% c("WEP", "EP"), "WEP.EP", result$treat3)
    )
  )
)

for (i in 1:nrow(result)) {
  if (result$treat1[i] == "C") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_richness[i] <- -result$Dif_richness[i]
  } else if (result$treat1[i] == "W" & result$treat2[i] == "WRP") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_richness[i] <- -result$Dif_richness[i]
  } else if (result$treat1[i] == "W" & result$treat2[i] == "WEP") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_richness[i] <- -result$Dif_richness[i]
  } else if (result$treat1[i] == "EP" & result$treat2[i] == "WEP") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_richness[i] <- -result$Dif_richness[i]
  } else if (result$treat1[i] == "RP" & result$treat2[i] == "WRP") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_richness[i] <- -result$Dif_richness[i]
  }
}

result$richness <- divindex.data[match(result$Dividend, rownames(divindex.data)), 1]

result$resistance <- 1 - (2 * abs(result$Dif_richness)) / (abs(result$Dif_richness) + result$richness)

result <- result[, c("resistance", "block1", "treat3", "plant1", "layer1")]

write.csv(result, paste0("alpha", "_", "resistance", "_", prefix.m, ".csv"))

# plot ####
library(ggplot2)
library(tidyverse)
library(gghalves)

prefix.m <- "Protist"
data.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\alpha\\resistance\\", "alpha", "_", "resistance", "_", prefix.m, ".csv")
result <- read.csv(data.file, row.names = 1, header = T, sep = ",")
result <- subset(result, plant1 == "TS")

# soil layers ####
## plot all ####
plot.data <- result
# plot.data <- subset(result, treat3 %in% c("RP", "EP", "W", "WRP", "WEP"))

plot.data$layer1 <- ifelse(plot.data$layer1 %in% c("L1", "L2"), "Topsoil", "Subsoil")
plot.data$layer1 <- factor(plot.data$layer1, levels = c("Topsoil", "Subsoil"))

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

# ggsave(paste0("all","_",prefix.m, "_", "all treatments.paired treatment.pdf"), width = 2.97, height = 2.81, units = "in")
ggsave(paste0("all", "_", prefix.m, "_", "unique treatments.paired treatment.pdf"), width = 2.97, height = 2.81, units = "in")

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
plot.data <- result
# plot.data <- subset(result, treat3 %in% c("RP", "EP", "W", "WRP", "WEP"))

plot.data$layer1 <- ifelse(plot.data$layer1 %in% c("L1", "L2"), "Topsoil", "Subsoil")
plot.data$layer1 <- factor(plot.data$layer1, levels = c("Topsoil", "Subsoil"))

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

ggsave(paste0("per treatments", "_", prefix.m, ".pdf"), width = 6.26, height = 5, units = "in")

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
write.csv(result.table, paste0(prefix.m, "_", "per treatment", "_", "t&wilcox", "_", "P value.csv"))
