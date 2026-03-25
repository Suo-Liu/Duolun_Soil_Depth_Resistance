setwd("C:/Users/True/OneDrive/桌面")
# network properties resistance under treatments among different steppes
# resistance calculation ####
cor.meth <- "pearson"
treat.file <- "treatment.csv"

prefix1 <- "Bac"
netindex.file <- "Bacteria_network_property.csv"

prefix <- paste0(prefix1, "_", cor.meth)

netindex <- read.csv(netindex.file, row.names = 1, header = T, sep = ",")
scale_cols <- c(
  "Average.clustering.coefficient", "Average.degree",
  "Density", "Krackhardt.Connectedness", "Total.links"
)

min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

netindex[, scale_cols] <- apply(netindex[, scale_cols], 2, min_max_scale)
netindex$complexity <- rowSums(netindex[, scale_cols], na.rm = TRUE)

netindex <- netindex[, "complexity", drop = F]
netindex <- netindex / length(scale_cols)

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")

netindex <- netindex[match(rownames(treat), rownames(netindex)), , drop = F]

result <- data.frame(Divisor = character(), Dividend = character(), Difference = numeric(), stringsAsFactors = FALSE)

for (i in 1:(nrow(netindex) - 1)) {
  for (j in (i + 1):nrow(netindex)) {
    result <- rbind(result, data.frame(
      Divisor = rownames(netindex)[i],
      Dividend = rownames(netindex)[j],
      Dif_index = netindex[i, 1] - netindex[j, 1]
    ))
  }
}

result$layer1 <- treat$Layer[match(result$Divisor, rownames(treat))]
result$layer2 <- treat$Layer[match(result$Dividend, rownames(treat))]

result$block1 <- treat$block[match(result$Divisor, rownames(treat))]
result$block2 <- treat$block[match(result$Dividend, rownames(treat))]

result$treat1 <- treat$combined_treat1[match(result$Divisor, rownames(treat))]
result$treat2 <- treat$combined_treat1[match(result$Dividend, rownames(treat))]

result <- result[result$block1 == result$block2 &
  result$layer1 == result$layer2, ]

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

    result$Dif_index[i] <- -result$Dif_index[i]
  } else if (result$treat1[i] == "W" & result$treat2[i] == "WRP") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_index[i] <- -result$Dif_index[i]
  } else if (result$treat1[i] == "W" & result$treat2[i] == "WEP") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_index[i] <- -result$Dif_index[i]
  } else if (result$treat1[i] == "EP" & result$treat2[i] == "WEP") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_index[i] <- -result$Dif_index[i]
  } else if (result$treat1[i] == "RP" & result$treat2[i] == "WRP") {
    a <- result$Dividend[i]
    result$Dividend[i] <- result$Divisor[i]
    result$Divisor[i] <- a

    b <- result$treat1[i]
    result$treat1[i] <- result$treat2[i]
    result$treat2[i] <- b

    result$Dif_index[i] <- -result$Dif_index[i]
  }
}

result$complexity <- netindex[match(result$Dividend, rownames(netindex)), 1]

result$resistance <- 1 - (2 * abs(result$Dif_index)) / (abs(result$Dif_index) + result$complexity)

result <- result[, c("resistance", "block1", "treat3", "layer1")]

write.csv(result, paste0("resistance", "_", prefix, ".csv"))

# plot ####
library(ggplot2)
library(tidyverse)
library(gghalves)

prefix.m <- "Bac"
cor.meth <- "pearson"
data.file <- paste0("resistance_", prefix.m, "_", cor.meth, ".csv")
result <- read.csv(data.file, row.names = 1, header = T, sep = ",")
result$layer1 <- ifelse(result$layer1 %in% c("L1", "L2"), "Topsoil", "Subsoil")
result$layer1 <- factor(result$layer1, levels = c("Topsoil", "Subsoil"))

plot.data <- result

ggplot(
  data = plot.data,
  aes(x = layer1, y = resistance, fill = layer1)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 2, color = "white") +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("Topsoil", "Subsoil")) +
  labs(y = "Resistance", x = NULL) +
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

ggsave(paste0("all", "_", prefix.m, "_", "all treatments.paired treatment.pdf"),
  width = 2.97, height = 2.81, units = "in"
)

# statistics test ####
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
