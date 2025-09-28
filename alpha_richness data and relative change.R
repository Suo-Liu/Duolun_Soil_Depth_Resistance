setwd("C:/Users/True/OneDrive/桌面")
treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

Bac.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\alpha\\diversity index\\alpha index_Bacteria.csv"
Fungi.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\alpha\\diversity index\\alpha index_Fungi.csv"
Pro.divindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\alpha\\diversity index\\alpha index_Protist.csv"

Bac.divindex <- read.csv(Bac.divindex.file, row.names = 1, header = T, sep = ",")
Fungi.divindex <- read.csv(Fungi.divindex.file, row.names = 1, header = T, sep = ",")
Pro.divindex <- read.csv(Pro.divindex.file, row.names = 1, header = T, sep = ",")

Bac.divindex <- Bac.divindex[, 1, drop = F]
Fungi.divindex <- Fungi.divindex[, 1, drop = F]
Pro.divindex <- Pro.divindex[, 1, drop = F]

divindex.data <- cbind(Bac.divindex, Fungi.divindex, Pro.divindex)
colnames(divindex.data) <- c("Bacteria", "Fungi", "Protist")
treat <- read.csv(treat.file, header = T, row.names = 1)
treat <- subset(treat, plant.type == "TS")

library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, divindex.data = divindex.data))
divindex.data <- sampc$divindex.data
treat <- sampc$treat

# richness top vs sub ####
treat$Layer <- ifelse(treat$Layer %in% c("L1", "L2"), "Topsoil", "Subsoil")

divindex.data$Layer <- treat$Layer
divindex.data$block <- treat$block
divindex.data$treat <- treat$combined_treat1

plot.data <- divindex.data
plot.data$ID <- rownames(plot.data)
library(reshape2)
library(tidyverse)
plot.data <- melt(plot.data, id = c("ID", "Layer", "block", "treat"))
## plot ####
library(ggplot2)
library(gghalves)
plot.data$Layer <- factor(plot.data$Layer, levels = c("Topsoil", "Subsoil"))
ggplot(
  data = plot.data,
  aes(x = Layer, y = value, fill = Layer)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.5) +
  geom_half_boxplot(
    side = "r", errorbar.draw = FALSE, width = 0.25,
    # position = position_dodge(width = 0.76),
    linewidth = 0.2
  ) +
  geom_half_point_panel(
    side = "l", shape = 21, size = 2,
    alpha = 0.7
    # position = position_dodge(width = 0.78)
    , color = "white"
  ) +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("Topsoil", "Subsoil")) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  labs(y = "Richness", x = NULL) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.position = "none",
    strip.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
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
ggsave("richness_top_sub.pdf", width = 7.17, height = 2.40, units = "in")

## statistics test ####
div.table <- data.frame(stringsAsFactors = F)

library(tidyverse)
library(magrittr)
plot.data %<>% pivot_wider(
  names_from = variable,
  values_from = value
)

divindex <- plot.data[, 5:7, drop = F]
treat.used <- plot.data[, -c(5:7)]
treat.used$layer <- ifelse(treat.used$Layer == "Topsoil", 0, 1)
divindex <- scale(divindex)
library(car)
library(lme4)
divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ layer + (1 | block) + (1 | treat), data = div)
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
  fm <- lmer(divtest ~ layer + (1 | block), data = div)
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
write.csv(div.table, paste0("richness", "_", "data", "_", "P value.csv"))

# env relative change top vs sub ####
## relative change calculation ####
pre.divindex.data <- divindex.data
result.table <- data.frame(stringsAsFactors = F)
for (k in c("Bacteria", "Fungi", "Protist")) {
  divindex.data <- pre.divindex.data[, k, drop = F]
  result <- data.frame(
    Divisor = character(), Dividend = character(),
    Difference = numeric(), stringsAsFactors = FALSE
  )

  for (i in 1:(nrow(divindex.data) - 1)) {
    for (j in (i + 1):nrow(divindex.data)) {
      result <- rbind(result, data.frame(
        Divisor = rownames(divindex.data)[i],
        Dividend = rownames(divindex.data)[j],
        Dif_data = divindex.data[i, 1] - divindex.data[j, 1]
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

      result$Dif_data[i] <- -result$Dif_data[i]
    } else if (result$treat1[i] == "W" & result$treat2[i] == "WRP") {
      a <- result$Dividend[i]
      result$Dividend[i] <- result$Divisor[i]
      result$Divisor[i] <- a

      b <- result$treat1[i]
      result$treat1[i] <- result$treat2[i]
      result$treat2[i] <- b

      result$Dif_data[i] <- -result$Dif_data[i]
    } else if (result$treat1[i] == "W" & result$treat2[i] == "WEP") {
      a <- result$Dividend[i]
      result$Dividend[i] <- result$Divisor[i]
      result$Divisor[i] <- a

      b <- result$treat1[i]
      result$treat1[i] <- result$treat2[i]
      result$treat2[i] <- b

      result$Dif_data[i] <- -result$Dif_data[i]
    } else if (result$treat1[i] == "EP" & result$treat2[i] == "WEP") {
      a <- result$Dividend[i]
      result$Dividend[i] <- result$Divisor[i]
      result$Divisor[i] <- a

      b <- result$treat1[i]
      result$treat1[i] <- result$treat2[i]
      result$treat2[i] <- b

      result$Dif_data[i] <- -result$Dif_data[i]
    } else if (result$treat1[i] == "RP" & result$treat2[i] == "WRP") {
      a <- result$Dividend[i]
      result$Dividend[i] <- result$Divisor[i]
      result$Divisor[i] <- a

      b <- result$treat1[i]
      result$treat1[i] <- result$treat2[i]
      result$treat2[i] <- b

      result$Dif_data[i] <- -result$Dif_data[i]
    }
  }

  result$data <- divindex.data[match(result$Dividend, rownames(divindex.data)), 1]
  result$index <- k
  result.table <- rbind(result.table, result)
}

result.table$relative_change <- abs(result.table$Dif_data) / result.table$data

result.table <- result.table[, c("relative_change", "index", "block1", "treat3", "layer1")]

write.csv(result.table, paste0("richness", "_", "relative_change", ".csv"))

## plot ####
library(ggplot2)
library(tidyverse)
library(gghalves)

data.file <- paste0("C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\alpha\\relative change\\richness_relative_change.csv")

plot.data <- read.csv(data.file, row.names = 1, header = T, sep = ",")

plot.data$layer1 <- ifelse(plot.data$layer1 %in% c("L1", "L2"), "Topsoil", "Subsoil")
plot.data$layer1 <- factor(plot.data$layer1, levels = c("Topsoil", "Subsoil"))

ggplot(
  data = plot.data,
  aes(x = layer1, y = relative_change, fill = layer1)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(
    side = "r", errorbar.draw = FALSE,
    width = 0.2, linewidth = 0.5,
    outlier.size = 1
  ) +
  geom_half_point_panel(
    side = "l", shape = 21,
    size = 1.5, color = "white"
  ) +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("Topsoil", "Subsoil")) +
  labs(y = "Relative change", x = NULL) +
  facet_wrap(~index, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.position = "none",
    strip.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
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

ggsave("richness_relative change_top_sub.pdf", width = 7.17, height = 2.40, units = "in")

## statistics test ####
div.table <- data.frame(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
plot.data <- plot.data[order(plot.data$index), ]
plot.data <- plot.data[order(plot.data$block1), ]
plot.data <- plot.data[order(plot.data$treat3), ]
plot.data <- plot.data[order(plot.data$layer1), ]
plot.data$ID <- rep(1:2, nrow(plot.data) / 2)
plot.data %<>% pivot_wider(
  names_from = index,
  values_from = relative_change
)

divindex <- plot.data[, 5:7, drop = F]
treat.used <- plot.data[, -c(5:7)]
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
write.csv(div.table, paste0("richness", "_", "relative change", "_", "P value.csv"))
