setwd("C:/Users/True/OneDrive/桌面")
# beta diversity resistance under treatments among different steppes
# adopt the resistance index
# reference: Peng, Z., van der Heijden, M.G.A., Liu, Y. et al. Agricultural subsoil microbiomes and functions exhibit lower resistance to global change than topsoils in Chinese agroecosystems. Nat FoodIF 23.6SCIEJCR (2025)
# calculation ####
library(dplyr)
library(vegan)
library(Rmisc)
library(reshape2)

treat.file <- "treatment.csv"
prefix.m <- "Bacteria"
com.file <- "bacteria_zotu.txt"

comm <- t(read.table(com.file,
  header = TRUE, sep = "\t", row.names = 1,
  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
  check.names = FALSE
))

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")

comm <- comm[match(rownames(treat), rownames(comm)), ]

dist.used <- vegdist(comm, binary = T)

table.used <- dist.3col(dist.used)

table.used$plot1 <- treat$plot[match(table.used$name1, rownames(treat))]
table.used$plot2 <- treat$plot[match(table.used$name2, rownames(treat))]
table.used$layer1 <- treat$Layer[match(table.used$name1, rownames(treat))]
table.used$layer2 <- treat$Layer[match(table.used$name2, rownames(treat))]
table.used$treat1 <- treat$combined_treat1[match(table.used$name1, rownames(treat))]
table.used$treat2 <- treat$combined_treat1[match(table.used$name2, rownames(treat))]
table.used$block1 <- treat$block[match(table.used$name1, rownames(treat))]
table.used$block2 <- treat$block[match(table.used$name2, rownames(treat))]

table.used <- table.used[!(table.used$treat1 == table.used$treat2), ]
table.used <- table.used[table.used$block1 == table.used$block2 &
  table.used$layer1 == table.used$layer2, ]

table.used <- table.used[
  (table.used$treat1 %in% c("EP", "C") & table.used$treat2 %in% c("EP", "C")) |
    (table.used$treat1 %in% c("WEP", "C") & table.used$treat2 %in% c("WEP", "C")) |
    (table.used$treat1 %in% c("WRP", "C") & table.used$treat2 %in% c("WRP", "C")) |
    (table.used$treat1 %in% c("W", "C") & table.used$treat2 %in% c("W", "C")) |
    (table.used$treat1 %in% c("RP", "C") & table.used$treat2 %in% c("RP", "C")) |
    (table.used$treat1 %in% c("WEP", "W") & table.used$treat2 %in% c("WEP", "W")) |
    (table.used$treat1 %in% c("WRP", "W") & table.used$treat2 %in% c("WRP", "W")) |
    (table.used$treat1 %in% c("WRP", "RP") & table.used$treat2 %in% c("WRP", "RP")) |
    (table.used$treat1 %in% c("WEP", "EP") & table.used$treat2 %in% c("WEP", "EP")),
]

table.used$treat3 <- ifelse(table.used$treat1 == "C", table.used$treat2, table.used$treat1)
table.used$treat3 <- ifelse(table.used$treat1 %in% c("W", "WRP") & table.used$treat2 %in% c("W", "WRP"),
  "WRP.W", ifelse(table.used$treat1 %in% c("W", "WEP") & table.used$treat2 %in% c("W", "WEP"), "WEP.W",
    ifelse(table.used$treat1 %in% c("WRP", "RP") & table.used$treat2 %in% c("WRP", "RP"), "WRP.RP",
      ifelse(table.used$treat1 %in% c("WEP", "EP") & table.used$treat2 %in% c("WEP", "EP"), "WEP.EP",
        table.used$treat3
      )
    )
  )
)

table.used$treat <- "Sorenson"

library(tidyr)

table.used <- table.used[, c("dis", "block1", "treat3", "treat", "layer1")]
table.used$dis <- 1 - table.used$dis
colnames(table.used)[1] <- "resistance"
table.used$group <- paste0(
  table.used$treat3, "_", table.used$treat
)
write.csv(table.used, paste0("beta", "_", "resistance", "_", prefix.m, ".csv"))

# plot ####
library(ggplot2)
library(tidyverse)
library(gghalves)

prefix.m <- "Bacteria"
data.file <- paste0("beta_", "resistance", "_", prefix.m, ".csv")
result <- read.csv(data.file, row.names = 1, header = T, sep = ",")

prefix.m.1 <- paste0(prefix.m, "_", "Sorenson")

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

ggsave(paste0("all", "_", prefix.m.1, "_", "all treatments treatment.pdf"), width = 2.97, height = 2.81, units = "in")

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
