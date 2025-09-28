setwd("C:/Users/True/OneDrive/桌面")

prefix1 <- "Bac"
prefix2 <- "Pro"
cor.meth <- "pearson"
maj.1 <- 0.75
maj.2 <- 0.25

treat.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance transplant Duolun\\data for use\\treatment.csv"

prefix <- paste0(cor.meth, "_", prefix1, "_", maj.1, "_", prefix2, "_", maj.2)
netindex.file <- "C:\\Users\\True\\OneDrive\\桌面\\research\\03252024 resistance depth Duolun\\network\\BP\\Bac_0.75_Pro_0.25_pearson.csv"

netindex <- read.csv(netindex.file, row.names = 1, header = T, sep = ",")
scale_cols <- c(
  "linkage.density", "Average.degree.all",
  "Density", "Connectance", "Total.links",
  "weighted.connectance"
)

netindex = netindex[,scale_cols]

treat <- read.csv(treat.file, row.names = 1, header = T, sep = ",")
treat = treat[, c("Layer", "block", "combined_treat1")]
library(ieggr)
sampc <- match.name(rn.list = list(treat = treat, netindex = netindex), silent = T)
netindex <- sampc$netindex
treat <- sampc$treat

result <- cbind(netindex, treat)
result$Layer <- ifelse(result$Layer %in% c("L1", "L2"),
                       "Topsoil", "Subsoil")

## plot ####
library(ggplot2)
library(tidyverse)
library(gghalves)
library(reshape2)

plot.data <- result
plot.data = melt(plot.data, id.vars = c("Layer", "block", "combined_treat1"))
plot.data$Layer <- factor(plot.data$Layer, levels = c("Topsoil", "Subsoil"))
plot.data = subset(plot.data, variable %in% c("linkage.density", "Average.degree.all",
                                              "Total.links"))
plot.data$variable = ifelse(plot.data$variable == "linkage.density", "Linkage density",
                            ifelse(plot.data$variable == "Average.degree.all", "Average degree",
                                   ifelse(plot.data$variable == "Total.links", "Link number", plot.data$variable)))
plot.data$variable <- factor(plot.data$variable, levels = c("Link number", "Average degree",
                                                            "Linkage density"))

ggplot(
  data = plot.data,
  aes(x = Layer, y = value, fill = Layer)
) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 2, color = "white") +
  scale_fill_manual(values = c("#CD2626", "#1874CD"), limits = c("Topsoil", "Subsoil")) +
  labs(y = "Value", x = NULL) +
  facet_wrap(~variable, scales = "free_y", ncol = 3)+
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
      hjust = 0.5,
      colour = "#000000"
    ),
    axis.text = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      colour = "#000000"
    )
  )

ggsave(paste0("complexity", "_", prefix, ".pdf"),
       width = 6.99, height = 2.22, units = "in"
)

## simple statistics test ####
div.table <- data.frame(stringsAsFactors = F)

divindex <- result[, 1:6, drop = F]
treat.used <- result[, -c(1:6)]
treat.used$Layer <- ifelse(treat.used$Layer == "Topsoil", 0, 1)
divindex <- scale(divindex)
library(car)
library(lme4)
divs1 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  div <- data.frame(divtest = divindex[, j], treat.used)
  fm <- lmer(divtest ~ Layer + (1 | block) + (1 | combined_treat1), data = div)
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
  fm <- lmer(divtest ~ Layer + (1 | block), data = div)
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
write.csv(div.table, paste0("complexity","_",prefix, "_", "P value.csv"))
