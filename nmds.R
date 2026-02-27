# Packages
#these are the ones Kevin used for nmds, Spearman, and CCA
library(tidyverse)
library(readr)
library(here)
here <- here()
#this is a shortcut to avoid typing the long file path to this folder,
#may not be necessary since it defaults to this data frame
library(lubridate)
library(vegan)
library(beepr)
library(openxlsx)
library(ape)
library(rotl)
library(picante)
library(taxize)
library(ggplot2)
library(ggvegan)
library(corrplot)
library(permute)
library(cluster)
library(plotrix)
library(dplyr)


# sampling data
samp_meta <- read_csv(file = 'data/sample_metadata.csv')
#some NA values from missing DOC data and parameters that don't apply to rainfall samples
#making this data numeric instead of character
samp_meta$DOC_uM <- as.numeric(samp_meta$DOC_uM)
samp_meta$DOC_mgL <- as.numeric(samp_meta$DOC_mgL)
samp_meta$DOCflux_mg <- as.numeric(samp_meta$DOCflux_mg)
samp_meta$DBH_cm <- as.numeric(samp_meta$DBH_cm)
samp_meta$BLDG_dist <- as.numeric(samp_meta$BLDG_dist)
samp_meta$ROAD_dist <- as.numeric(samp_meta$ROAD_dist)
samp_meta$TRAIN_dist <- as.numeric(samp_meta$TRAIN_dist)
samp_meta$DOC_uM <- as.numeric(samp_meta$DOC_uM)
samp_meta$BLDG_sin <- as.numeric(samp_meta$BLDG_sin)
samp_meta$BLDG_cos <- as.numeric(samp_meta$BLDG_cos)
samp_meta$ROAD_sin <- as.numeric(samp_meta$ROAD_sin)
samp_meta$ROAD_cos <- as.numeric(samp_meta$ROAD_cos)
samp_meta$TRAIN_sin <- as.numeric(samp_meta$TRAIN_sin)
samp_meta$TRAIN_cos <- as.numeric(samp_meta$TRAIN_cos)
samp_meta$x <- as.numeric(samp_meta$x)
samp_meta$y <- as.numeric(samp_meta$y)
samp_meta$OL_dist <- as.numeric(samp_meta$OL_dist)
samp_meta$OL_sin <- as.numeric(samp_meta$OL_sin)
samp_meta$OL_cos <- as.numeric(samp_meta$OL_cos)
samp_meta$GL_dist <- as.numeric(samp_meta$GL_dist)
samp_meta$GL_sin <- as.numeric(samp_meta$GL_sin)
samp_meta$GL_cos <- as.numeric(samp_meta$GL_cos)
samp_meta$CR_dist <- as.numeric(samp_meta$CR_dist)
samp_meta$CR_sin <- as.numeric(samp_meta$CR_sin)
samp_meta$CR_cos <- as.numeric(samp_meta$CR_cos)
#this used to be after making tf_meta if anything changes
samp_meta$Tf_depth_mm <- samp_meta$Tf_volume_mL*1000/29460
#29460 is bucket opening size in mm^2
precip_bucket_avg_mm <- mean(samp_meta$Tf_depth_mm[samp_meta$Sample_type == "Precipitation"])
precip_bucket_avg_mL <- mean(samp_meta$Tf_volume_mL[samp_meta$Sample_type == "Precipitation"])
samp_meta$Tf_depth_pct_buckets <- samp_meta$Tf_depth_mm*100 / precip_bucket_avg_mm
samp_meta$Tf_depth_pct_gauge <- samp_meta$Tf_depth_mm*100 / 28.956
#rain gauge on NU campus measured total 1.14 in or 28.956 mm of rain

tf_meta <- samp_meta[samp_meta$Sample_type == "Throughfall", ]
#Used to be samp_meta[-(54:58), ] if anything changes
precip_meta <- samp_meta[samp_meta$Sample_type == "Precipitation", ]
mean(precip_meta$Tf_depth_mm)
std.error(precip_meta$Tf_depth_mm)
min(precip_meta$Tf_depth_mm)
max(precip_meta$Tf_depth_mm)

mean(precip_meta$DOC_mgL)
std.error(precip_meta$DOC_mgL)
mean(precip_meta$TDN_mgL)
std.error(precip_meta$TDN_mgL)

mean(tf_meta$Tf_depth_mm)
std.error(tf_meta$Tf_depth_mm)
min(tf_meta$Tf_depth_mm)
max(tf_meta$Tf_depth_mm)

mean(tf_meta$Tf_depth_pct_buckets)
std.error(tf_meta$Tf_depth_pct_buckets)
min(tf_meta$Tf_depth_pct_buckets)
max(tf_meta$Tf_depth_pct_buckets)
median(tf_meta$Tf_depth_pct_buckets)

mean(tf_meta$Tf_depth_pct_gauge)
std.error(tf_meta$Tf_depth_pct_gauge)
min(tf_meta$Tf_depth_pct_gauge)
max(tf_meta$Tf_depth_pct_gauge)
median(tf_meta$Tf_depth_pct_gauge)

mean(tf_meta$DOC_mgL, na.rm = TRUE)
std.error(tf_meta$DOC_mgL, na.rm = TRUE)
min(tf_meta$DOC_mgL, na.rm = TRUE)
max(tf_meta$DOC_mgL, na.rm = TRUE)

mean(tf_meta$TDN_mgL)
std.error(tf_meta$TDN_mgL)
min(tf_meta$TDN_mgL)
max(tf_meta$TDN_mgL)

mean_precip_DOC_flux <- mean(precip_meta$DOCflux_mg)
tf_meta$DOC_enrich <- tf_meta$DOCflux_mg / mean_precip_DOC_flux
mean_precip_TDN_flux <- mean(precip_meta$TDNflux_mg)
tf_meta$TDN_enrich <- tf_meta$TDNflux_mg / mean_precip_TDN_flux

mean(tf_meta$DOC_enrich, na.rm = TRUE)
std.error(tf_meta$DOC_enrich, na.rm = TRUE)
mean(tf_meta$TDN_enrich)
std.error(tf_meta$TDN_enrich)

tf_meta %>%
  group_by(Species) %>%
  summarise_at(vars(Tf_depth_mm), list(name = mean, std.error))|> print(n=36)

mean(tf_meta$Tf_depth_mm[tf_meta$Class == "Pinopsida"])
std.error(tf_meta$Tf_depth_mm[tf_meta$Class == "Pinopsida"])
mean(tf_meta$Tf_depth_mm[tf_meta$Class == "Dicotyledoneae"])
std.error(tf_meta$Tf_depth_mm[tf_meta$Class == "Dicotyledoneae"])

set.seed(0)
t.test(tf_meta$Tf_depth_mm[tf_meta$Class == "Pinopsida"],
       tf_meta$Tf_depth_mm[tf_meta$Class == "Dicotyledoneae"])
set.seed(0)
t.test(tf_meta$Tf_depth_mm[tf_meta$Leaf_type == "scales"],
       tf_meta$Tf_depth_mm[tf_meta$Leaf_type == "needles"])

mean(tf_meta$Tf_depth_mm[tf_meta$Leaf_type == "needles"])
std.error(tf_meta$Tf_depth_mm[tf_meta$Leaf_type == "needles"])
mean(tf_meta$Tf_depth_mm[tf_meta$Leaf_type == "scales"])
std.error(tf_meta$Tf_depth_mm[tf_meta$Leaf_type == "scales"])

set.seed(0)
t.test(tf_meta$Tf_depth_mm[tf_meta$Order == "Magnoliales"],
       tf_meta$Tf_depth_mm[tf_meta$Order != "Magnoliales"])

set.seed(0)
t.test(tf_meta$Tf_depth_mm[tf_meta$Genus == "Gleditsia"],
       tf_meta$Tf_depth_mm[tf_meta$Genus != "Gleditsia"])

set.seed(0)
t.test(tf_meta$Tf_depth_mm[tf_meta$Genus == "Gleditsia"],
       tf_meta$Tf_depth_mm[tf_meta$Class != "Pinopsida" & tf_meta$Genus != "Gleditsia"])
#small dense broadleaves aren't significantly like needles

#DOC and TDN conc and flux values by biological grouping
tf_meta %>%
  group_by(Species) %>%
  summarise_at(vars(DOC_mgL), list(name = mean, std.error), na.rm = TRUE)|> print(n=36)

set.seed(0)
t.test(tf_meta$DOC_mgL[tf_meta$Class == "Pinopsida"],
       tf_meta$DOC_mgL[tf_meta$Class == "Dicotyledoneae"])

mean(tf_meta$DOC_mgL[tf_meta$Class == "Pinopsida"], na.rm = TRUE)
std.error(tf_meta$DOC_mgL[tf_meta$Class == "Pinopsida"], na.rm = TRUE)
mean(tf_meta$DOC_mgL[tf_meta$Class == "Dicotyledoneae"], na.rm = TRUE)
std.error(tf_meta$DOC_mgL[tf_meta$Class == "Dicotyledoneae"], na.rm = TRUE)

set.seed(0)
t.test(tf_meta$DOCflux_mg[tf_meta$Class == "Pinopsida"],
       tf_meta$DOCflux_mg[tf_meta$Class == "Dicotyledoneae"])

tf_meta %>%
  group_by(Species) %>%
  summarise_at(vars(TDN_mgL), list(name = mean, std.error))|> print(n=36)

set.seed(0)
t.test(tf_meta$TDN_mgL[tf_meta$Class == "Pinopsida"],
       tf_meta$TDN_mgL[tf_meta$Class == "Dicotyledoneae"])

mean(tf_meta$TDN_mgL[tf_meta$Class == "Pinopsida"], na.rm = TRUE)
std.error(tf_meta$TDN_mgL[tf_meta$Class == "Pinopsida"], na.rm = TRUE)
mean(tf_meta$TDN_mgL[tf_meta$Class == "Dicotyledoneae"], na.rm = TRUE)
std.error(tf_meta$TDN_mgL[tf_meta$Class == "Dicotyledoneae"], na.rm = TRUE)

set.seed(0)
t.test(tf_meta$TDNflux_mg[tf_meta$Class == "Pinopsida"],
       tf_meta$TDNflux_mg[tf_meta$Class == "Dicotyledoneae"])

set.seed(0)
t.test(tf_meta$DOC_mgL[tf_meta$Order == "Magnoliales"],
       tf_meta$DOC_mgL[tf_meta$Order != "Magnoliales"])

set.seed(0)
t.test(tf_meta$DOC_mgL[tf_meta$Genus == "Acer"],
       tf_meta$DOC_mgL[tf_meta$Genus != "Acer"])

mean(tf_meta$DOC_mgL[tf_meta$Order == "Magnoliales"], na.rm = TRUE)
std.error(tf_meta$DOC_mgL[tf_meta$Order == "Magnoliales"], na.rm = TRUE)
mean(tf_meta$DOC_mgL[tf_meta$Order != "Magnoliales"], na.rm = TRUE)
std.error(tf_meta$DOC_mgL[tf_meta$Order != "Magnoliales"], na.rm = TRUE)

set.seed(0)
t.test(tf_meta$DOCflux_mg[tf_meta$Order == "Magnoliales"],
       tf_meta$DOCflux_mg[tf_meta$Order != "Magnoliales"])

mean(tf_meta$DOCflux_mg[tf_meta$Order == "Magnoliales"], na.rm = TRUE)
std.error(tf_meta$DOCflux_mg[tf_meta$Order == "Magnoliales"], na.rm = TRUE)
mean(tf_meta$DOCflux_mg[tf_meta$Order != "Magnoliales"], na.rm = TRUE)
std.error(tf_meta$DOCflux_mg[tf_meta$Order != "Magnoliales"], na.rm = TRUE)

set.seed(0)
t.test(tf_meta$TDN_mgL[tf_meta$Order == "Magnoliales"],
       tf_meta$TDN_mgL[tf_meta$Order != "Magnoliales"])

set.seed(0)
t.test(tf_meta$TDNflux_mg[tf_meta$Order == "Magnoliales"],
       tf_meta$TDNflux_mg[tf_meta$Order != "Magnoliales"])

cor(tf_meta$DBH_cm, tf_meta$DOCflux_mg, use = "complete.obs")
cor.test(tf_meta$DBH_cm, tf_meta$DOCflux_mg, use = "complete.obs")

#Figure 1: dilution curve plot
#first separate plots for DOC and TDN side by side
par(mfrow = c(1, 2))
#DOC on left
DOC_linear_fit <- lm(log(DOC_mgL) ~ Tf_depth_pct_buckets, data = tf_meta)
a_start_DOC <- exp(coef(DOC_linear_fit)[1])
b_start_DOC <- coef(DOC_linear_fit)[2]

DOC_exp_fit <- nls(DOC_mgL ~ a * exp(b * Tf_depth_pct_buckets),
                   data = tf_meta,
                   start = list(a = a_start_DOC, b = b_start_DOC))
DOC_smooth_x <- seq(5.59, 115.57, length.out = 300)
DOC_smooth_y <- predict(DOC_exp_fit, newdata = data.frame(Tf_depth_pct_buckets = DOC_smooth_x))

plot(x = (tf_meta$Tf_depth_pct_buckets), y = tf_meta$DOC_mgL,
     xlim = c(0, 120), ylim = c(0, 140),
     xlab = "Throughfall volume (% of rainfall)",
     ylab = expression(paste("DOC concentration (mg-C L"^"-1",")")))

lines(DOC_smooth_x, DOC_smooth_y, col = "blue", lwd = 1.5)

#legend
DOC_coef <- coef(DOC_exp_fit)
#a=67.15, b=-0.036
#r squared value
DOC_resid <- sum(residuals(DOC_exp_fit)^2)
DOC_res_tot <- sum((tf_meta$DOC_mgL - mean(tf_meta$DOC_mgL, na.rm = TRUE))^2, na.rm = TRUE)
DOC_r2 <- 1 - (DOC_resid / DOC_res_tot)
#r2 = 0.217
#p value
#77 samples have DOC data
DOC_null <- lm(DOC_mgL ~ 1, data = tf_meta)
DOC_f_stat <- ((DOC_res_tot - DOC_resid) / 1) /
  (DOC_resid / (77-2))
DOC_p_val <- pf(DOC_f_stat, df1 = 1, df2 = 77-2, lower.tail = FALSE)
#p=1.95*10^-5

legend("topright", legend = c(expression(paste("y = 67.15 * e"^"-0.036x","")),
                              expression(paste("r"^"2"," = 0.217")),
                              expression(paste("p = 1.95 * 10"^"-5",""))),
       inset = c(-0.4, 0), y.intersp = 0.4,
       col = c("blue", NA, NA), lwd = c(1.5, NA, NA), bty = "n", cex = 0.75)
text(x=5, y=133, labels = "a", cex = 2)

#TDN on right
TDN_linear_fit <- lm(log(TDN_mgL) ~ Tf_depth_pct_buckets, data = tf_meta)
a_start_TDN <- exp(coef(TDN_linear_fit)[1])
b_start_TDN <- coef(TDN_linear_fit)[2]

TDN_exp_fit <- nls(TDN_mgL ~ a * exp(b * Tf_depth_pct_buckets),
                   data = tf_meta,
                   start = list(a = a_start_TDN, b = b_start_TDN))
TDN_smooth_x <- seq(5.59, 115.57, length.out = 300)
TDN_smooth_y <- predict(TDN_exp_fit, newdata = data.frame(Tf_depth_pct_buckets = TDN_smooth_x))

plot(x = (tf_meta$Tf_depth_pct_buckets), y = tf_meta$TDN_mgL,
     xlim = c(0, 120), ylim = c(0, 42),
     xlab = "Throughfall volume (% of rainfall)",
     ylab = expression(paste("TDN concentration (mg-N L"^"-1",")")))

lines(TDN_smooth_x, TDN_smooth_y, col = "blue", lwd = 1.5)

#legend
TDN_coef <- coef(TDN_exp_fit)
#a=16.81, b=-0.055
#r squared value
TDN_resid <- sum(residuals(TDN_exp_fit)^2)
TDN_res_tot <- sum((tf_meta$TDN_mgL - mean(tf_meta$TDN_mgL))^2)
TDN_r2 <- 1 - (TDN_resid / TDN_res_tot)
#r2 = 0.218
#p value
TDN_null <- lm(TDN_mgL ~ 1, data = tf_meta)
TDN_f_stat <- ((TDN_res_tot - TDN_resid) / 1) /
  (TDN_resid / (100-2))
TDN_p_val <- pf(TDN_f_stat, df1 = 1, df2 = 100-2, lower.tail = FALSE)
#p=9.71*10^-7

legend("topright", legend = c(expression(paste("y = 16.81 * e"^"-0.055x","")),
                              expression(paste("r"^"2"," = 0.218")),
                              expression(paste("p = 9.71 * 10"^"-7",""))),
       inset = c(-0.4, 0), y.intersp = 0.4,
       col = c("blue", NA, NA), lwd = c(1.5, NA, NA), bty = "n", cex = 0.75)
text(x=5, y=40, labels = "b", cex = 2)
#looks good
#need to add a and b labels to them if using this


#then both on one plot
par(mfrow = c(1, 1))
par(mar = c(5, 5, 3, 5))
#start by plotting DOC
plot(x = (tf_meta$Tf_depth_pct_buckets), y = tf_meta$DOC_mgL,
     pch = 1, col = "blue",
     xlim = c(0, 120), ylim = c(0, 140),
     xlab = "Throughfall volume (% of rainfall)",
     ylab = expression(paste("DOC concentration (mg-C L"^"-1",")")))
lines(DOC_smooth_x, DOC_smooth_y, col = "blue", lwd = 1.5)
#TDN over it now
par(new = TRUE)
plot(x = (tf_meta$Tf_depth_pct_buckets), y = tf_meta$TDN_mgL,
     pch = 2, col = "red",
     axes = FALSE, xlab = "", ylab = "",
     xlim = c(0, 120), ylim = c(0, 42))
lines(TDN_smooth_x, TDN_smooth_y, col = "red", lwd = 1.5)
axis(side = 4)
mtext(expression(paste("TDN concentration (mg-N L"^"-1",")")),
      side = 4, line = 3)

legend("topright",
       legend = c(expression(paste("DOC: y = 67.15 * e"^"-0.036x","")),
                  expression(paste("TDN: y = 16.81 * e"^"-0.055x",""))),
       inset = c(-0.15, 0),
       pch = c(1, 2), col = c("blue", "red"),
       lwd = 1.5, bty = "n")



#correlations between all numeric variables
ms_bulk <- read_csv(file = 'data/MS_bulk_parameters.csv')
ms_bulk_tf <- ms_bulk[-(54:58), ]
tf_meta_w_ms_bulk <- merge(tf_meta, ms_bulk_tf)

tf_meta_numeric <- tf_meta_w_ms_bulk[sapply(tf_meta_w_ms_bulk, is.numeric)]
correlations <- cor(tf_meta_numeric, use = "everything")
write.csv(correlations, "correlations.csv", row.names = TRUE)
#doing DOC separately since it has NA values
correlations_DOC <- cor(tf_meta_numeric, use = "complete.obs")
write.csv(correlations_DOC, "correlations_DOC.csv", row.names = TRUE)
#testing for significance to report
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$TDN_mgL, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Molecular.mass_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$N_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$S_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$P_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Heteroatoms_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$CHO_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$H.C_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$O.C_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$AI_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$AImod_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$DBE_avg, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Arom_O_poor_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$High_unsat_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$High_unsat_O_rich_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$High_unsat_O_poor_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Unsat_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Unsat_O_rich_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Unsat_O_poor_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Unsat_with_N_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Sat_O_rich_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Sat_O_poor_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Sugar_like_pct, use = "complete.obs")
cor.test(tf_meta_numeric$DOC_mgL, tf_meta_numeric$Peptide_like_pct, use = "complete.obs")

set.seed(0)
t.test(tf_meta_w_ms_bulk$Molecular.mass_avg[tf_meta$Class == "Pinopsida"],
       tf_meta_w_ms_bulk$Molecular.mass_avg[tf_meta$Class == "Dicotyledoneae"])

cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Molecular.mass_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$N_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$S_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$P_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Heteroatoms_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$CHO_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$H.C_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$O.C_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$AI_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$AImod_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$DBE_avg)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Cond_arom_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Arom_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Arom_O_rich_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Arom_O_poor_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$High_unsat_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$High_unsat_O_rich_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$High_unsat_O_poor_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Unsat_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Unsat_O_rich_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Unsat_O_poor_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Unsat_with_N_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Sat_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Sat_O_rich_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Sat_O_poor_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Sugar_like_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$Peptide_like_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$CHON_pct)
cor.test(tf_meta_numeric$TDN_mgL, tf_meta_numeric$CHONS_pct)
cor.test(tf_meta_numeric$N_avg, tf_meta_numeric$S_avg)

# FTICRMS data# FTICRMS data# FTICRMS data
peak_meta <- read.csv(file = 'data/peak_metadata.csv')
#if want to do peak filtering, do it here
#add condensed aromatic rows
peak_meta$Condensed.aromatic <- NA
peak_meta$Condensed.aromatic[peak_meta$AI.mod > 0.67] <- 1
peak_meta$Condensed.aromatic[peak_meta$AI.mod <= 0.67] <- 0

peak_meta$Condensed.aromatic.O_rich <- NA
peak_meta$Condensed.aromatic.O_rich[peak_meta$Condensed.aromatic == 1 & peak_meta$Aromatic.O_rich == 1] <- 1
peak_meta$Condensed.aromatic.O_rich[peak_meta$Condensed.aromatic != 1 | peak_meta$Aromatic.O_rich != 1] <- 0

peak_meta$Condensed.aromatic.O_poor <- NA
peak_meta$Condensed.aromatic.O_poor[peak_meta$Condensed.aromatic == 1 & peak_meta$Aromatic.O_poor == 1] <- 1
peak_meta$Condensed.aromatic.O_poor[peak_meta$Condensed.aromatic != 1 | peak_meta$Aromatic.O_poor != 1] <- 0

peak_meta$Uncondensed.aromatic <- NA
peak_meta$Uncondensed.aromatic[peak_meta$Condensed.aromatic == 0 & peak_meta$Aromatic == 1] <- 1
peak_meta$Condensed.aromatic[peak_meta$Condensed.aromatic != 0 | peak_meta$Aromatic != 1] <- 0

peak_meta$Uncondensed.aromatic.O_rich <- NA
peak_meta$Uncondensed.aromatic.O_rich[peak_meta$Condensed.aromatic.O_rich == 0 & peak_meta$Aromatic.O_rich == 1] <- 1
peak_meta$Uncondensed.aromatic.O_rich[peak_meta$Condensed.aromatic.O_rich != 0 | peak_meta$Aromatic.O_rich != 1] <- 0

peak_meta$Uncondensed.aromatic.O_poor <- NA
peak_meta$Uncondensed.aromatic.O_poor[peak_meta$Condensed.aromatic.O_poor == 0 & peak_meta$Aromatic.O_poor == 1] <- 1
peak_meta$Uncondensed.aromatic.O_poor[peak_meta$Condensed.aromatic.O_poor != 0 | peak_meta$Aromatic.O_poor != 1] <- 0
#should be good, check output and make sure it makes sense and matches what was done in excel

intens_w_meta <- read_csv(file = 'data/intensities_trimmed.csv')
#this has formula chemistry parameters in there as well
#peaks not detected are in there as NA, not 0

intensities <- read_csv(file = 'data/intensities_only.csv')
#this is only mz, formula, and intensities

norm_intensities <- read.table(file = 'data/intensities_norm.csv', sep = ",",
                               header = TRUE, row.names = 1)
#this has metadata removed, undetected peaks are 0s
#In case I need them to be NA instead
#norm_intensities[norm_intensities == 0] <- NA

#write.table(norm_transp_intensities,
            file = "norm_transp_intensities.csv", sep = ",", row.names = FALSE)
#can use this to make sure it looks good at any point

#kevin's method
all_samples <- samp_meta %>% pull(Sample_ID)
#often don't want precipitation samples included
tf_samples <- samp_meta %>% 
  filter(Sample_type !="Precipitation") %>% pull(Sample_ID)
precip_samples <- samp_meta %>% 
  filter(Sample_type =="Precipitation") %>% pull(Sample_ID)
#for use later when hybrid species will be excluded
#don't end up using this
tf_nonhyb_samples <- samp_meta %>% 
  filter(PhyloName !="NA") %>% pull(Sample_ID)
#for use later when samples with NA DOC values will be excluded
tf_samples_w_DOC <- samp_meta %>% 
  filter(DOC_uM !="#N/A" & Sample_type != "Precipitation") %>% pull(Sample_ID)

norm_intens_t <- norm_intensities %>% 
  column_to_rownames("mz") %>%
  select(all_of(tf_samples)) %>%
  # mutate(across(everything(), ~
  #               ~(. - min(., na.rm=TRUE))/
  #                 range(., na.rm = TRUE))) %>% 
  #                 # (max(.) - min(.)))) %>% 
  # na_if(., 0))) %>% #zeros skew z-score
  as.matrix() %>% t()
#nmds needs species/mass data as columns across rows(samples), this got it in that format
#it doesn't need to be a symmetric input
# scales::rescale(., to = c(0,max(expr_norm2)))

#check my data format is good
is.matrix(norm_intens_t)
#True, that's good
rownames(norm_intens_t)
#row names are sample IDs- correct
#doesn't include mz or formula, I can go back and include these in tf_samples if needed
colnames(norm_intens_t)
#mz values- good
norm_intens_t[1, ]
#showing mz and the rel intensity in first sample- good
str(norm_intens_t)
#100 rows (sample IDs) by 10799 columns (mzs)- good
isSymmetric(norm_intens_t)
#false, that should be fine

#NMDS function
set.seed(3) #for reproduction, function uses random starts
nmds <- metaMDS(comm = norm_intens_t, distance = "bray", k = 3,
                autotransform = TRUE, na.rm = TRUE)
#did 20 runs, all had stress<0.1
## autotransform chose: Square root transformation and Wisconsin double standardization
## could not converge with 2 dimensions
nmds_plot <- plot(nmds)
stressplot <- stressplot(nmds) #low scatter indicates nmds ordination good representation of original data
#linear fit r2=.987, non-metric fit r2=.996
stress <- nmds$stress
#stress should be below 0.1; it's 0.062

#Get ordination data
#nmds$points contains the positions of each sample in nmds dimensions
#scores() vegan function that gets $points
ord_samp <- scores(nmds, "sites") %>% as_tibble(rownames = "Sample_ID") 
ord_mass <- scores(nmds, "species") %>%
  as_tibble(rownames = "mz") %>% mutate(across(mz, as.numeric))

#add ancillary environmental variables
ord_samp2 <- left_join(ord_samp, samp_meta, by = "Sample_ID")
ord_mass2 <- left_join(ord_mass, peak_meta, by = "mz")

#plot samples in NMDS space
p_samp <- ggplot()+
  geom_point(data= ord_samp2,
             aes(x=NMDS1, y=NMDS2, color = Leaf_type)) +
  theme(legend.position = "none")
p_samp

#plot masses in NMDS space
p_mass <- ggplot()+
  geom_point(data= ord_mass2,
             aes(x=NMDS1, y=NMDS2), color= "grey", alpha = 0.5)+
  # scale_color_distiller()+
  theme_bw()
# coord_fixed()
p_mass

#bray-curtis distance calculation
#use transposed intensities, otherwise it will calculate distances between pairs of peaks
#hybrid samples excluded for use w/ phylogenetic data
bray_input_nonhyb <- norm_intens_t[tf_nonhyb_samples,]
peak_dist_nonhyb <- vegdist(bray_input_nonhyb, method = "bray")
peak_dist_matrix_nonhyb <- as.matrix(peak_dist_nonhyb)
#all tf samples
bray_input_all <- norm_intens_t[tf_samples,]
peak_dist_all <- vegdist(bray_input_all, method = "bray")
peak_dist_matrix_all <- as.matrix(peak_dist_all)
#only tf samples with DOC
bray_input_DOC <- norm_intens_t[tf_samples_w_DOC, ]
peak_dist_DOC <- vegdist(bray_input_DOC, method = "bray")
peak_dist_matrix_DOC <- as.matrix(peak_dist_DOC)


#phylogenetic relationship distance matrix
species <- samp_meta$`Latin`[samp_meta$`Latin` != "Precipitation"]
taxa <- tnrs_match_names(species)
#this has problems with the two hybrid species
non_hyb_taxa <- taxa[-c(18, 30),]
subtree <- tol_induced_subtree(ott_ids = non_hyb_taxa$ott_id)
phylo_dist_matrix <- cophenetic.phylo(subtree)
phylo_dist_matrix
plot(subtree)
#done, tree looks good

#now trying with hybrids as average of the parents
hybrid_parents <- c("Platanus orientalis", "Platanus occidentalis", "Magnolia denudata", "Magnolia liliiflora")
species_with_hyb_parents <- species
species_with_hyb_parents[101:104] <- hybrid_parents
species_hyb_as_parents <- species_with_hyb_parents[-c(48, 49, 50, 82, 83, 84)]
taxa_hyb_as_parents <- tnrs_match_names(species_hyb_as_parents)
tree_hyb_as_parents <- tol_induced_subtree(ott_ids = taxa_hyb_as_parents$ott_id)
plot(tree_hyb_as_parents)
phylo_dist_matrix_hyb_parents <- cophenetic.phylo(tree_hyb_as_parents)
phylo_dist_matrix_hyb_parents
write.csv(phylo_dist_matrix_hyb_parents, "phylo_dist_matrix_hyb_parents.csv", row.names = TRUE)
#looks good

#next, add hybrids back in and calculate them as average of parents
expanded_matrix <- matrix(0, nrow = 40, ncol = 40)
expanded_row_names <- c(rownames(phylo_dist_matrix_hyb_parents), "Platanus_acerifolia", "Magnolia_x_soulangeana")
rownames(expanded_matrix) <- expanded_row_names
colnames(expanded_matrix) <- expanded_row_names
known_species <- rownames(phylo_dist_matrix_hyb_parents)
expanded_matrix[known_species, known_species] <- phylo_dist_matrix_hyb_parents
expanded_matrix[39, ] <- (expanded_matrix[24, ] + expanded_matrix[25, ]) / 2
expanded_matrix[, 39] <- (expanded_matrix[, 24] + expanded_matrix[, 25]) / 2
expanded_matrix[40, ] <- (expanded_matrix[26, ] + expanded_matrix[28, ]) / 2
expanded_matrix[, 40] <- (expanded_matrix[, 26] + expanded_matrix[, 28]) / 2
expanded_matrix[39, 39] <- 0
expanded_matrix[40, 40] <- 0
isSymmetric(expanded_matrix)
write.csv(expanded_matrix, "expanded_matrix.csv", row.names = TRUE)

#looks good, now can remove parents
phylo_dist_matrix_final_species <- expanded_matrix[-c(24, 25, 26, 28), -c(24, 25, 26, 28)]
write.csv(phylo_dist_matrix_final_species, "phylo_dist_matrix_final_species.csv", row.names = TRUE)
#everything has copied over right so far
#create clustering diagram
phylo_dist <- as.dist(phylo_dist_matrix_final_species)
cluster_nj_phylo <- nj(phylo_dist)
plot(cluster_nj_phylo, type = "phylogram")
#didn't work

hclust_phylo_dist <- hcphylo_disthclust_phylo_dist <- hclust(phylo_dist, method = "ward.D")
plot(hclust_phylo_dist)
#looks terrible
#get rid of tiny exponents for transfer to jmp
phylo_dist_for_jmp <- phylo_dist_matrix_final_species*10^300
phylo_dist_for_jmp_2 <- phylo_dist_for_jmp*10^12
phylo_dist_jmp <- as.dist(phylo_dist_for_jmp_2)
hclust_phylo_dist_jmp <- hclust(phylo_dist_jmp, method = "average")
plot(hclust_phylo_dist_jmp, cex = 0.7)
cladogram_labels <- c("Koelreuteria paniculata", "Acer palmatum", "Acer griseum", "Acer rubrum", "Acer platanoides",
                      "Tilia cordata", "Ulmus americana", "Zelkova serrata", "Pyrus calleryana", "Gleditsia triacanthos",
                      "Betula papyrifera", "Betula nigra", "Quercus frainetto", "Quercus rubra", "Quercus palustris", 
                      "Quercus acutissima", "Fagus sylvatica", "Liquidambar styraciflua", "Cercidiphyllum japonicum", "Ilex opaca",
                      "Cornus kousa", "Cornus florida", "Nyssa sylvatica", "Magnolia kobus", "Magnolia stellata",
                      "Magnolia virginiana", "Juniperus virginiana", "Callitropsis nootkatensis", "Chamaecyparis obtusa", "Thuja occidentalis",
                      "Metasequoia glyptostroboides", "Pinus resinosa", "Pinus strobus", "Abies concolor", "Platanus acerifolia", "Magnolia x soulangeana")
cluster_nj_phylo2 <- nj(phylo_dist_jmp)
cluster_nj_phylo2$tip.label <- cladogram_labels
plot(cluster_nj_phylo2, type = "cladogram")
plot(cluster_nj_phylo2, type = "tidy")
write.csv(phylo_dist_for_jmp, "phylo_dist_for_jmp.csv", row.names = TRUE)

#now need to expand it to the same size as sample distance matrix
phylo_dist_matrix_final_sample <- matrix(0, nrow = 100, ncol = 100)
rownames(phylo_dist_matrix_final_sample) <- tf_samples
colnames(phylo_dist_matrix_final_sample) <- tf_samples
tf_meta$PhyloName[48:50] <- "Platanus_acerifolia"
tf_meta$PhyloName[82:84] <- "Magnolia_x_soulangeana"
for (i in 1:100) {
  sample_i <- tf_samples[i]
  species_i <- tf_meta$PhyloName[tf_meta$Sample_ID == sample_i]
  for (j in 1:100) {
    sample_j <- tf_samples[j]
    species_j <- tf_meta$PhyloName[tf_meta$Sample_ID == sample_j]
    phylo_dist_matrix_final_sample[i,j] <- phylo_dist_matrix_final_species[species_i, species_j]
  }
}
write.csv(phylo_dist_matrix_final_sample, "phylo_dist_matrix_final_sample.csv", row.names = TRUE)
isSymmetric(phylo_dist_matrix_final_sample)

#now need to expand it to the same size as sample distance matrix
#which is 94x94 once hybrids are removed
large_phylo_matrix <- matrix(0, nrow = 94, ncol = 94)
rownames(large_phylo_matrix) <- tf_nonhyb_samples
colnames(large_phylo_matrix) <- tf_nonhyb_samples
for (i in 1:94) {
  sample_i <- tf_nonhyb_samples[i]
  species_i <- samp_meta$PhyloName[samp_meta$Sample_ID == sample_i]
  for (j in 1:94) {
    sample_j <- tf_nonhyb_samples[j]
    species_j <- samp_meta$PhyloName[samp_meta$Sample_ID == sample_j]
    large_phylo_matrix[i,j] <- phylo_dist_matrix[species_i, species_j]
  }
}

#rda on norm intensities
#analyzes peaks, can look for peaks associated w/ env variables
#need samples as rows, peaks and parameters as columns
#peaks same input as bray, just need to remove precip from metadata
#columns wanted: tf volume, DOC and TDN conc, DBH,
#and distance, sine, and cosine for buildings, roads, and trains
rda_metadata <- as.data.frame(samp_meta[-(54:58), c(14, 15, 18, 22, 43, 45, 46, 47, 49, 50, 51, 53, 54)])
#fix that to be names instead of numbers
#shouldn't need these because did them above
rda_metadata$DBH_cm <- as.numeric(rda_metadata$DBH_cm)
rda_metadata$BLDG_dist <- as.numeric(rda_metadata$BLDG_dist)
rda_metadata$ROAD_dist <- as.numeric(rda_metadata$ROAD_dist)
rda_metadata$TRAIN_dist <- as.numeric(rda_metadata$TRAIN_dist)
rda_metadata$DOC_uM <- as.numeric(rda_metadata$DOC_uM)
rda_metadata$BLDG_sin <- as.numeric(rda_metadata$BLDG_sin)
rda_metadata$BLDG_cos <- as.numeric(rda_metadata$BLDG_cos)
rda_metadata$ROAD_sin <- as.numeric(rda_metadata$ROAD_sin)
rda_metadata$ROAD_cos <- as.numeric(rda_metadata$ROAD_cos)
rda_metadata$TRAIN_sin <- as.numeric(rda_metadata$TRAIN_sin)
rda_metadata$TRAIN_cos <- as.numeric(rda_metadata$TRAIN_cos)

rownames(rda_metadata) <- tf_samples
rda_peaks <- rda(bray_input_all ~ ., data = rda_metadata, na.action = na.exclude)
summary(rda_peaks)
summary(rda_peaks)$cont
anova(rda_peaks)
anova(rda_peaks, by = "terms")
plot(rda_peaks)

#db-rda on sample distance matrix
#analyzes samples for similarities/differences between them
dbrda <- dbrda(peak_dist_matrix_all ~ ., data = rda_metadata, na.action = na.exclude)
summary(dbrda)
summary(dbrda)$cont
plot(dbrda)
dbrda$CCA$tot.chi / dbrda$tot.chi
dbrda$CA$tot.chi / dbrda$tot.chi
#these say 29.5% of variance constrained, 70.5% unconstrained
#with angles added in as sin and cos, 36.8% constrained, 63.2% unconstrained
anova(dbrda, permutations = 999)
anova(dbrda, by = "terms", permutations = 999)
#DOC and TDN significant at P<.001 level
#tf volume and dbh significant at P<.05 level
#building distance on the edge of significance, sometimes reaches .05 level, sometimes not
#adding angles didn't change this, angles aren't significant at all
anova(dbrda, by = "axis", permutations = 999)
#first two axes are significant at P<.001 level, others not significant

#plot labeled with sample ID
plot(dbrda, type = "n")
points(dbrda, display = "sites", pch = 19, col = "blue", cex = 0.8)
text(dbrda, display = "sites", labels = rownames(rda_metadata), pos = 3, cex = 0.6, offset = 0.3)
text(dbrda, display = "bp", col = "red", cex = 0.8)

rda_meta_clean <- rda_metadata
clean_colnames <- c("Volume", "[DOC]", "[TDN]", "Tree size",
                    "Building distance", "Building sin", "Building cos",
                    "Road distance", "Road sin", "Road cos",
                    "Train distance", "Train sin", "Train cos")
colnames(rda_meta_clean) <- clean_colnames
dbrda_clean <- dbrda(peak_dist_matrix_all ~ ., data = rda_meta_clean, na.action = na.exclude)
plot(dbrda_clean, type = "n")
text(dbrda_clean, display = "bp", col = "blue", cex = 0.8)

#dbrda just on potentially significant variables
rda_meta_signif <- rda_meta_clean[, 1:5]
dbrda_signif <- dbrda(peak_dist_matrix_all ~ ., data = rda_meta_signif, na.action = na.exclude)
plot(dbrda_signif, type = "n")
text(dbrda_signif, display = "bp", col = "blue", cex = 0.8)
points(dbrda_signif, display = "sites", pch = 20, col = "black", cex = 0.8)
#want to project the rainwater samples into the same space
norm_intens_t_w_precip <- norm_intensities %>% 
  column_to_rownames("mz") %>%
  select(all_of(all_samples)) %>%
  # mutate(across(everything(), ~
  #               ~(. - min(., na.rm=TRUE))/
  #                 range(., na.rm = TRUE))) %>% 
  #                 # (max(.) - min(.)))) %>% 
  # na_if(., 0))) %>% #zeros skew z-score
  as.matrix() %>% t()
peak_dist_w_precip <- vegdist(norm_intens_t_w_precip, method = "bray")
peak_dist_matrix_w_precip <- as.matrix(peak_dist_w_precip)
tf_scores <- scores(dbrda_signif)
all_scores <- capscale(peak_dist_w_precip ~ 1)
positions_all <- scores(all_scores, display = "sites", choices = c(1, 2))
positions_precip <- positions_all[precip_samples, ]
plot(dbrda_signif, type = "n", xlim = c(-3, 3), ylim = c(-2., 2))
text(dbrda_signif, display = "bp", col = "blue", cex = 0.8)
points(dbrda_signif, display = "sites", pch = 20, col = "black", cex = 0.8)
points(positions_precip, pch = 17, col = "red", cex = 0.8)
legend("topright",
       legend = c("Throughfall", "Rain"),
       col = c("black", "red"), pch = c(20, 17),
       cex = 0.8, bty = "n", inset = c(0.08, 0), y.intersp = 0.6)

#want to see if conifers/magnolias cluster uniquely
conifers <- samp_meta %>% 
  filter(Class =="Pinopsida") %>% pull(Sample_ID)
positions_conifer <- positions_all[conifers, ]
magnolias <- samp_meta %>% 
  filter(Order == "Magnoliales") %>% pull(Sample_ID)
positions_magnolia <- positions_all[magnolias, ]
other_tf <- samp_meta %>% 
  filter(Class !="Pinopsida" & Order != "Magnoliales" & Sample_type != "Precipitation") %>% pull(Sample_ID)
positions_other <- positions_all[other_tf, ]

order_groups <- ifelse(tf_meta$Order == "Pinales", "Pinales",
                       ifelse(tf_meta$Order == "Magnoliales", "Magnoliales",
                              "Other"))
order_colors <- c("Pinales" = "forestgreen", "Magnoliales" = "purple", "Other" = "gray40")
order_shapes <- c("Pinales" = 22, "Magnoliales" = 25, "Other" = 21)
plot(dbrda_signif, type = "n", xlim = c(-3, 3), ylim = c(-2., 1.9))
text(dbrda_signif, display = "bp", col = "blue", cex = 0.8, lwd = 1.5)
for(group in c("Other", "Pinales", "Magnoliales")) {
  indices <- which(order_groups == group)
  points(dbrda_signif, display = "sites", select = indices,
         col = order_colors[group], bg = order_colors[group],
         pch = order_shapes[group], cex = 0.8)
}
points(positions_precip, pch = 17, col = "red", cex = 0.8)
legend("topright",
       legend = c("Pinales", "Magnoliales", "Other", "Rain"),
       pch = c(22, 25, 21, 17),
       col = c("forestgreen", "purple", "gray40", "red"),
       pt.bg = c("forestgreen", "purple", "gray40", "red"),
       cex = 0.8, bty = "n", inset = c(0.08, 0), y.intersp = 0.6)
#adding density ellipses around the four groups
ordiellipse(dbrda_signif, groups = order_groups,
            col = c("purple", "gray40", "forestgreen"), border = c("purple", "gray40", "forestgreen"),
            draw = "polygon", alpha = 50, lwd = 1,
            kind = "sd", conf = 0.95)

all_scores_ord <- rbind(scores(dbrda_signif, display = "sites", choices = c(1, 2)),
                        positions_precip)
precip_group <- c(order_groups, rep("Precipitation", 5))
ordiellipse(all_scores_ord, groups = precip_group, show.groups = "Precipitation",
            col = "red", border = "red",
            draw = "polygon", alpha = 50, lwd = 1,
            kind = "sd", conf = 0.95)
#display is being weird and not filling in this ellipse when I go to export it
#so I'll make it and draw it manually
precip_ellipse <- ordiellipse(all_scores_ord, groups = precip_group, show.groups = "Precipitation",
                              draw = "none", kind = "sd", conf = 0.95)
precip_ellipse_coords <- precip_ellipse[[1]]
polygon(precip_ellipse_coords,
        col = adjustcolor("red", alpha.f = 0.2),
        border = "red", lwd = 1)
#still not working, I'll deal with it later


precip_meta_rda_clean <- precip_meta[, c("Tf_volume_mL", "DOC_uM", "TDN_uM", "DBH_cm", "BLDG_dist")]
colnames(precip_meta_rda_clean) <- c("Volume", "[DOC]", "[TDN]", "Tree size", "Building distance")
precip_meta_rda_clean$`Tree size` <- as.numeric(precip_meta_rda_clean$`Tree size`)
precip_meta_rda_clean$`Building distance` <- as.numeric(precip_meta_rda_clean$`Building distance`)
precip_scores <- predict(dbrda_signif, newdata = precip_meta_rda_clean, type = "lc", na.action = na.exclude)


anova(dbrda_signif, permutations = 999)
anova(dbrda_signif, by = "terms", permutations = 999)
anova(dbrda_signif, by = "axis", permutations = 999)
#same result as with all variables
#DOC and TDN significant at P<.001 level
#tf volume and dbh significant at P<.05 level
#building distance almost significant but 0.5<P<0.1
dbrda_signif$CCA$tot.chi / dbrda_signif$tot.chi
#constrains 27.5% of variance

#test individual environmental variables for significance and variation explained
no_NA_vars <- colnames(rda_metadata)[!colnames(rda_metadata) %in% "DOC_uM"]
individual_results <- list()
for(var in no_NA_vars) {
  formula <- as.formula(paste("peak_dist_matrix_all ~", var))
  single_dbrda <- dbrda(formula, data = rda_metadata)
  sig_test <- anova(single_dbrda, permutations = 999)
  var_explained <- single_dbrda$CCA$tot.chi / single_dbrda$tot.chi
  individual_results[[var]] <- list(
    variable = var,
    variance_explained = var_explained,
    p_value = sig_test$`Pr(>F)`[1],
    F_statistic = sig_test$F[1]
  )
}
explanatory_vars <- do.call(rbind, lapply(individual_results, data.frame))
explanatory_vars <- explanatory_vars[order(explanatory_vars$variance_explained, decreasing = TRUE), ]
print(explanatory_vars)
#this is showing only TDN being significant on its own, DBH is close (.058)
#now running dbrda on DOC alone because NA values interfere
DOC_test <- !is.na(rda_metadata[["DOC_uM"]])
DOC_samples <- rda_metadata[DOC_test, ]
dist_mx_w_DOC <- peak_dist_matrix_all[DOC_test, DOC_test]
formula_DOC <- as.formula(paste("dist_mx_w_DOC ~", "DOC_uM"))
dbrda_DOC <- dbrda(formula_DOC, data = DOC_samples)
sig_test_DOC <- anova(dbrda_DOC, permutations = 999)
var_explained_DOC <- dbrda_DOC$CCA$tot.chi / dbrda_DOC$tot.chi
print(paste("DOC - Variance Explained:", var_explained_DOC))
print(paste("p-value:", sig_test_DOC$`Pr(>F)`[1]))
print(paste("Samples Used:", sum(DOC_test), "out of", nrow(rda_metadata)))
#significant at .001 level, explains 11% of variation
#this might not be right thought because I just pulled numbers from the full distance matrix
#using the DOC-only distance matrix from early on
formula_DOC_2 <- as.formula(paste("peak_dist_matrix_DOC ~", "DOC_uM"))
dbrda_DOC_2 <- dbrda(formula_DOC_2, data = DOC_samples)
sig_test_DOC_2 <- anova(dbrda_DOC_2, permutations = 999)
var_explained_DOC_2 <- dbrda_DOC_2$CCA$tot.chi / dbrda_DOC_2$tot.chi
print(paste("DOC - Variance Explained:", var_explained_DOC_2))
print(paste("p-value:", sig_test_DOC_2$`Pr(>F)`[1]))
#exact same result as above: 11.01% of variance explained, significant at .001

#Figure 2: db-rda for environmental variables explaining mass spec similarities
#want precip samples in there too so can project them onto dbrda plot
norm_intens_t_w_precip <- norm_intensities %>% 
  column_to_rownames("mz") %>%
  select(all_of(all_samples)) %>%
  # mutate(across(everything(), ~
  #               ~(. - min(., na.rm=TRUE))/
  #                 range(., na.rm = TRUE))) %>% 
  #                 # (max(.) - min(.)))) %>% 
  # na_if(., 0))) %>% #zeros skew z-score
  as.matrix() %>% t()
peak_dist_w_precip <- vegdist(norm_intens_t_w_precip, method = "bray")
peak_dist_matrix_w_precip <- as.matrix(peak_dist_w_precip)
tf_scores <- scores()




#test for correlations between variables
cor_matrix <- cor(rda_metadata, use = "complete.obs")
cor_matrix
#this says only variable pair with strong correlation is train sin and cos (-0.95)
#the rest are 0.43 or less
corrplot(cor_matrix, method = "color", type = "upper",
         order = "hclust", tl.cex = 0.8, addCoef.col = "black", number.cex = 0.6)

#skipping VIF test for now because packages aren't working
#vif_dummy <- rnorm(nrow(rda_metadata))
#vif_model <- lm(vif_dummy ~ ., data = rda_metadata)
#vif_values <- vif(vif_model)
#ols_vif_tol(vif_model)
#VIF in jmp agreed with correlation matrix

#condition number as another test of collinearity
cond_numb <- kappa(cor(rda_metadata, use = "complete.obs"))
print(paste("Condition number:", cond_numb))
#>30 indicates serious collinearity, my condition number is 99.89


#analyzing phylogeny data for significance
#excluding phylum because it's no broader than class
#there are 2 different classes, 12 orders, 16 families, 23 genera, and 36 species
phenolog_meta <- as.data.frame(samp_meta[-(54:58), (5:9)], as.factor)
rownames(phenolog_meta) <- tf_samples
phenolog_vars <- colnames(samp_meta)[5:9]
phenolog_signif <- list()
for(var in phenolog_vars) {
  formula_ph <- as.formula(paste("peak_dist_matrix_all ~", var))
  single_dbrda_ph <- dbrda(formula_ph, data = phenolog_meta)
  sig_test_ph <- anova(single_dbrda_ph, permutations = 999)
  var_explained_ph <- single_dbrda_ph$CCA$tot.chi / single_dbrda_ph$tot.chi
  phenolog_signif[[var]] <- list(
    variable = var,
    variance_explained = var_explained_ph,
    p_value = sig_test_ph$`Pr(>F)`[1],
    F_statistic = sig_test_ph$F[1],
    df = sig_test_ph$Df[1]
  )
}
phenolog_signif_results <- do.call(rbind, lapply(phenolog_signif, data.frame))
phenolog_signif_results <- phenolog_signif_results[order(phenolog_signif_results$variance_explained, decreasing = TRUE), ]
print(phenolog_signif_results)

norm_intens_t_meta <- cbind(phenolog_meta, norm_intens_t)
write.csv(norm_intens_t_meta, "norm_intens_t_meta.csv", row.names = TRUE)
#permanova on phenological data- this shows the same results as dbrda above
for(var in phenolog_vars) {
  formula <- as.formula(paste("peak_dist_matrix_all ~", var))
  permanova_result <- adonis2(formula, data = phenolog_meta, permutations = 999)
  print(paste("Variable:", var))
  print(permanova_result)
  cat("\n")
}
#all 5 levels separately are significant at .001 level
#class explains 5.7% of variation, order 26.0%, family 33.2%, genus 40.5%, and species 52.9%
#there's a good linear relationship between # of categories in a phenological level and % of variance explained

#now look at them sequentially instead of individually
phenolog_meta[phenolog_vars] <- lapply(phenolog_meta[phenolog_vars], as.factor)
sequential_phenolog <- list()
remaining_var <- 1.0
for(i in 1:5) {
  current_var <- phenolog_vars[i]
  if(i == 1) {
    formula_phenolog <- as.formula(paste("peak_dist_matrix_all ~", current_var))
    dbrda_phenolog <- dbrda(formula_phenolog, data = phenolog_meta)
  } else {
    conditioning_vars <- paste(phenolog_vars[1:(i-1)], collapse = " + ")
    formula_phenolog <- as.formula(paste("peak_dist_matrix_all ~", current_var, "+ Condition(", conditioning_vars, ")"))
    dbrda_phenolog <- dbrda(formula_phenolog, data = phenolog_meta)
  }
  if(i == 1) {
    var_explained_total <- dbrda_phenolog$CCA$tot.chi / dbrda_phenolog$tot.chi
    var_explained_remaining <- var_explained_total
  } else {
    var_explained_total <- dbrda_phenolog$CCA$tot.chi / dbrda_phenolog$tot.chi
    var_explained_remaining <- var_explained_total / remaining_var
  }
  remaining_var <- remaining_var - var_explained_total
  sig_test_phenolog <- anova(dbrda_phenolog, permutations = 999)
  sequential_phenolog[[current_var]] <- list(
    variable = current_var,
    variance_explained_total = var_explained_total,
    variance_explained_remaining = var_explained_remaining,
    remaining_variance_after = remaining_var,
    p_value = sig_test$`Pr(>F)`[1],
    F_statistic = sig_test_phenolog$F[1]
  )
  cat("Level:", current_var, "\n")
  cat("Variance explained (of total):", round(var_explained_total, 4), "\n")
  cat("Variance explained (of remaining):", round(var_explained_remaining, 4), "\n")
  cat("P-value:", sig_test$`Pr(>F)`[1], "\n")
  cat("Remaining variance:", round(remaining_var, 4), "\n\n")
}
sequential_df <- do.call(rbind, lapply(sequential_phenolog, data.frame))
print(sequential_df)
#this shows the same results as doing them independently except that species increased (15% instead of 12.4%)
#also the p values are odd, all showing the same (0.721)
###I think that's because I used sig_test instead of sig_test_phenolog in a couple of places- check

#restricted permutations to get meaningful p values, though with reduced statistical power
seq_phenolog_restr <- list()
for (i in 1:5) {
  current_var_rest <- phenolog_vars[i]
  if(i == 1) {
    formula_rest <- as.formula(paste("peak_dist_matrix_all ~", current_var_rest))
    dbrda_rest <- dbrda(formula_rest, data = phenolog_meta)
    sig_test_rest <- anova(dbrda_rest, permutations = 999)
  } else {
    conditioning_vars_rest <- paste(phenolog_vars[1:(i-1)], collapse = " + ")
    formula_rest <- as.formula(paste("peak_dist_matrix_all ~", current_var_rest, "+ Condition(", conditioning_vars_rest, ")"))
    dbrda_rest <- dbrda(formula_rest, data = phenolog_meta)
    if(i == 2) {
      perm_design <- how(within = Within(type = "free"),
                         blocks = phenolog_meta[[phenolog_vars[1]]])
    } else {
      perm_design <- how(within = Within(type = "free"),
                         blocks = phenolog_meta[[phenolog_vars[i-1]]])
    }
    sig_test_rest <- anova(dbrda_rest, permutations = perm_design)
  }
  var_explained_total_rest <- dbrda_rest$CCA$tot.chi / dbrda_rest$tot.chi
  seq_phenolog_restr[[current_var_rest]] <- list(
    variable = current_var_rest,
    variance_explained = var_explained_total_rest,
    p_value = sig_test_rest$`Pr(>F)`[1],
    F_statistic <- sig_test_rest$F[1]
  )
  cat("Level:", current_var_rest, "\n")
  cat("Variance explained:", round(var_explained_total_rest, 4), "\n")
  cat("P-value:", sig_test_rest$`Pr(>F)`[1], "\n\n")
}
#result: class explains 5.76% of variation, p-value .002; order explains 20.23%, p .005; 
#family explains 7.21%, p .05; genus explains 7.34%, p .055; species explains 15%, p .06
#these are percentages of total variance explained
#of remaining variance, it's 5.76%, 21.47%, 9.74%, 10.99%, and 25.24%
#remaining variance after is 94.23%, 74.00%, 66.79%, 59.45%, and 44.44%

plot(dbrda_phenolog, type = "n")
points(dbrda_phenolog, display = "sites", pch = 19, col = Order, cex = 0.8)
text(dbrda_phenolog, display = "sites", labels = rownames(rda_metadata), pos = 3, cex = 0.6, offset = 0.3)
text(dbrda_phenolog, display = family, col = "red", cex = 0.8)

#exporting dbrda for use in jmp
site_scores <- scores(dbrda, display = "sites")
site_scores_df <- data.frame(SampleID = rownames(site_scores), site_scores)
site_scores_full <- cbind(site_scores_df, tf_meta)
write.csv(site_scores_full, "dbrda_scores.csv", row.names = FALSE)

#exporting phylogeny dbrda for use in jmp
site_scores_phenolog <- scores(dbrda_phenolog, display = "sites")
site_scores_phenolog_df <- data.frame(SampleID = rownames(site_scores_phenolog), site_scores_phenolog)
site_scores_phenolog_full <- cbind(site_scores_phenolog_df, tf_meta)
write.csv(site_scores_phenolog_full, "dbrda_scores_phenolog.csv", row.names = FALSE)


#analysis of env variables on residuals after phylogeny
species_formula <- as.formula(paste("peak_dist_matrix_all ~", paste(phenolog_meta, collapse = " + ")))
dbrda_species <- dbrda(species_formula, data = phenolog_meta)
residual_after_phylo <- residuals(dbrda_species)
residual_dist <- as.dist(residual_after_phylo)
residual_matrix <- as.matrix(residual_dist)

#analysis of other factors after accounting for DOC
DOC_residual <- residuals(dbrda_DOC)
DOC_residual_dist <- as.dist(DOC_residual)
DOC_residual_matrix <- as.matrix(DOC_residual_dist)
#comparing to #2
DOC_residual_2 <- residuals(dbrda_DOC_2)
DOC_residual_matrix_2 <- as.matrix(DOC_residual_2)
write.csv(DOC_residual_matrix, "residuals_after_DOC.csv", row.names = FALSE)
write.csv(DOC_residual_matrix_2, "residuals_after_DOC_2.csv", row.names = FALSE)
#same output for both- 77x77 matrix of distance scores, they're identical
env_vars <- colnames(rda_metadata)
env_vars_after_DOC <- colnames(rda_metadata[, -(2)])
NA_DOC <- is.na(rda_metadata$DOC_uM)
env_meta_after_DOC <- rda_metadata[NA_DOC == FALSE, env_vars_after_DOC]
env_after_DOC <- list()
for(var in env_vars_after_DOC) {
  formula <- as.formula(paste("DOC_residual_matrix ~", var))
  single_dbrda <- dbrda(formula, data = env_meta_after_DOC)
  sig_test <- anova(single_dbrda, permutations = 999)
  var_explained_residual <- single_dbrda$CCA$tot.chi / single_dbrda$tot.chi
  original_unexplained <- dbrda_DOC$CCA$tot.chi / dbrda_DOC$tot.chi
  var_explained_total <- var_explained_residual * original_unexplained
  env_after_DOC[[var]] <- list(
    variable = var,
    variance_explained_of_residual = var_explained_residual,
    variance_explained_of_total = var_explained_total,
    p_value = sig_test$`Pr(>F)`[1],
    F_statistic = sig_test$F[1]
  )
  cat("Variable:", var, "\n")
  cat("Variance explained (of residual):", round(var_explained_residual, 4), "\n")
  cat("Variance explained (of total original):", round(var_explained_total, 4), "\n")
  cat("P-value:", sig_test$`Pr(>F)`[1], "\n\n")
}

#looking into distance matrix and residual properties
full_var <- sum(peak_dist_matrix_all^2) / (2 * nrow(peak_dist_matrix_all))
DOC_var <- sum(peak_dist_matrix_DOC^2) / (2 * nrow(peak_dist_matrix_DOC))
#variance is 11.04 for full matrix, 4.58 for DOC-only matrix

#trying partial db-rda instead
#conditions on DOC and avoids need to calculate residual variance
phylo_meta_with_DOC <- phenolog_meta[DOC_test, ]
for(var in env_vars_after_DOC) {
  DOC_condition <- paste(DOC_uM, collapse = " + ")
  formula <- as.formula(paste("peak_dist_matrix_DOC ~", var, "+ Condition(", DOC_condition, ")"))
  partial_dbrda_DOC <- dbrda(formula, data = env_meta_after_DOC)
  sig_test <- anova(partial_dbrda_DOC, permutations = 999)
  var_explained_additional <- partial_dbrda_DOC$CCA$tot.chi / partial_dbrda_DOC$tot.chi
  cat("Variable:", var, "\n")
  cat("Additional variance explained (beyond DOC):", round(var_explained_additional, 4), "\n")
  cat("P-value:", sig_test$`Pr(>F)`[1], "\n\n")
}
#not working, trying from scratch
#each level individually, conditioned on DOC
tf_meta_after_DOC <- tf_meta[DOC_test, ]
tf_meta_after_DOC$DOC_uM <- as.numeric(tf_meta_after_DOC$DOC_uM)
rownames(tf_meta_after_DOC) <- tf_samples_w_DOC
tf_meta_after_DOC[phenolog_vars] <- lapply(tf_meta_after_DOC[phenolog_vars], as.factor)
phylo_cond_on_DOC <- list()
for(var in phenolog_vars) {
  formula_DOC_cond <- as.formula(paste("peak_dist_matrix_DOC ~", var, "+ Condition(DOC_uM)"))
  partial_dbrda_DOC <- dbrda(formula_DOC_cond, data = tf_meta_after_DOC)
  sig_test_partial_DOC <- anova(partial_dbrda_DOC, permutations = 999)
  var_explained_partial <- partial_dbrda_DOC$CCA$tot.chi / partial_dbrda_DOC$tot.chi
  doc_variance <- partial_dbrda_DOC$pCCA$tot.chi / partial_dbrda_DOC$tot.chi
  phylo_cond_on_DOC[[var]] <- list(
    taxonomic_level = var,
    add_variance_beyond_DOC = var_explained_partial,
    DOC_variance_conditioned = doc_variance,
    p_value = sig_test$`Pr(>F)`[1],
    f_statistic = sig_test$F[1]
  )
  cat("Taxonomic level:", var, "\n")
  cat("Additional variance beyond DOC:", round(var_explained_partial, 4), "\n")
  cat("DOC variance conditioned out:", round(doc_variance, 4), "\n")
  cat("P-value:", sig_test$`Pr(>F)`[1], "\n\n")
}
phylo_partial_df <- do.call(rbind, lapply(phylo_cond_on_DOC, data.frame))
print(phylo_partial_df)
#conditioned out the 11% of variance covered by DOC
#additional 8.7%, 28.8%, 41.1%, 46.7%, and 58.4% of variance 
#explained by class, order, family, genus, and species respectively
#p values not meaningful again

#testing levels sequentially, conditioned on DOC
sequential_cond_on_DOC <- list()
for(i in 1:5) {
  current_var <- phenolog_vars[i]
  if(i == 1) {
    formula_seq_DOC <- as.formula(paste("peak_dist_matrix_DOC ~", current_var, "+ Condition(DOC_uM)"))
  } else {
    previous_vars <- paste(phenolog_vars[1:(i-1)], collapse = " + ")
    formula_seq_DOC <- as.formula(paste("peak_dist_matrix_DOC ~", current_var,
                                        "+ Condition(DOC_uM +", previous_vars, ")"))
  }
  partial_dbrda_seq <- dbrda(formula_seq_DOC, data = tf_meta_after_DOC)
  sig_test_cond_seq <- anova(partial_dbrda_seq, permutations = 999)
  var_explained_cond_seq <- partial_dbrda_seq$CCA$tot.chi / partial_dbrda_seq$tot.chi
  conditioned_variance <- partial_dbrda_seq$pCCA$tot.chi / partial_dbrda_seq$tot.chi
  sequential_cond_on_DOC[[current_var]] <- list(
    taxonomic_level = current_var,
    additional_variance = var_explained_cond_seq,
    conditioned_variance = conditioned_variance,
    p_value = sig_test_cond_seq$`Pr(>F)`[1]
  )
  cat("Level:", current_var, "\n")
  cat("Additional variance:", round(var_explained_cond_seq, 4), "\n")
  cat("P-value:", sig_test_cond_seq$`Pr(>F)`[1], "\n\n")
}
#sequentially, class explains 8.7% (p=.001), order 20.1% (p=.002),
#family 12.3% (p=.001), genus 5.6% (p=.168), species 13.6% (p=.035)
#again sequential percentages are just separating the components of looking at them individually
#except for species, which is 13.6% when 11.6% would be expected from individual values

#testing which axes are significant
anova(partial_dbrda_DOC, by = "axis", permutations = 999)
anova(partial_dbrda_seq, by = "axis", permutations = 999)
#neither is working, first because of negative eigenvalues
eigenvals(partial_dbrda_DOC, constrained = TRUE, unconstrained = TRUE)
#there are 32 axes, last one is negative so don't want it
scores_DOC_ind <- scores(partial_dbrda_DOC, display = "sites", choices = 1:31)
#axis significance testing
positive_axes_ind <- partial_dbrda_DOC
positive_axes_ind$CCA$eig <- partial_dbrda_DOC$CCA$eig[1:31]
positive_axes_ind$CCA$u <- partial_dbrda_DOC$CCA$u[1:31]
positive_axes_ind$CCA$v <- partial_dbrda_DOC$CCA$v[1:31]
positive_axes_ind$CCA$rank <- 31
anova(positive_axes_ind, by = "axis", permutations = 999)

#axis explanatory power testing
positive_eig_DOC_ind <- partial_dbrda_DOC$CCA$eig[partial_dbrda_DOC$CCA$eig > 0]
total_inertia_DOC_ind <- partial_dbrda_DOC$tot.chi
constrained_inertia_DOC_ind <- sum(positive_eig_DOC_ind)
axes_DOC_ind_var <- data.frame(
  Axis = paste0("dbRDA", 1:31),
  Eigenvalue = positive_eig_DOC_ind,
  Prop_Total_Variance = positive_eig_DOC_ind / total_inertia_DOC_ind,
  Prop_Constrained_Variance = positive_eig_DOC_ind / constrained_inertia_DOC_ind,
  Cumulative_Constrained = cumsum(positive_eig_DOC_ind) / constrained_inertia_DOC_ind,
  Important_5pct = (positive_eig_DOC_ind / constrained_inertia_DOC_ind) > 0.05,
  Important_1pct = (positive_eig_DOC_ind / constrained_inertia_DOC_ind) > 0.01
)
write.csv(axes_DOC_ind_var, "dbrda_axes_after_DOC.csv", row.names = FALSE)
#first 2 axes explain at least 10% of constrained variance each (28.8% and 17.3%)
#first 5 explain 5%, first 17 explain 1%

#exporting information for k means clustering
#adding weighting by sqrt of eigenvalue and proportion of variation explained




#alternative is to extract scores directly from dbrda object
#should give the same result
scores_DOC_ind_2 <- partial_dbrda_DOC$CCA$u[, 1:31] %*% diag(sqrt(positive_eig_DOC_seq))
colnames(scores_DOC_ind_2) <- paste0("dbRDA", 1:31)


scores_DOC_ind_df <- data.frame(Sample_ID = rownames(scores_DOC_ind), scores_DOC_ind)
scores_DOC_ind_export <- cbind(scores_DOC_ind_df, tf_meta_after_DOC)
write.csv(scores_DOC_ind_export, "scores_DOC_ind.csv", row.names = FALSE)

#now looking at sequential dbrda
eigenvals(partial_dbrda_seq, constrained = TRUE, unconstrained = TRUE)
#13 axes, all positive
scores_DOC_seq <- scores(partial_dbrda_seq, display = "sites")
scores_DOC_seq_df <- data.frame(Sample_ID = rownames(scores_DOC_seq), scores_DOC_seq)
scores_DOC_seq_export <- cbind(scores_DOC_seq_df, tf_meta_after_DOC)
write.csv(scores_DOC_seq_export, "scores_DOC_seq.csv", row.names = FALSE)
#these look very different
#independent looks like it hasn't removed variation by phylo level
#sequential looks like it has



#alternative way of accounting for DOC enrichment
#filter out samples non significantly enriched in DOC
precip_meta <- samp_meta[samp_meta$Sample_type == "Precipitation", ]
mean_precip_DOC <- mean(precip_meta$DOC_uM)
double_precip_DOC <- DOC_samples[DOC_samples$DOC_uM > (2*mean_precip_DOC), ]
double_precip_DOC_samples <- rownames(double_precip_DOC)
triple_precip_DOC <- DOC_samples[DOC_samples$DOC_uM > (3*mean_precip_DOC), ]
triple_precip_DOC_samples <- rownames(triple_precip_DOC)
#77 samples have DOC data, 74 are 2x precip, 63 are 3x precip
peak_dist_matrix_2xDOC <- peak_dist_matrix_all[double_precip_DOC_samples, double_precip_DOC_samples]
peak_dist_matrix_3xDOC <- peak_dist_matrix_all[triple_precip_DOC_samples, triple_precip_DOC_samples]

#phylo levels independently, 2x precip DOC minimum
phylo_meta_2xDOC <- phenolog_meta[double_precip_DOC_samples, ]
phylo_signif_2xDOC <- list()
for(var in phenolog_vars) {
  formula_double <- as.formula(paste("peak_dist_matrix_2xDOC ~", var))
  single_dbrda_double <- dbrda(formula_double, data = phylo_meta_2xDOC)
  sig_test_double <- anova(single_dbrda_double, permutations = 999)
  var_explained_double <- single_dbrda_double$CCA$tot.chi / single_dbrda_double$tot.chi
  phylo_signif_2xDOC[[var]] <- list(
    variable = var,
    variance_explained = var_explained_double,
    p_value = sig_test_double$`Pr(>F)`[1],
    F_statistic = sig_test_double$F[1],
    df = sig_test_double$Df[1]
  )
}
phylo_signif_2xDOC_results <- do.call(rbind, lapply(phylo_signif_2xDOC, data.frame))
phylo_signif_2xDOC_results <- phylo_signif_2xDOC_results[order(phylo_signif_2xDOC_results$variance_explained, decreasing = TRUE), ]
print(phylo_signif_2xDOC_results)

#permanova on phylogeny data- this shows the same results as dbrda above
for(var in phenolog_vars) {
  formula <- as.formula(paste("peak_dist_matrix_2xDOC ~", var))
  permanova_result <- adonis2(formula, data = phylo_meta_2xDOC, permutations = 999)
  print(paste("Variable:", var))
  print(permanova_result)
  cat("\n")
}

#phylo levels independently, 3x precip DOC minimum
phylo_meta_3xDOC <- phenolog_meta[triple_precip_DOC_samples, ]
phylo_signif_3xDOC <- list()
for(var in phenolog_vars) {
  formula_triple <- as.formula(paste("peak_dist_matrix_3xDOC ~", var))
  single_dbrda_triple <- dbrda(formula_triple, data = phylo_meta_3xDOC)
  sig_test_triple <- anova(single_dbrda_triple, permutations = 999)
  var_explained_triple <- single_dbrda_triple$CCA$tot.chi / single_dbrda_triple$tot.chi
  phylo_signif_3xDOC[[var]] <- list(
    variable = var,
    variance_explained = var_explained_triple,
    p_value = sig_test_triple$`Pr(>F)`[1],
    F_statistic = sig_test_triple$F[1],
    df = sig_test_triple$Df[1]
  )
}
phylo_signif_3xDOC_results <- do.call(rbind, lapply(phylo_signif_3xDOC, data.frame))
phylo_signif_3xDOC_results <- phylo_signif_3xDOC_results[order(phylo_signif_3xDOC_results$variance_explained, decreasing = TRUE), ]
print(phylo_signif_3xDOC_results)

#permanova on phylogeny data- this shows the same results as dbrda above
for(var in phenolog_vars) {
  formula <- as.formula(paste("peak_dist_matrix_3xDOC ~", var))
  permanova_result <- adonis2(formula, data = phylo_meta_3xDOC, permutations = 999)
  print(paste("Variable:", var))
  print(permanova_result)
  cat("\n")
}

#min 2x precip DOC samples: class explains 15.0% of variation,
#order 36.4%, family 51.3%, genus 58.0%, species 69.9%
#all significant at .001 level
#min 3x precip DOC samples: class explains 15.3% of variation,
#order 37.2%, family 53.7%, genus 61.0%, species 73.8%
#all significant at .001 level

#now repeat but do it sequentially
#phylo levels sequentially, 2x precip DOC minimum
sequential_phylo_2xDOC <- list()
remaining_var <- 1.0
for(i in 1:5) {
  current_var <- phenolog_vars[i]
  if(i == 1) {
    formula_double_seq <- as.formula(paste("peak_dist_matrix_2xDOC ~", current_var))
    dbrda_double_seq <- dbrda(formula_double_seq, data = phylo_meta_2xDOC)
  } else {
    conditioning_vars <- paste(phenolog_vars[1:(i-1)], collapse = " + ")
    formula_double_seq <- as.formula(paste("peak_dist_matrix_2xDOC ~", current_var, "+ Condition(", conditioning_vars, ")"))
    dbrda_double_seq <- dbrda(formula_double_seq, data = phylo_meta_2xDOC)
  }
  if(i == 1) {
    var_explained_total <- dbrda_double_seq$CCA$tot.chi / dbrda_double_seq$tot.chi
    var_explained_remaining <- var_explained_total
  } else {
    var_explained_total <- dbrda_double_seq$CCA$tot.chi / dbrda_double_seq$tot.chi
    var_explained_remaining <- var_explained_total / remaining_var
  }
  remaining_var <- remaining_var - var_explained_total
  sig_test_double_seq <- anova(dbrda_double_seq, permutations = 999)
  sequential_phylo_2xDOC[[current_var]] <- list(
    variable = current_var,
    variance_explained_total = var_explained_total,
    variance_explained_remaining = var_explained_remaining,
    remaining_variance_after = remaining_var,
    p_value = sig_test_double_seq$`Pr(>F)`[1],
    F_statistic = sig_test_double_seq$F[1]
  )
  cat("Level:", current_var, "\n")
  cat("Variance explained (of total):", round(var_explained_total, 4), "\n")
  cat("Variance explained (of remaining):", round(var_explained_remaining, 4), "\n")
  cat("P-value:", sig_test_double_seq$`Pr(>F)`[1], "\n")
  cat("Remaining variance:", round(remaining_var, 4), "\n\n")
}
sequential_2xDOC_df <- do.call(rbind, lapply(sequential_phylo_2xDOC, data.frame))
print(sequential_2xDOC_df)

#phylo levels sequentially, 3x precip DOC minimum
sequential_phylo_3xDOC <- list()
remaining_var <- 1.0
for(i in 1:5) {
  current_var <- phenolog_vars[i]
  if(i == 1) {
    formula_triple_seq <- as.formula(paste("peak_dist_matrix_3xDOC ~", current_var))
    dbrda_triple_seq <- dbrda(formula_triple_seq, data = phylo_meta_3xDOC)
  } else {
    conditioning_vars <- paste(phenolog_vars[1:(i-1)], collapse = " + ")
    formula_triple_seq <- as.formula(paste("peak_dist_matrix_3xDOC ~", current_var, "+ Condition(", conditioning_vars, ")"))
    dbrda_triple_seq <- dbrda(formula_triple_seq, data = phylo_meta_3xDOC)
  }
  if(i == 1) {
    var_explained_total <- dbrda_triple_seq$CCA$tot.chi / dbrda_triple_seq$tot.chi
    var_explained_remaining <- var_explained_total
  } else {
    var_explained_total <- dbrda_triple_seq$CCA$tot.chi / dbrda_triple_seq$tot.chi
    var_explained_remaining <- var_explained_total / remaining_var
  }
  remaining_var <- remaining_var - var_explained_total
  sig_test_triple_seq <- anova(dbrda_triple_seq, permutations = 999)
  sequential_phylo_3xDOC[[current_var]] <- list(
    variable = current_var,
    variance_explained_total = var_explained_total,
    variance_explained_remaining = var_explained_remaining,
    remaining_variance_after = remaining_var,
    p_value = sig_test_triple_seq$`Pr(>F)`[1],
    F_statistic = sig_test_triple_seq$F[1]
  )
  cat("Level:", current_var, "\n")
  cat("Variance explained (of total):", round(var_explained_total, 4), "\n")
  cat("Variance explained (of remaining):", round(var_explained_remaining, 4), "\n")
  cat("P-value:", sig_test_triple_seq$`Pr(>F)`[1], "\n")
  cat("Remaining variance:", round(remaining_var, 4), "\n\n")
}
sequential_3xDOC_df <- do.call(rbind, lapply(sequential_phylo_3xDOC, data.frame))
print(sequential_3xDOC_df)
plot(sequential_cond_on_DOC)

#min 2x precip DOC samples: class explains 15.0% of total variation
#order 21.4%, family 14.9%, genus 6.7%, species 13.2%
#of remaining variation: order 25.2%, family 23.4%, genus 13.7%, species 31.5%
#class, order, and family significant at .001 level
#genus and species not significant because of limited sample size

#min 3x precip DOC samples: class explains 15.3% of total variation
#order 21.9%, family 16.5%, genus 7.2%, species 12.8%
#of remaining variation: order 25.8%, family 26.3%, genus 15.7%, species 32.9%
#class significant at .001 level, order .003, family .001
#genus and species not significant because of limited sample size

#silhouette analysis to look at hierarchical clustering
#forcing taxonomic groups to be clusters
#on all samples
clust_matr_all <- read.csv(file = 'data/clustering_matrix_all.csv')

clustering_class_all <- as.numeric(as.factor(tf_meta$Class))
silh_class_all <- silhouette(clustering_class_all, clust_matr_all)

clustering_order_all <- as.numeric(as.factor(tf_meta$Order))
silh_order_all <- silhouette(clustering_order_all, clust_matr_all)

clustering_family_all <- as.numeric(as.factor(tf_meta$Family))
silh_family_all <- silhouette(clustering_family_all, clust_matr_all)

clustering_genus_all <- as.numeric(as.factor(tf_meta$Genus))
silh_genus_all <- silhouette(clustering_genus_all, clust_matr_all)

clustering_species_all <- as.numeric(as.factor(tf_meta$Species))
silh_species_all <- silhouette(clustering_species_all, clust_matr_all)

silh_output_all <- cbind(silh_class_all, silh_order_all, silh_family_all, silh_genus_all, silh_species_all)
write.csv(silh_output_all, "silh_output_all.csv")

#on samples with >3x precip DOC
clust_matr_3xDOC <- read.csv(file = 'data/clustering_matrix_3xDOC.csv')
tf_meta_rownames <- tf_meta
rownames(tf_meta_rownames) <- tf_meta$Sample_ID
tf_meta_3xDOC <- tf_meta_rownames[triple_precip_DOC_samples, ]

clustering_class_3xDOC <- as.numeric(as.factor(tf_meta_3xDOC$Class))
silh_class_3xDOC <- silhouette(clustering_class_3xDOC, clust_matr_3xDOC)

clustering_order_3xDOC <- as.numeric(as.factor(tf_meta_3xDOC$Order))
silh_order_3xDOC <- silhouette(clustering_order_3xDOC, clust_matr_3xDOC)

clustering_family_3xDOC <- as.numeric(as.factor(tf_meta_3xDOC$Family))
silh_family_3xDOC <- silhouette(clustering_family_3xDOC, clust_matr_3xDOC)

clustering_genus_3xDOC <- as.numeric(as.factor(tf_meta_3xDOC$Genus))
silh_genus_3xDOC <- silhouette(clustering_genus_3xDOC, clust_matr_3xDOC)

clustering_species_3xDOC <- as.numeric(as.factor(tf_meta_3xDOC$Species))
silh_species_3xDOC <- silhouette(clustering_species_3xDOC, clust_matr_3xDOC)

silh_output_3xDOC <- cbind(silh_class_3xDOC, silh_order_3xDOC, silh_family_3xDOC, silh_genus_3xDOC, silh_species_3xDOC)
write.csv(silh_output_3xDOC, "silh_output_3xDOC.csv")





#mantel
large_phylo_dist <- as.dist(large_phylo_matrix)
#do 999 permutations for speed, 9999 for final published result
mantel <- mantel(large_phylo_dist, peak_dist,
                 method = "spearman", permutations = 999)
#r=0.01525, significance = 0.37-- not significant
#try again on updated matrices
mantel_phylo <- as.dist(phylo_dist_matrix_final_sample)
mantel_output <- mantel(mantel_phylo, peak_dist_all,
                        method = "spearman", permutations = 999)
#r=-.007334, significance = 0.579-- not significant


#PCoA axis creation
phylo_pcoa <- pcoa(large_phylo_matrix)
plot(phylo_pcoa, x = Axis.1, y = Axis.2)
#having problems, looking at eigenvalues
#eig <- eigen(scaled_large_phylo, symmetric=TRUE)$values
#print("First 10 eigenvalues:")
#print(eig[1:min(10, length(eig))])
#print(paste("Number of positive eigenvalues:", sum(eig > 0)))
#eigenvalues are too small to process correctly
#can just multiply whole matrix by a constant to scale it up
#max in matrix is 1.76*10^-310, so I'll multiply all by 10^309
#power of 309 is too big, R treats it as infinity
scale <- 10 ^ 103
scaled_large_phylo <- large_phylo_matrix * scale * scale * scale
phylo_pcoa <- pcoa(scaled_large_phylo)
#looks like there are negative eigenvalues here, may want to go back and apply a correction
phylo_axes <- phylo_pcoa$vectors
cumulative_var <- cumsum(phylo_pcoa$values$Relative_eig)
num_axes <- which(cumulative_var >= 0.8)[1]
phylo_axes_reduced <- phylo_axes[, 1:num_axes]
#80% of cumulative variance explained by 3 axes, so I'll use 3

#cca
nonhyb_data <- norm_intens_t[-c(48, 49, 50, 82, 83, 84),]
cca <- cca(nonhyb_data ~ ., data = as.data.frame(phylo_axes_reduced))
summary(cca)
anova(cca)
#p value .001
anova(cca, by = "axis")
#vector too large to allocate (7.7Mb)
cca_r2 <- RsquareAdj(cca)$r.squared
cca_r2_adj <- RsquareAdj(cca)$adj.r.squared
print(paste("Phylogeny explains", round(cca_r2_adj * 100, 2),
            "% of variation in FTICRMS data (adjusted r2)"))
#[1] "Phylogeny explains 3.3 % of variation in FTICRMS data (adjusted r2)"
plot(cca, scaling = 2)
points(cca, display = "sites", pch = 16, scaling = 2)
env_data <- as.data.frame(samp_meta_nonhyb)
env_arrows <- envfit(cca, env_data, na.rm = TRUE, perm = 999)

#rda
rda <- rda(nonhyb_data ~ ., data = as.data.frame(phylo_axes_reduced))
summary(rda)
anova(rda)
#p value .002
anova(rda, by = "axis")
#vector too large again
rda_r2 <- RsquareAdj(rda)$r.squared
rda_r2_adj <- RsquareAdj(rda)$adj.r.squared
print(paste("Phylogeny explains", round(rda_r2_adj * 100, 2),
            "% of variation in FTICRMS data (adjusted r2)"))
#[1] "Phylogeny explains 3.16 % of variation in FTICRMS data (adjusted r2)"
plot(rda, scaling = 2)

#db-rda
db_rda <- capscale(peak_dist_matrix ~ ., data = as.data.frame(phylo_axes_reduced))
summary(db_rda)
anova(db_rda)
#p value .004
anova(db_rda, by = "axis")
#this one ran successfully- axis 1 p=.02, axis 2 p=.235, axis 3 p=.621
db_rda_r2 <- RsquareAdj(db_rda)$r.squared
db_rda_r2_adj <- RsquareAdj(db_rda)$adj.r.squared
print(paste("Phylogeny explains", round(db_rda_r2_adj * 100, 2),
            "% of variation in FTICRMS data (adjusted r2)"))
#[1] "Phylogeny explains 5.23 % of variation in FTICRMS data (adjusted r2)"

plot(db_rda, axes = CAP1 & CAP3, scaling = 2)


#first make later things easier by trimming samp meta to desired rows
samp_meta_nonhyb <- samp_meta[-c(48, 49, 50, 54, 55, 56, 57, 58, 87, 88, 89),]
factor <- as.factor(samp_meta_nonhyb$PhyloName)
names(factor) <- samp_meta_nonhyb$Sample_ID
dbrda <- vegan::capscale(peak_dist_matrix ~ factor)
anova(dbrda)
#p<.001, so there is a significant species effect at least
#next create PCoA axes from phylogenetic distances
phylo_pcoa_large <- ape::pcoa(large_phylo_matrix)
phylo_pcoa <- ape::pcoa(phylo_dist_matrix, correction = "cailliez")
if (inherits(phylo_pcoa, "try-error")) {
  phylo_pcoa <- pcoa(phylo_dist_matrix, correction = "lingoes")
}
phylo_pcoa <- ape::pcoa(phylo_dist_matrix, correction = "cailliez")
if (inherits(phylo_pcoa, "try-error")) {
  cmd_result <- cmdscale(phylo_dist_matrix, k = min(10, 33), eig = TRUE)
}




#save data as .Rds (object) file
#this way they can be used in different scripts using readRDS() command
saveRDS(intensities,'data/objects/intensities.Rds')
saveRDS(norm_intensities, 'data/objects/norm_intensities.Rds')
saveRDS(norm_transp_intensities, 'data/objects/norm_transp_intensities.Rds')
