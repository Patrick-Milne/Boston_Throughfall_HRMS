# Packages
library(tidyverse)
library(readr)
library(here)
library(lubridate)
library(vegan)
library(beepr)
library(openxlsx)
library(ape)
library(rotl)
library(picante)
library(taxize)
library(ggplot2)
library(ggtree)
library(ggvegan)
library(corrplot)
library(permute)
library(cluster)
library(plotrix)
library(dplyr)
library(data.tree)
library(DiagrammeR)
here <- here()


# sampling data
samp_meta <- read_csv(file = 'data/sample_metadata.csv')
#some NA values from missing DOC data and parameters that don't apply to rainfall samples
#making this data numeric instead of character
samp_meta$DOC_uM <- as.numeric(samp_meta$DOC_uM)
samp_meta$DOC_mgL <- as.numeric(samp_meta$DOC_mgL)
samp_meta$DOCflux_mg <- as.numeric(samp_meta$DOCflux_mg)
samp_meta$PPL_extraction_eff <- as.numeric(samp_meta$PPL_extraction_eff)
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
samp_meta$Tf_depth_mm <- samp_meta$Tf_volume_mL*1000/29460
#29460 is bucket opening size in mm^2
samp_meta$DOCflux_mg_m2 <- samp_meta$DOCflux_mg/.02946
samp_meta$TDNflux_mg_m2 <- samp_meta$TDNflux_mg/.02946
precip_bucket_avg_mm <- mean(samp_meta$Tf_depth_mm[samp_meta$Sample_type == "Precipitation"])
precip_bucket_avg_mL <- mean(samp_meta$Tf_volume_mL[samp_meta$Sample_type == "Precipitation"])
samp_meta$Tf_depth_pct_buckets <- samp_meta$Tf_depth_mm*100 / precip_bucket_avg_mm
samp_meta$Tf_depth_pct_gauge <- samp_meta$Tf_depth_mm*100 / 28.956
#rain gauge on NU campus measured total 1.14 in or 28.956 mm of rain
mean_precip_DOC_flux_m2 <- mean(samp_meta$DOCflux_mg_m2[samp_meta$Sample_type == "Precipitation"])
samp_meta$DOC_enrich <- samp_meta$DOCflux_mg_m2 / mean_precip_DOC_flux_m2
mean_precip_TDN_flux_m2 <- mean(samp_meta$TDNflux_mg_m2[samp_meta$Sample_type == "Precipitation"])
samp_meta$TDN_enrich <- samp_meta$TDNflux_mg_m2 / mean_precip_TDN_flux_m2
samp_meta$DOC.TDN_ratio <- samp_meta$DOC_uM / samp_meta$TDN_uM


all_samples <- samp_meta %>% pull(Sample_ID)
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

tf_meta <- samp_meta[samp_meta$Sample_type == "Throughfall", ]
precip_meta <- samp_meta[samp_meta$Sample_type == "Precipitation", ]

mean(precip_meta$Tf_depth_mm)
std.error(precip_meta$Tf_depth_mm)
min(precip_meta$Tf_depth_mm)
max(precip_meta$Tf_depth_mm)

mean(precip_meta$DOC_mgL)
std.error(precip_meta$DOC_mgL)
mean(precip_meta$TDN_mgL)
std.error(precip_meta$TDN_mgL)

mean(precip_meta$DOCflux_mg_m2)
std.error(precip_meta$DOCflux_mg_m2)
mean(precip_meta$TDNflux_mg_m2)
std.error(precip_meta$TDNflux_mg_m2)

mean(precip_meta$PPL_extraction_eff)
std.error(precip_meta$PPL_extraction_eff)
min(precip_meta$PPL_extraction_eff)
max(precip_meta$PPL_extraction_eff)

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

mean(tf_meta$DOCflux_mg_m2, na.rm = TRUE)
std.error(tf_meta$DOCflux_mg_m2, na.rm = TRUE)
min(tf_meta$DOCflux_mg_m2, na.rm = TRUE)
max(tf_meta$DOCflux_mg_m2, na.rm = TRUE)

mean(tf_meta$TDNflux_mg_m2)
std.error(tf_meta$TDNflux_mg_m2)
min(tf_meta$TDNflux_mg_m2)
max(tf_meta$TDNflux_mg_m2)

mean(tf_meta$DOC_enrich, na.rm = TRUE)
std.error(tf_meta$DOC_enrich, na.rm = TRUE)
min(tf_meta$DOC_enrich, na.rm = TRUE)
max(tf_meta$DOC_enrich, na.rm = TRUE)

mean(tf_meta$TDN_enrich)
std.error(tf_meta$TDN_enrich)
min(tf_meta$TDN_enrich)
max(tf_meta$TDN_enrich)

mean(tf_meta$DOC.TDN_ratio, na.rm = TRUE)
std.error(tf_meta$DOC.TDN_ratio, na.rm = TRUE)
min(tf_meta$DOC.TDN_ratio, na.rm = TRUE)
max(tf_meta$DOC.TDN_ratio, na.rm = TRUE)

mean(tf_meta$PPL_extraction_eff, na.rm = TRUE)
std.error(tf_meta$PPL_extraction_eff, na.rm = TRUE)
min(tf_meta$PPL_extraction_eff, na.rm = TRUE)
max(tf_meta$PPL_extraction_eff, na.rm = TRUE)

tf_meta %>%
  group_by(Species) %>%
  summarise_at(vars(Tf_depth_mm), list(name = mean, std.error))|> print(n=36)

tf_meta %>%
  group_by(Species) %>%
  summarise_at(vars(Tf_depth_pct_buckets), list(name = mean, std.error))|> print(n=36)

mean(tf_meta$Tf_depth_mm[tf_meta$Class == "Pinopsida"])
std.error(tf_meta$Tf_depth_mm[tf_meta$Class == "Pinopsida"])
mean(tf_meta$Tf_depth_mm[tf_meta$Class == "Dicotyledoneae"])
std.error(tf_meta$Tf_depth_mm[tf_meta$Class == "Dicotyledoneae"])

mean(tf_meta$Tf_depth_pct_buckets[tf_meta$Class == "Pinopsida"])
std.error(tf_meta$Tf_depth_pct_buckets[tf_meta$Class == "Pinopsida"])
mean(tf_meta$Tf_depth_pct_buckets[tf_meta$Class == "Dicotyledoneae"])
std.error(tf_meta$Tf_depth_pct_buckets[tf_meta$Class == "Dicotyledoneae"])

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
mean(tf_meta$DOCflux_mg_m2[tf_meta$Order == "Magnoliales"], na.rm = TRUE)
std.error(tf_meta$DOCflux_mg_m2[tf_meta$Order == "Magnoliales"], na.rm = TRUE)
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

set.seed(0)
t.test(tf_meta$DBH_cm[tf_meta$Class == "Pinopsida"],
       tf_meta$DBH_cm[tf_meta$Class != "Pinopsida"])

set.seed(0)
t.test(tf_meta$DBH_cm[tf_meta$Order == "Magnoliales"],
       tf_meta$DBH_cm[tf_meta$Order != "Magnoliales"])

set.seed(0)
t.test(tf_meta$DBH_cm[tf_meta$Order == "Pinales"],
       tf_meta$DBH_cm[tf_meta$Order == "Magnoliales"])

set.seed(0)
t.test(tf_meta$DOC.TDN_ratio[tf_meta$Class == "Pinopsida"],
       tf_meta$DOC.TDN_ratio[tf_meta$Class == "Dicotyledoneae"])

mean(tf_meta$DOC.TDN_ratio[tf_meta$Class == "Pinopsida"], na.rm = TRUE)
std.error(tf_meta$DOC.TDN_ratio[tf_meta$Class == "Pinopsida"], na.rm = TRUE)
mean(tf_meta$DOC.TDN_ratio[tf_meta$Class == "Dicotyledoneae"], na.rm = TRUE)
std.error(tf_meta$DOC.TDN_ratio[tf_meta$Class == "Dicotyledoneae"], na.rm = TRUE)


cor(tf_meta$BLDG_dist, tf_meta$Tf_depth_mm)
cor.test(tf_meta$BLDG_dist, tf_meta$Tf_depth_mm)
plot(x = tf_meta$BLDG_dist, y = tf_meta$Tf_depth_mm)

#Figure S4: dilution curve plot
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

plot(x = tf_meta$Tf_depth_pct_buckets, y = tf_meta$DOC_mgL,
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
                              expression(paste("p<0.0001"))),
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

plot(x = tf_meta$Tf_depth_pct_buckets, y = tf_meta$TDN_mgL,
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
                              expression(paste("p<0.0001"))),
       inset = c(-0.4, 0), y.intersp = 0.4,
       col = c("blue", NA, NA), lwd = c(1.5, NA, NA), bty = "n", cex = 0.75)
text(x=5, y=40, labels = "b", cex = 2)
#looks good


#then both on one plot
#start by plotting DOC
pdf("dilution_plot.pdf", width = 8, height = 8)
par(mfrow = c(1, 1))
par(mar = c(7, 5, 3, 5))
plot(x = tf_meta$Tf_depth_pct_buckets, y = tf_meta$DOC_mgL,
     pch = 1, col = "blue",
     xlim = c(0, 120), ylim = c(0, 140), xaxs = "i", yaxs = "i",
     xlab = "",
     ylab = "")
mtext(expression(paste("DOC concentration (mg-C L"^"-1",")")), 
      side = 2, line = 3, col = "blue")
pct_axis <- axis(1, labels = FALSE)
max_depth_mm <- 120 * precip_bucket_avg_mm / 100
depth_axis_labels <- seq(0, ceiling(max_depth_mm / 10) * 10, by = 10)
depth_axis <- depth_axis_labels * 100 / precip_bucket_avg_mm
axis(side = 1, at = depth_axis, labels = depth_axis_labels, line = 3.1)
mtext("Throughfall depth (% of rainfall)", side = 1, line = 1.9)
mtext("Throughfall depth (mm)", side = 1, line = 5)
axis(side = 2, col = "blue", col.axis = "blue")
lines(DOC_smooth_x, DOC_smooth_y, col = "blue", lwd = 1.5)
#TDN over it now
par(new = TRUE)
plot(x = tf_meta$Tf_depth_pct_buckets, y = tf_meta$TDN_mgL,
     pch = 2, col = "red3",
     axes = FALSE, xlab = "", ylab = "",
     xlim = c(0, 120), ylim = c(0, 42), xaxs = "i", yaxs = "i")
lines(TDN_smooth_x, TDN_smooth_y, col = "red3", lwd = 1.5)
axis(side = 4, col = "red3", col.axis = "red3")
segments(x0 = par("usr")[2], y0 = 0, y1 = 42, col = "red3", xpd = TRUE)
mtext(expression(paste("TDN concentration (mg-N L"^"-1",")")),
      side = 4, line = 3, col = "red3")
legend("topright",
       legend = c(expression(paste("[DOC] = 67.15 * e"^"-0.036x","")),
                  expression(paste("r"^"2"," = 0.217; p<0.0001")),
                  expression(paste("[TDN] = 16.81 * e"^"-0.055x","")),
                  expression(paste("r"^"2"," = 0.218; p<0.0001"))),
       inset = c(0.01, 0.01), y.intersp = c(1.9, 1.5, 1.5, 1.46),
       pch = c(1, NA, 2, NA), col = c("blue", NA, "red3", NA),
       lwd = 1.5, bty = "n")
dev.off()
#this is final dilution curve plot, Figure S4


#FTICRMS data
FTICRMS_raw <- read.csv(file = 'data/ICBM_OCEAN_output.csv')
#add columns with desired info
#condensed aromatics
FTICRMS_raw$Condensed.aromatic <- NA
FTICRMS_raw$Condensed.aromatic[FTICRMS_raw$AI.mod >= 0.67] <- 1
FTICRMS_raw$Condensed.aromatic[FTICRMS_raw$AI.mod < 0.67] <- 0

FTICRMS_raw$Condensed.aromatic.O_rich <- NA
FTICRMS_raw$Condensed.aromatic.O_rich[FTICRMS_raw$Condensed.aromatic == 1 & FTICRMS_raw$Aromatic.O_rich == 1] <- 1
FTICRMS_raw$Condensed.aromatic.O_rich[FTICRMS_raw$Condensed.aromatic != 1 | FTICRMS_raw$Aromatic.O_rich != 1] <- 0

FTICRMS_raw$Condensed.aromatic.O_poor <- NA
FTICRMS_raw$Condensed.aromatic.O_poor[FTICRMS_raw$Condensed.aromatic == 1 & FTICRMS_raw$Aromatic.O_poor == 1] <- 1
FTICRMS_raw$Condensed.aromatic.O_poor[FTICRMS_raw$Condensed.aromatic != 1 | FTICRMS_raw$Aromatic.O_poor != 1] <- 0

#remove condensed from the rest of the aromatics
FTICRMS_raw$Uncondensed.aromatic <- NA
FTICRMS_raw$Uncondensed.aromatic[FTICRMS_raw$Condensed.aromatic == 0 & FTICRMS_raw$Aromatic == 1] <- 1
FTICRMS_raw$Uncondensed.aromatic[FTICRMS_raw$Condensed.aromatic != 0 | FTICRMS_raw$Aromatic != 1] <- 0

FTICRMS_raw$Uncondensed.aromatic.O_rich <- NA
FTICRMS_raw$Uncondensed.aromatic.O_rich[FTICRMS_raw$Condensed.aromatic.O_rich == 0 & FTICRMS_raw$Aromatic.O_rich == 1] <- 1
FTICRMS_raw$Uncondensed.aromatic.O_rich[FTICRMS_raw$Condensed.aromatic.O_rich != 0 | FTICRMS_raw$Aromatic.O_rich != 1] <- 0

FTICRMS_raw$Uncondensed.aromatic.O_poor <- NA
FTICRMS_raw$Uncondensed.aromatic.O_poor[FTICRMS_raw$Condensed.aromatic.O_poor == 0 & FTICRMS_raw$Aromatic.O_poor == 1] <- 1
FTICRMS_raw$Uncondensed.aromatic.O_poor[FTICRMS_raw$Condensed.aromatic.O_poor != 0 | FTICRMS_raw$Aromatic.O_poor != 1] <- 0

#groupings by elements present
FTICRMS_raw <- FTICRMS_raw %>%
  mutate(non_CHO_atoms = rowSums(select(., c(N, S, P))))

FTICRMS_raw$CHO <- NA
FTICRMS_raw$CHO[FTICRMS_raw$non_CHO_atoms == 0] <- 1
FTICRMS_raw$CHO[FTICRMS_raw$non_CHO_atoms >0] <- 0

FTICRMS_raw$CHON <- NA
FTICRMS_raw$CHON[FTICRMS_raw$N > 0 & FTICRMS_raw$S == 0 & FTICRMS_raw$P == 0] <- 1
FTICRMS_raw$CHON[FTICRMS_raw$N == 0 | FTICRMS_raw$S > 0 | FTICRMS_raw$P > 0] <- 0

FTICRMS_raw$CHOS <- NA
FTICRMS_raw$CHOS[FTICRMS_raw$S > 0 & FTICRMS_raw$N == 0 & FTICRMS_raw$P == 0] <- 1
FTICRMS_raw$CHOS[FTICRMS_raw$S == 0 | FTICRMS_raw$N > 0 | FTICRMS_raw$P > 0] <- 0

FTICRMS_raw$CHOP <- NA
FTICRMS_raw$CHOP[FTICRMS_raw$P > 0 & FTICRMS_raw$S == 0 & FTICRMS_raw$N == 0] <- 1
FTICRMS_raw$CHOP[FTICRMS_raw$P == 0 | FTICRMS_raw$S > 0 | FTICRMS_raw$N > 0] <- 0

FTICRMS_raw$CHONS <- NA
FTICRMS_raw$CHONS[FTICRMS_raw$N > 0 & FTICRMS_raw$S > 0 & FTICRMS_raw$P == 0] <- 1
FTICRMS_raw$CHONS[FTICRMS_raw$N == 0 | FTICRMS_raw$S == 0 | FTICRMS_raw$P > 0] <- 0
#looks good

#sugar-like and peptide-like
FTICRMS_raw$sugar.like <- NA
FTICRMS_raw$sugar.like[FTICRMS_raw$O.C >= 0.9] <- 1
FTICRMS_raw$sugar.like[FTICRMS_raw$O.C < 0.9] <- 0

FTICRMS_raw$peptide.like <- NA
FTICRMS_raw$peptide.like[FTICRMS_raw$H.C >= 1.5 & FTICRMS_raw$H.C < 2 &
                           FTICRMS_raw$O.C < 0.9 & FTICRMS_raw$N > 0] <- 1
FTICRMS_raw$peptide.like[FTICRMS_raw$H.C < 1.5 | FTICRMS_raw$H.C >= 2 |
                           FTICRMS_raw$O.C >= 0.9 | FTICRMS_raw$N == 0] <- 0

#zero intensity was read in as NA- replace with 0, just in sample columns
FTICRMS_no_NA <- FTICRMS_raw %>%
  mutate(across(57:171, ~replace_na(., 0)))

#remove peaks detected in blanks- potential contaminants
FTICRMS_no_NA <- FTICRMS_no_NA %>%
  mutate(blank_detect = rowSums(select(., c(B1, B2, B3, B4, B5, IB1, IB2, IB3, IB4, IB5))))
FTICRMS_no_blank_peaks <- FTICRMS_no_NA[FTICRMS_no_NA$blank_detect == 0,]
#looks good, 611 peaks removed, 17094 peaks remain

#remove peaks with non-dominant isotopes- duplicates
FTICRMS_no_blank_peaks <- FTICRMS_no_blank_peaks %>%
  mutate(isotopes = rowSums(select(., c(C13, O18, N15, S34))))
FTICRMS_no_isotopes <- FTICRMS_no_blank_peaks[FTICRMS_no_blank_peaks$isotopes == 0,]
#looks good, 4644 peaks removed, 12450 peaks remain

#remove peaks with O=0- likely misassignments, would be insoluble
FTICRMS_positive_O <- FTICRMS_no_isotopes[FTICRMS_no_isotopes$O > 0,]
#looks good, 48 peaks removed, 12402 peaks remain

#remove peaks w/ systematic misassignment 1
FTICRMS_positive_O$misassignment_1 <- NA
FTICRMS_positive_O$misassignment_1[FTICRMS_positive_O$mz > 600 & FTICRMS_positive_O$N >= 2] <- 1
FTICRMS_positive_O$misassignment_1[FTICRMS_positive_O$mz <= 600 | FTICRMS_positive_O$N < 2] <- 0
FTICRMS_no_misassignment_1 <- FTICRMS_positive_O[FTICRMS_positive_O$misassignment_1 == 0,]
#looks good, 1203 peaks removed, 11199 peaks remain

#remove peaks w/ systematic misassignment 2
FTICRMS_no_misassignment_1$misassignment_2 <- NA
FTICRMS_no_misassignment_1$misassignment_2[FTICRMS_no_misassignment_1$N >= 1 & FTICRMS_no_misassignment_1$S >= 1 &
                                             FTICRMS_no_misassignment_1$H.C < 0.75 & FTICRMS_no_misassignment_1$O.C < 0.25] <- 1
FTICRMS_no_misassignment_1$misassignment_2[FTICRMS_no_misassignment_1$N < 1 | FTICRMS_no_misassignment_1$S < 1 |
                                             FTICRMS_no_misassignment_1$H.C >= 0.75 | FTICRMS_no_misassignment_1$O.C >= 0.25] <- 0
FTICRMS_no_misassignment_2 <- FTICRMS_no_misassignment_1[FTICRMS_no_misassignment_1$misassignment_2 == 0,]
#looks good, 72 peaks removed, 11127 peaks remain

#remove peaks w/ 4 N atoms- likely misassignment
FTICRMS_no_4N <- FTICRMS_no_misassignment_2[FTICRMS_no_misassignment_2$N < 4,]
#looks good, 66 peaks removed, 11061 peaks remain

#remove peaks w/ P and another non-CHO atom
FTICRMS_no_4N$P_and_another <- NA
FTICRMS_no_4N$P_and_another[FTICRMS_no_4N$P == 1 & FTICRMS_no_4N$non_CHO > 1] <- 1
FTICRMS_no_4N$P_and_another[FTICRMS_no_4N$P != 1 | FTICRMS_no_4N$non_CHO <= 1] <- 0
FTICRMS_no_P_and_another <- FTICRMS_no_4N[FTICRMS_no_4N$P_and_another == 0,]
#looks good, 123 peaks removed, 10938 peaks remain

#remove peaks w/ alternative formula provided
FTICRMS_no_alt_form <- FTICRMS_no_P_and_another[is.na(FTICRMS_no_P_and_another$alternative_formula),]
#looks good, 71 peaks removed, 10867 peaks remain

#remove peaks w/ duplicate masses
#either mz (one peak/measured mass was assigned two different formulae)
#or reference (two different peaks/measured masses were assigned the same formula)
FTICRMS_no_alt_form$dup.mass <- 0
FTICRMS_no_alt_form$dup.mass[duplicated(FTICRMS_no_alt_form$reference) | duplicated(FTICRMS_no_alt_form$reference, fromLast = TRUE) |
                               duplicated(FTICRMS_no_alt_form$mz) | duplicated(FTICRMS_no_alt_form$mz, fromLast = TRUE)] <- 1
FTICRMS_no_dupes <- FTICRMS_no_alt_form[FTICRMS_no_alt_form$dup.mass == 0,]
#looks good, 70 peaks removed, 10797 peaks remain
#this is final set of 10797 peaks for analysis
#looks good

#separate metadata from MS intensities
colnames(FTICRMS_no_dupes)
peak_meta_new <- FTICRMS_no_dupes[, c(1:17, 22:55, 172:185)]
#reorder columns
colnames(peak_meta_new)
peak_meta_new <- peak_meta_new[, c(1:17, 58:63, 18:24, 52:57, 25:37, 64:65, 38:51)]
peak_meta_t <- as.data.frame(t(peak_meta_new))

intensities <- FTICRMS_no_dupes[, c(57:61, 67:90, 96:171)]
rownames(intensities) <- FTICRMS_no_dupes$mz

#normalize intensities
intens_t <- t(intensities)
total_intensities <- rowSums(intens_t)
norm_intens_t_all <- intens_t / total_intensities
write.csv(norm_intens_t_all, file = 'norm_intens_t_all.csv', row.names = TRUE)
#this will be used to make sample MS distance matrix
#and peak based clustering diagram in jmp

#reattach metadata and un-transpose
colnames(peak_meta_t) <- colnames(norm_intens_t_all)
FTICRMS_norm_t <- as.data.frame(rbind(peak_meta_t, norm_intens_t_all))
FTICRMS_norm <- as.data.frame(t(FTICRMS_norm_t))
FTICRMS <- type.convert(FTICRMS_norm, as.is = TRUE)
write.csv(FTICRMS, file = 'FTICRMS.csv', row.names = TRUE)
#this is Table S3




#groupings and their averages
#use transposed version so samples are rows, add metadata
#extract rows where condition is met, average down all peak columns for those extracted rows
#make that a new row or separate object with averages for each peak
MS_w_samp_meta <- cbind(samp_meta, norm_intens_t_all)
avg_tf_MS <- as.data.frame(colMeans(MS_w_samp_meta[MS_w_samp_meta$Sample_type == "Throughfall", -(1:63)]))
avg_precip_all_MS <- as.data.frame(colMeans(MS_w_samp_meta[MS_w_samp_meta$Sample_type == "Precipitation", -(1:63)]))
avg_precip_MS <- as.data.frame(colMeans(MS_w_samp_meta[c("P2", "P3", "P4"), -(1:63)]))
colnames(avg_precip_MS) <- "precip_avg"
colnames(avg_precip_all_MS) <- "precip_all_avg"
colnames(avg_tf_MS) <- "throughfall_avg"
#these three are the most distinctively rain-like, they separate from throughfall samples based on the db-rda analysis
#so just these three will be used for calculating enrichment

#enrichment as subtraction of rainwater signal
avg_precip_105col <- do.call("cbind", replicate(105, avg_precip_MS, simplify = FALSE))
enrich_sample_subt_wo_meta <- FTICRMS[, 66:170] - avg_precip_105col
enrich_sample_subt <- cbind(peak_meta_new, enrich_sample_subt_wo_meta)
write.csv(enrich_sample_subt, file = 'enrich_sample_subt.csv', row.names = TRUE)

#enrichment as division by rainwater signal
enrich_sample_div_wo_meta <- FTICRMS[, 66:170] / avg_precip_105col
#zero over a number gives zero (good),
#a number over 0 gives infinity- I want it to be a number, I'll use the maximum of the non-infinity values in that sample group
#zero over zero gives NaN- I want it to be zero because it wasn't detected in throughfall samples, so wasn't enriched at all
enrich_sample_div_wo_meta[is.na(enrich_sample_div_wo_meta)] <- 0
enrich_sample_div_no_NA <- enrich_sample_div_wo_meta %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.), max(.[is.finite(.)]), .)))
enrich_sample_div <- cbind(peak_meta_new, enrich_sample_div_no_NA)
write.csv(enrich_sample_div, file = 'enrich_sample_div.csv', row.names = TRUE)

avg_precip_w_meta <- cbind(peak_meta_new, avg_precip_MS)


#averages by class
avg_pinopsida_MS <- as.data.frame(colMeans(MS_w_samp_meta[MS_w_samp_meta$Class == "Pinopsida", -(1:63)]))
avg_dicotyledoneae_MS <- as.data.frame(colMeans(MS_w_samp_meta[MS_w_samp_meta$Class == "Dicotyledoneae", -(1:63)]))
avg_classes <- cbind(avg_pinopsida_MS, avg_dicotyledoneae_MS)
colnames(avg_classes) <- c("Pinopsida", "Dicotyledoneae")
MS_classes <- cbind(peak_meta_new, avg_classes)

#enrichment by class (relative to rainwater)
avg_precip_2col <- cbind(avg_precip_MS, avg_precip_MS)
enrich_class_subt <- avg_classes - avg_precip_2col
enrich_class_div <- avg_classes / avg_precip_2col
enrich_class_div[is.na(enrich_class_div)] <- 0
enrich_class_div$Pinopsida[sapply(enrich_class_div$Pinopsida, is.infinite)] <- 
  max(enrich_class_div$Pinopsida[sapply(enrich_class_div$Pinopsida, is.finite)])
enrich_class_div$Dicotyledoneae[sapply(enrich_class_div$Dicotyledoneae, is.infinite)] <- 
  max(enrich_class_div$Dicotyledoneae[sapply(enrich_class_div$Dicotyledoneae, is.finite)])
colnames(enrich_class_subt) <- c("Pinopsida_subt", "Dicotyledoneae_subt")
colnames(enrich_class_div) <- c("Pinopsida_div", "Dicotyledoneae_div")
enrich_class_w_meta <- cbind(peak_meta_new, avg_classes, enrich_class_subt, enrich_class_div)
write.csv(enrich_class_w_meta, file = 'enrich_class_w_meta.csv', row.names = TRUE)


#averages by order
avg_orders_t <- MS_w_samp_meta %>%
  filter(Order != "Precipitation") %>%
  group_by(Order) %>%
  summarise(across(-(1:62), ~mean(.)))
avg_orders <- as.data.frame(t(avg_orders_t[, -1]))
colnames(avg_orders) <- avg_orders_t$Order
MS_orders <- cbind(peak_meta_new, avg_orders)

#enrichment by order
avg_precip_12col <- do.call("cbind", replicate(12, avg_precip_MS, simplify = FALSE))
enrich_order_subt <- avg_orders - avg_precip_12col
enrich_order_div <- avg_orders / avg_precip_12col
enrich_order_div[is.na(enrich_order_div)] <- 0
enrich_order_div <- enrich_order_div %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.), max(.[is.finite(.)]), .)))
enrich_order_w_meta <- cbind(
  peak_meta_new,
  avg_orders %>% rename_with(~ paste0(., "_avg")),
  enrich_order_subt %>% rename_with(~ paste0(., "_subt")),
  enrich_order_div %>% rename_with(~ paste0(., "_div")))
write.csv(enrich_order_w_meta, file = 'enrich_order_w_meta.csv', row.names = TRUE)


#averages by family
avg_families_t <- MS_w_samp_meta %>%
  filter(Family != "Precipitation") %>%
  group_by(Family) %>%
  summarise(across(-(1:62), ~mean(.)))
avg_families <- as.data.frame(t(avg_families_t[, -1]))
colnames(avg_families) <- avg_families_t$Family

#enrichment by family
avg_precip_16col <- do.call("cbind", replicate(16, avg_precip_MS, simplify = FALSE))
enrich_family_subt <- avg_families - avg_precip_16col
enrich_family_div <- avg_families / avg_precip_16col
enrich_family_div[is.na(enrich_family_div)] <- 0
enrich_family_div <- enrich_family_div %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.), max(.[is.finite(.)]), .)))
enrich_family_w_meta <- cbind(
  peak_meta_new,
  avg_families %>% rename_with(~ paste0(., "_avg")),
  enrich_family_subt %>% rename_with(~ paste0(., "_subt")),
  enrich_family_div %>% rename_with(~ paste0(., "_div")))
write.csv(enrich_family_w_meta, file = 'enrich_family_w_meta.csv', row.names = TRUE)


#averages by genus
avg_genera_t <- MS_w_samp_meta %>%
  filter(Genus != "Precipitation") %>%
  group_by(Genus) %>%
  summarise(across(-(1:62), ~mean(.)))
avg_genera <- as.data.frame(t(avg_genera_t[, -1]))
colnames(avg_genera) <- avg_genera_t$Genus

#enrichment by genus
avg_precip_24col <- do.call("cbind", replicate(24, avg_precip_MS, simplify = FALSE))
enrich_genus_subt <- avg_genera - avg_precip_24col
enrich_genus_div <- avg_genera / avg_precip_24col
enrich_genus_div[is.na(enrich_genus_div)] <- 0
enrich_genus_div <- enrich_genus_div %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.), max(.[is.finite(.)]), .)))
enrich_genus_w_meta <- cbind(
  peak_meta_new,
  avg_genera %>% rename_with(~ paste0(., "_avg")),
  enrich_genus_subt %>% rename_with(~ paste0(., "_subt")),
  enrich_genus_div %>% rename_with(~ paste0(., "_div")))
write.csv(enrich_genus_w_meta, file = 'enrich_genus_w_meta.csv', row.names = TRUE)


#averages by species
avg_species_t <- MS_w_samp_meta %>%
  filter(Species != "Precipitation") %>%
  group_by(Species) %>%
  summarise(across(-(1:62), ~mean(.)))
avg_species <- as.data.frame(t(avg_species_t[, -1]))
colnames(avg_species) <- avg_species_t$Species

#enrichment by species
avg_precip_36col <- do.call("cbind", replicate(36, avg_precip_MS, simplify = FALSE))
enrich_species_subt <- avg_species - avg_precip_36col
enrich_species_div <- avg_species / avg_precip_36col
enrich_species_div[is.na(enrich_species_div)] <- 0
enrich_species_div <- enrich_species_div %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.), max(.[is.finite(.)]), .)))
enrich_species_w_meta <- cbind(
  peak_meta_new,
  avg_species %>% rename_with(~ paste0(., "_avg")),
  enrich_species_subt %>% rename_with(~ paste0(., "_subt")),
  enrich_species_div %>% rename_with(~ paste0(., "_div")))
write.csv(enrich_species_w_meta, file = 'enrich_species_w_meta.csv', row.names = TRUE)
write.csv(avg_classes, file = 'avg_classes.csv', row.names = TRUE)
write.csv(avg_orders, file = 'avg_orders.csv', row.names = TRUE)
#these are used to make van Krevelens in Figure 3



#number of formulae in different sample groups
table(sign(avg_precip_all_MS$precip_all_avg))
table(sign(avg_tf_MS$throughfall_avg))

MS_w_samp_meta$formula_counts <- rowSums(MS_w_samp_meta[, -(1:63)] > 0)
plot(x = MS_w_samp_meta$DOC_mgL, y = MS_w_samp_meta$formula_counts)
cor(MS_w_samp_meta$DOC_mgL, MS_w_samp_meta$formula_counts, use = "complete.obs")
cor.test(MS_w_samp_meta$DOC_mgL, MS_w_samp_meta$formula_counts, use = "complete.obs", method = "pearson")
cor.test(MS_w_samp_meta$DOC_mgL, MS_w_samp_meta$formula_counts, use = "complete.obs", method = "kendall")
cor.test(MS_w_samp_meta$DOC_mgL, MS_w_samp_meta$formula_counts, use = "complete.obs", method = "spearman")
cor(MS_w_samp_meta$DOCflux_mg_m2, MS_w_samp_meta$formula_counts, use = "complete.obs")
cor.test(MS_w_samp_meta$DOCflux_mg_m2, MS_w_samp_meta$formula_counts, use = "complete.obs", method = "pearson")
cor.test(MS_w_samp_meta$DOCflux_mg_m2, MS_w_samp_meta$formula_counts, use = "complete.obs", method = "kendall")
cor.test(MS_w_samp_meta$DOCflux_mg_m2, MS_w_samp_meta$formula_counts, use = "complete.obs", method = "spearman")
#linear correlations between DOC and formula counts are not significant
#but kendall and spearman are significant

mean(MS_w_samp_meta$formula_counts[MS_w_samp_meta$Sample_type == "Precipitation"])
std.error(MS_w_samp_meta$formula_counts[MS_w_samp_meta$Sample_type == "Precipitation"])
min(MS_w_samp_meta$formula_counts[MS_w_samp_meta$Sample_type == "Precipitation"])
max(MS_w_samp_meta$formula_counts[MS_w_samp_meta$Sample_type == "Precipitation"])

mean(MS_w_samp_meta$formula_counts[MS_w_samp_meta$Sample_type == "Throughfall"])
std.error(MS_w_samp_meta$formula_counts[MS_w_samp_meta$Sample_type == "Throughfall"])
min(MS_w_samp_meta$formula_counts[MS_w_samp_meta$Sample_type == "Throughfall"])
max(MS_w_samp_meta$formula_counts[MS_w_samp_meta$Sample_type == "Throughfall"])
MS_w_samp_meta$formula_counts[MS_w_samp_meta$formula_counts < 2978.6 & MS_w_samp_meta$Sample_type == "Throughfall"]
MS_w_samp_meta$formula_counts[MS_w_samp_meta$formula_counts < 3968 & MS_w_samp_meta$Sample_type == "Throughfall"]

avg_formulae_by_species <- MS_w_samp_meta %>%
  filter(Sample_type == "Throughfall") %>%
  group_by(Species) %>%
  summarise(
    n = n(),
    mean = mean(formula_counts),
    std_error = ifelse(n > 1, sd(formula_counts) / sqrt(n), NA)
  )
print(avg_formulae_by_species, n = 36)
min(avg_formulae_by_species$mean)
max(avg_formulae_by_species$mean)

grouped_formulae_by_species <- colSums(avg_species > 0)
grouped_formulae_by_species
min(grouped_formulae_by_species)
max(grouped_formulae_by_species)
mean(grouped_formulae_by_species)
std.error(grouped_formulae_by_species)
grouped_formulae_by_species[grouped_formulae_by_species<4979]

#test if condensed aromatics are large molecules
set.seed(0)
t.test(FTICRMS$reference[FTICRMS$Condensed.aromatic == 1],
       FTICRMS$reference[FTICRMS$Condensed.aromatic == 0])
set.seed(0)
t.test(FTICRMS$reference[FTICRMS$Condensed.aromatic == 1],
       FTICRMS$reference[FTICRMS$Condensed.aromatic == 0 & FTICRMS$Aromatic == 1])
#no, they're actually smaller than average and smaller than other aromatics



#calculate % of intensity in each molecular class for each sample
pct.cond.arom <- FTICRMS %>%
  filter(Condensed.aromatic == 1) %>%
  summarise(across(66:170, sum))
pct.cond.arom.O.rich <- FTICRMS %>%
  filter(Condensed.aromatic.O_rich == 1) %>%
  summarise(across(66:170, sum))
pct.cond.arom.O.poor <- FTICRMS %>%
  filter(Condensed.aromatic.O_poor == 1) %>%
  summarise(across(66:170, sum))
pct.uncond.arom <- FTICRMS %>%
  filter(Uncondensed.aromatic == 1) %>%
  summarise(across(66:170, sum))
pct.uncond.arom.O.rich <- FTICRMS %>%
  filter(Uncondensed.aromatic.O_rich == 1) %>%
  summarise(across(66:170, sum))
pct.uncond.arom.O.poor <- FTICRMS %>%
  filter(Uncondensed.aromatic.O_poor == 1) %>%
  summarise(across(66:170, sum))
pct.arom <- FTICRMS %>%
  filter(Aromatic == 1) %>%
  summarise(across(66:170, sum))
pct.arom.O.rich <- FTICRMS %>%
  filter(Aromatic.O_rich == 1) %>%
  summarise(across(66:170, sum))
pct.arom.O.poor <- FTICRMS %>%
  filter(Aromatic.O_poor == 1) %>%
  summarise(across(66:170, sum))
pct.high.unsat <- FTICRMS %>%
  filter(Highly.unsaturated == 1) %>%
  summarise(across(66:170, sum))
pct.high.unsat.O.rich <- FTICRMS %>%
  filter(Highly.unsaturated.O_rich == 1) %>%
  summarise(across(66:170, sum))
pct.high.unsat.O.poor <- FTICRMS %>%
  filter(Highly.unsaturated.O_poor == 1) %>%
  summarise(across(66:170, sum))
pct.unsat <- FTICRMS %>%
  filter(Unsaturated == 1) %>%
  summarise(across(66:170, sum))
pct.unsat.O.rich <- FTICRMS %>%
  filter(Unsaturated.O_rich == 1) %>%
  summarise(across(66:170, sum))
pct.unsat.O.poor <- FTICRMS %>%
  filter(Unsaturated.O_poor == 1) %>%
  summarise(across(66:170, sum))
pct.unsat.with.N <- FTICRMS %>%
  filter(Unsaturated.with.N == 1) %>%
  summarise(across(66:170, sum))
pct.saturated <- FTICRMS %>%
  filter(Saturated == 1) %>%
  summarise(across(66:170, sum))
pct.saturated.O.rich <- FTICRMS %>%
  filter(Saturated.O_rich == 1) %>%
  summarise(across(66:170, sum))
pct.saturated.O.poor <- FTICRMS %>%
  filter(Saturated.O_poor == 1) %>%
  summarise(across(66:170, sum))
pct.sugar.like <- FTICRMS %>%
  filter(sugar.like == 1) %>%
  summarise(across(66:170, sum))
pct.peptide.like <- FTICRMS %>%
  filter(peptide.like == 1) %>%
  summarise(across(66:170, sum))
pct.O.rich <- FTICRMS %>%
  filter(O.C > 0.5) %>%
  summarise(across(66:170, sum))
pct.O.poor <- FTICRMS %>%
  filter(O.C <= 0.5) %>%
  summarise(across(66:170, sum))

pct.CHO <- FTICRMS %>%
  filter(CHO == 1) %>%
  summarise(across(66:170, sum))
pct.CHON <- FTICRMS %>%
  filter(CHON == 1) %>%
  summarise(across(66:170, sum))
pct.CHOS <- FTICRMS %>%
  filter(CHOS == 1) %>%
  summarise(across(66:170, sum))
pct.CHOP <- FTICRMS %>%
  filter(CHOP == 1) %>%
  summarise(across(66:170, sum))
pct.CHONS <- FTICRMS %>%
  filter(CHONS == 1) %>%
  summarise(across(66:170, sum))

mean.MW <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * reference)))
mean.H.C <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * H.C)))
mean.O.C <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * O.C)))
mean.AI <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * AI)))
mean.AI.mod <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * AI.mod)))
mean.DBE <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * DBE)))

mean.C <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * C)))
mean.H <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * H)))
mean.O <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * O)))
mean.N <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * N)))
mean.S <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * S)))
mean.P <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * P)))
mean.non.CHO.atoms <- FTICRMS %>%
  summarise(across(66:170, ~ sum(. * non_CHO_atoms)))
mean.N.C <- mean.N / mean.C
mean.S.C <- mean.S / mean.C
mean.P.C <- mean.P / mean.C

#combine into one data frame
FTICRMS_means <- rbind(mean.MW, mean.H.C, mean.O.C, mean.AI, mean.AI.mod, mean.DBE,
                          mean.C, mean.H, mean.O, mean.N, mean.S, mean.P,
                          mean.non.CHO.atoms, mean.N.C, mean.S.C, mean.P.C,
                          pct.CHO, pct.CHON, pct.CHOS, pct.CHOP, pct.CHONS,
                          pct.cond.arom, pct.cond.arom.O.rich, pct.cond.arom.O.poor,
                          pct.uncond.arom, pct.uncond.arom.O.rich, pct.uncond.arom.O.poor,
                          pct.arom, pct.arom.O.rich, pct.arom.O.poor,
                          pct.high.unsat, pct.high.unsat.O.rich, pct.high.unsat.O.poor,
                          pct.unsat, pct.unsat.O.rich, pct.unsat.O.poor, pct.unsat.with.N,
                          pct.saturated, pct.saturated.O.rich, pct.saturated.O.poor,
                          pct.sugar.like, pct.peptide.like, pct.O.rich, pct.O.poor)
chemical_parameters <- c("mean.MW", "mean.H.C", "mean.O.C", "mean.AI", "mean.AI.mod", "mean.DBE",
                         "mean.C", "mean.H", "mean.O", "mean.N", "mean.S", "mean.P",
                         "mean.non.CHO.atoms", "mean.N.C", "mean.S.C", "mean.P.C",
                         "pct.CHO", "pct.CHON", "pct.CHOS", "pct.CHOP", "pct.CHONS",
                         "pct.cond.arom", "pct.cond.arom.O.rich", "pct.cond.arom.O.poor",
                         "pct.uncond.arom", "pct.uncond.arom.O.rich", "pct.uncond.arom.O.poor",
                         "pct.arom", "pct.arom.O.rich", "pct.arom.O.poor",
                         "pct.high.unsat", "pct.high.unsat.O.rich", "pct.high.unsat.O.poor",
                         "pct.unsat", "pct.unsat.O.rich", "pct.unsat.O.poor", "pct.unsat.with.N",
                         "pct.saturated", "pct.saturated.O.rich", "pct.saturated.O.poor",
                         "pct.sugar.like", "pct.peptide.like", "pct.O.rich", "pct.O.poor")
rownames(FTICRMS_means) <- chemical_parameters
#looks good

#also want standard deviations for SI table- for ones where that's applicable
#which are ones that calculate mean value instead of % of intensity in a category
weighted_se <- function(x, w) {
  weighted_mean <- sum(w * x)
  weighted_var <- sum(w * (x - weighted_mean)^2)
  sqrt(weighted_var) / sqrt(length(x))
}

se.MW <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(reference, .)))
se.H.C <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(H.C, .)))
se.O.C <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(O.C, .)))
se.AI <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(AI, .)))
se.AI.mod <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(AI.mod, .)))
se.DBE <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(DBE, .)))

se.C <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(C, .)))
se.H <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(H, .)))
se.O <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(O, .)))
se.N <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(N, .)))
se.S <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(S, .)))
se.P <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(P, .)))
se.non.CHO.atoms <- FTICRMS %>%
  summarise(across(66:170, ~ weighted_se(non_CHO_atoms, .)))

FTICRMS_means_and_SEs <- rbind(mean.MW, se.MW, mean.H.C, se.H.C, mean.O.C, se.O.C,
                              mean.AI, se.AI, mean.AI.mod, se.AI.mod, mean.DBE, se.DBE,
                              mean.C, se.C, mean.H, se.H, mean.O, se.O,
                              mean.N, se.N, mean.S, se.S, mean.P, se.P,
                              mean.non.CHO.atoms, se.non.CHO.atoms, mean.N.C, mean.S.C, mean.P.C,
                              pct.CHO, pct.CHON, pct.CHOS, pct.CHOP, pct.CHONS,
                              pct.cond.arom, pct.cond.arom.O.rich, pct.cond.arom.O.poor,
                              pct.uncond.arom, pct.uncond.arom.O.rich, pct.uncond.arom.O.poor,
                              pct.arom, pct.arom.O.rich, pct.arom.O.poor,
                              pct.high.unsat, pct.high.unsat.O.rich, pct.high.unsat.O.poor,
                              pct.unsat, pct.unsat.O.rich, pct.unsat.O.poor, pct.unsat.with.N,
                              pct.saturated, pct.saturated.O.rich, pct.saturated.O.poor,
                              pct.sugar.like, pct.peptide.like, pct.O.rich, pct.O.poor)
rownames(FTICRMS_means_and_SEs) <- c("mean.MW", "se.MW", "mean.H.C", "se.H.C", "mean.O.C", "se.O.C",
                                     "mean.AI", "se.AI", "mean.AI.mod", "se.AI.mod", "mean.DBE", "se.DBE",
                                     "mean.C", "se.C", "mean.H", "se.H", "mean.O", "se.O",
                                     "mean.N", "se.N", "mean.S", "se.S", "mean.P", "se.P",
                                     "mean.non.CHO.atoms", "se.non.CHO.atoms", "mean.N.C", "mean.S.C", "mean.P.C",
                                     "pct.CHO", "pct.CHON", "pct.CHOS", "pct.CHOP", "pct.CHONS",
                                     "pct.cond.arom", "pct.cond.arom.O.rich", "pct.cond.arom.O.poor",
                                     "pct.uncond.arom", "pct.uncond.arom.O.rich", "pct.uncond.arom.O.poor",
                                     "pct.arom", "pct.arom.O.rich", "pct.arom.O.poor",
                                     "pct.high.unsat", "pct.high.unsat.O.rich", "pct.high.unsat.O.poor",
                                     "pct.unsat", "pct.unsat.O.rich", "pct.unsat.O.poor", "pct.unsat.with.N",
                                     "pct.saturated", "pct.saturated.O.rich", "pct.saturated.O.poor",
                                     "pct.sugar.like", "pct.peptide.like", "pct.O.rich", "pct.O.poor")
bulk_all <- as.data.frame(t(FTICRMS_means_and_SEs))
bulk_all <- bulk_all[c(54:58, 1:53, 59:105),]
write.csv(bulk_all, file = 'bulk_parameters_all_samples.csv', row.names = TRUE)


#report parameters by sample before grouping
sample_stats_tf <- function(row) {
  x <- unlist(FTICRMS_means[row, -(54:58)])
  data.frame(
    min = min(x),
    max = max(x),
    mean = mean(x),
    std_error = std.error(x)
  )
}
sample_stats_tf_results <- do.call(rbind, lapply(1:44, sample_stats_tf))
rownames(sample_stats_tf_results) <- chemical_parameters
sample_stats_tf_results

sample_stats_precip <- function(row) {
  x <- unlist(FTICRMS_means[row, 54:58])
  data.frame(
    min = min(x),
    max = max(x),
    mean = mean(x),
    std_error = std.error(x)
  )
}
sample_stats_precip_results <- do.call(rbind, lapply(1:44, sample_stats_precip))
rownames(sample_stats_precip_results) <- chemical_parameters
sample_stats_precip_results

#tests for significant differences between throughfall and precipitation
sample_stats_t_test <- apply(FTICRMS_means, 1, function(row) {
  tf <- as.numeric(row[-(54:58)])
  precip <- as.numeric(row[54:58])
  test <- t.test(tf, precip)
  data.frame(
    p_value = test$p.value,
    t_stat = test$statistic,
    mean_tf = test$estimate[1],
    mean_precip = test$estimate[2]
  )
})
sample_stats_t_test_results <- do.call(rbind, sample_stats_t_test)


#add metadata for groupings
FTICRMS_means_t <- as.data.frame(t(FTICRMS_means))
means_w_samp_meta <- cbind(samp_meta, FTICRMS_means_t)
write.csv(means_w_samp_meta, file = 'means_w_samp_meta.csv', row.names = FALSE)
#this is used to make clustering heatmap in jmp (Figure S3)

#groupings- averages and st errors by species, genus, etc.
bulk_avgs_class <- means_w_samp_meta %>%
  group_by(Class) %>%
  summarise(across(all_of(chemical_parameters),
                   list(mean = ~mean(.),
                        se = ~sd(.) / sqrt(sum(!is.na(.))))))
bulk_avgs_class <- as.data.frame(bulk_avgs_class)
bulk_avgs_order <- means_w_samp_meta %>%
  group_by(Order) %>%
  summarise(across(all_of(chemical_parameters),
                   list(mean = ~mean(.),
                        se = ~sd(.) / sqrt(sum(!is.na(.))))))
bulk_avgs_order <- as.data.frame(bulk_avgs_order)
bulk_avgs_family <- means_w_samp_meta %>%
  group_by(Family) %>%
  summarise(across(all_of(chemical_parameters),
                   list(mean = ~mean(.),
                        se = ~sd(.) / sqrt(sum(!is.na(.))))))
bulk_avgs_family <- as.data.frame(bulk_avgs_family)
bulk_avgs_genus <- means_w_samp_meta %>%
  group_by(Genus) %>%
  summarise(across(all_of(chemical_parameters),
                   list(mean = ~mean(.),
                        se = ~sd(.) / sqrt(sum(!is.na(.))))))
bulk_avgs_genus <- as.data.frame(bulk_avgs_genus)
bulk_avgs_species <- means_w_samp_meta %>%
  group_by(Species) %>%
  summarise(across(all_of(chemical_parameters),
                   list(mean = ~mean(.),
                        se = ~sd(.) / sqrt(sum(!is.na(.))))))
bulk_avgs_species <- as.data.frame(bulk_avgs_species)
bulk_avgs_all_tf <- means_w_samp_meta %>%
  filter(Sample_type == "Throughfall") %>%
  summarise(across(all_of(chemical_parameters),
                   list(mean = ~mean(.),
                        se = ~sd(.) / sqrt(sum(!is.na(.))))))
bulk_avgs_precip <- means_w_samp_meta %>%
  filter(Sample_type == "Precipitation") %>%
  summarise(across(all_of(chemical_parameters),
                   list(mean = ~mean(.),
                        se = ~sd(.) / sqrt(sum(!is.na(.))))))
write.csv(bulk_avgs_class, file = 'bulk_avgs_class.csv')
write.csv(bulk_avgs_order, file = 'bulk_avgs_order.csv')
write.csv(bulk_avgs_family, file = 'bulk_avgs_family.csv')
write.csv(bulk_avgs_genus, file = 'bulk_avgs_genus.csv')
write.csv(bulk_avgs_species, file = 'bulk_avgs_species.csv')
write.csv(bulk_avgs_all_tf, file = 'bulk_avgs_all_tf.csv')
write.csv(bulk_avgs_precip, file = 'bulk_avgs_precip.csv')
#species, all_tf, and precip are the data in Table 1
#and all of these go into Table S1

#want enrichment file of everything to export too
rownames(bulk_avgs_class) <- bulk_avgs_class[, 1]
rownames(bulk_avgs_order) <- bulk_avgs_order[, 1]
rownames(bulk_avgs_family) <- bulk_avgs_family[, 1]
rownames(bulk_avgs_genus) <- bulk_avgs_genus[, 1]
rownames(bulk_avgs_species) <- bulk_avgs_species[, 1]
bulk_avgs_combined <- rbind(bulk_avgs_class[, -1], bulk_avgs_order[, -1], bulk_avgs_family[, -1], bulk_avgs_genus[, -1], bulk_avgs_species[, -1])
bulk_avgs_all_tf_95row <- do.call("rbind", replicate(95, bulk_avgs_all_tf, simplify = FALSE))
bulk_avgs_enrich <- bulk_avgs_combined - bulk_avgs_all_tf_95row
bulk_avgs_pct_enrich <- bulk_avgs_enrich / bulk_avgs_all_tf_95row
write.csv(bulk_avgs_pct_enrich, file = 'bulk_avgs_pct_enrich.csv', row.names = TRUE)
#also want tf average, precip average, and species averages for SI
SI_table_species_bulk <- rbind(bulk_avgs_all_tf, bulk_avgs_precip, bulk_avgs_species[-7, -1])
write.csv(SI_table_species_bulk, file = 'SI_table_species_bulk.csv')

#test if conifers vs broadleaves have different FTICRMS properties
set.seed(0)
t.test(means_w_samp_meta$mean.N.C[means_w_samp_meta$Class == "Pinopsida"],
       means_w_samp_meta$mean.N.C[means_w_samp_meta$Class == "Dicotyledoneae"])
mean(means_w_samp_meta$mean.N.C[means_w_samp_meta$Class == "Pinopsida"])
std.error(means_w_samp_meta$mean.N.C[means_w_samp_meta$Class == "Pinopsida"])
mean(means_w_samp_meta$mean.N.C[means_w_samp_meta$Class == "Dicotyledoneae"])
std.error(means_w_samp_meta$mean.N.C[means_w_samp_meta$Class == "Dicotyledoneae"])

set.seed(0)
t.test(means_w_samp_meta$mean.S.C[means_w_samp_meta$Class == "Pinopsida"],
       means_w_samp_meta$mean.S.C[means_w_samp_meta$Class == "Dicotyledoneae"])
mean(means_w_samp_meta$mean.S.C[means_w_samp_meta$Class == "Pinopsida"])
std.error(means_w_samp_meta$mean.S.C[means_w_samp_meta$Class == "Pinopsida"])
mean(means_w_samp_meta$mean.S.C[means_w_samp_meta$Class == "Dicotyledoneae"])
std.error(means_w_samp_meta$mean.S.C[means_w_samp_meta$Class == "Dicotyledoneae"])

set.seed(0)
t.test(means_w_samp_meta$mean.P.C[means_w_samp_meta$Class == "Pinopsida"],
       means_w_samp_meta$mean.P.C[means_w_samp_meta$Class == "Dicotyledoneae"])
mean(means_w_samp_meta$mean.P.C[means_w_samp_meta$Class == "Pinopsida"])
std.error(means_w_samp_meta$mean.P.C[means_w_samp_meta$Class == "Pinopsida"])
mean(means_w_samp_meta$mean.P.C[means_w_samp_meta$Class == "Dicotyledoneae"])
std.error(means_w_samp_meta$mean.P.C[means_w_samp_meta$Class == "Dicotyledoneae"])

set.seed(0)
t.test(means_w_samp_meta$mean.MW[means_w_samp_meta$Class == "Pinopsida"],
       means_w_samp_meta$mean.MW[means_w_samp_meta$Class == "Dicotyledoneae"])
set.seed(0)
t.test(means_w_samp_meta$mean.H.C[means_w_samp_meta$Class == "Pinopsida"],
       means_w_samp_meta$mean.H.C[means_w_samp_meta$Class == "Dicotyledoneae"])
set.seed(0)
t.test(means_w_samp_meta$mean.O.C[means_w_samp_meta$Class == "Pinopsida"],
       means_w_samp_meta$mean.O.C[means_w_samp_meta$Class == "Dicotyledoneae"])
set.seed(0)
t.test(means_w_samp_meta$mean.AI.mod[means_w_samp_meta$Class == "Pinopsida"],
       means_w_samp_meta$mean.AI.mod[means_w_samp_meta$Class == "Dicotyledoneae"])

#test if the only deciduous conifer species is different from other conifers
set.seed(0)
t.test(means_w_samp_meta$mean.N.C[means_w_samp_meta$Species == "glyptostroboides"],
       means_w_samp_meta$mean.N.C[means_w_samp_meta$Class == "Pinopsida" & means_w_samp_meta$Species != "glyptostroboides"])
set.seed(0)
t.test(means_w_samp_meta$mean.N[means_w_samp_meta$Species == "glyptostroboides"],
       means_w_samp_meta$mean.N[means_w_samp_meta$Class == "Pinopsida" & means_w_samp_meta$Species != "glyptostroboides"])
set.seed(0)
t.test(means_w_samp_meta$mean.P.C[means_w_samp_meta$Species == "glyptostroboides"],
       means_w_samp_meta$mean.P.C[means_w_samp_meta$Class == "Pinopsida" & means_w_samp_meta$Species != "glyptostroboides"])
set.seed(0)
t.test(means_w_samp_meta$mean.P[means_w_samp_meta$Species == "glyptostroboides"],
       means_w_samp_meta$mean.P[means_w_samp_meta$Class == "Pinopsida" & means_w_samp_meta$Species != "glyptostroboides"])
set.seed(0)
t.test(means_w_samp_meta$mean.S.C[means_w_samp_meta$Species == "glyptostroboides"],
       means_w_samp_meta$mean.S.C[means_w_samp_meta$Class == "Pinopsida" & means_w_samp_meta$Species != "glyptostroboides"])
set.seed(0)
t.test(means_w_samp_meta$mean.S[means_w_samp_meta$Species == "glyptostroboides"],
       means_w_samp_meta$mean.S[means_w_samp_meta$Class == "Pinopsida" & means_w_samp_meta$Species != "glyptostroboides"])
set.seed(0)
t.test(tf_meta$DOC.TDN_ratio[tf_meta$Species == "glyptostroboides"],
       tf_meta$DOC.TDN_ratio[tf_meta$Class == "Pinopsida" & tf_meta$Species != "glyptostroboides"])
set.seed(0)
t.test(tf_meta$TDN_mgL[tf_meta$Species == "glyptostroboides"],
       tf_meta$TDN_mgL[tf_meta$Class == "Pinopsida" & tf_meta$Species != "glyptostroboides"])
set.seed(0)
t.test(tf_meta$DOC_mgL[tf_meta$Species == "glyptostroboides"],
       tf_meta$DOC_mgL[tf_meta$Class == "Pinopsida" & tf_meta$Species != "glyptostroboides"])

#or if only the evergreen broadleaf species is different from other broadleaves
set.seed(0)
t.test(means_w_samp_meta$mean.N.C[means_w_samp_meta$Species == "opaca"],
       means_w_samp_meta$mean.N.C[means_w_samp_meta$Class == "Dicotyledoneae" & means_w_samp_meta$Species != "opaca"])
set.seed(0)
t.test(means_w_samp_meta$mean.N[means_w_samp_meta$Species == "opaca"],
       means_w_samp_meta$mean.N[means_w_samp_meta$Class == "Dicotyledoneae" & means_w_samp_meta$Species != "opaca"])
set.seed(0)
t.test(means_w_samp_meta$mean.P.C[means_w_samp_meta$Species == "opaca"],
       means_w_samp_meta$mean.P.C[means_w_samp_meta$Class == "Dicotyledoneae" & means_w_samp_meta$Species != "opaca"])
set.seed(0)
t.test(means_w_samp_meta$mean.P[means_w_samp_meta$Species == "opaca"],
       means_w_samp_meta$mean.P[means_w_samp_meta$Class == "Dicotyledoneae" & means_w_samp_meta$Species != "opaca"])
set.seed(0)
t.test(means_w_samp_meta$mean.S.C[means_w_samp_meta$Species == "opaca"],
       means_w_samp_meta$mean.S.C[means_w_samp_meta$Class == "Dicotyledoneae" & means_w_samp_meta$Species != "opaca"])
set.seed(0)
t.test(means_w_samp_meta$mean.S[means_w_samp_meta$Species == "opaca"],
       means_w_samp_meta$mean.S[means_w_samp_meta$Class == "Dicotyledoneae" & means_w_samp_meta$Species != "opaca"])
set.seed(0)
t.test(tf_meta$DOC.TDN_ratio[tf_meta$Species == "opaca"],
       tf_meta$DOC.TDN_ratio[tf_meta$Class == "Dicotyledoneae" & tf_meta$Species != "opaca"])
set.seed(0)
t.test(tf_meta$TDN_mgL[tf_meta$Species == "opaca"],
       tf_meta$TDN_mgL[tf_meta$Class == "Dicotyledoneae" & tf_meta$Species != "opaca"])
set.seed(0)
t.test(tf_meta$DOC_mgL[tf_meta$Species == "opaca"],
       tf_meta$DOC_mgL[tf_meta$Class == "Dicotyledoneae" & tf_meta$Species != "opaca"])
#answer is no to both, they look like their class
#at least in terms of N, S, and P content


#test for correlations between environmental variables and bulk chemical parameters
#this can be used to explain the dbrda results
#don't want rainwater samples included, and can't be run on non-numeric variables
corr_input_w_precip <- means_w_samp_meta[sapply(means_w_samp_meta, is.numeric)]
corr_input <- corr_input_w_precip[-(54:58), ]
corr <- as.data.frame(cor(corr_input, method = "pearson", use = "everything"))
#doing DOC separately since it has NA values
corr_DOC <- as.data.frame(cor(corr_input, method = "pearson", use = "complete.obs"))
write.csv(corr, file = "corr.csv", row.names = TRUE)
write.csv(corr_DOC, file = "corr_DOC.csv", row.names = TRUE)

#calculate significance for chemical parameters with r>0.25
#focusing on the environmental variables that are significant in the dbrda
#DOC correlations
rownames(corr_DOC[corr_DOC$DOC_mgL > 0.25 | corr_DOC$DOC_mgL < -0.25,])
cor.test(corr_input$DOC_mgL, corr_input$TDN_mgL, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$mean.MW, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$mean.O.C, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$mean.DBE, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$mean.N, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$mean.MW, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$mean.non.CHO.atoms, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.CHO, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.CHON, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.high.unsat, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.high.unsat.O.rich, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.high.unsat.O.poor, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.unsat.O.rich, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.unsat.with.N, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.saturated, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.saturated.O.rich, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.saturated.O.poor, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.sugar.like, use = "complete.obs")
cor.test(corr_input$DOC_mgL, corr_input$pct.peptide.like, use = "complete.obs")
#all are significant below p<0.05, many below <0.0001

#TDN correlations
rownames(corr[corr$TDN_mgL > 0.25 | corr$TDN_mgL < -0.25,])
cor.test(corr_input$TDN_mgL, corr_input$mean.H.C)
cor.test(corr_input$TDN_mgL, corr_input$mean.AI)
cor.test(corr_input$TDN_mgL, corr_input$mean.AI.mod)
cor.test(corr_input$TDN_mgL, corr_input$mean.DBE)
cor.test(corr_input$TDN_mgL, corr_input$mean.N)
cor.test(corr_input$TDN_mgL, corr_input$mean.S)
cor.test(corr_input$TDN_mgL, corr_input$mean.P)
cor.test(corr_input$TDN_mgL, corr_input$mean.non.CHO.atoms)
cor.test(corr_input$TDN_mgL, corr_input$pct.CHO)
cor.test(corr_input$TDN_mgL, corr_input$pct.CHOS)
cor.test(corr_input$TDN_mgL, corr_input$pct.CHOP)
cor.test(corr_input$TDN_mgL, corr_input$pct.CHONS)
cor.test(corr_input$TDN_mgL, corr_input$pct.cond.arom)
cor.test(corr_input$TDN_mgL, corr_input$pct.cond.arom.O.poor)
cor.test(corr_input$TDN_mgL, corr_input$pct.uncond.arom)
cor.test(corr_input$TDN_mgL, corr_input$pct.uncond.arom.O.poor)
cor.test(corr_input$TDN_mgL, corr_input$pct.arom)
cor.test(corr_input$TDN_mgL, corr_input$pct.arom.O.poor)
cor.test(corr_input$TDN_mgL, corr_input$pct.unsat.with.N)
cor.test(corr_input$TDN_mgL, corr_input$pct.saturated.O.rich)
cor.test(corr_input$TDN_mgL, corr_input$pct.peptide.like)
#all are significant below p<0.01, many below p<0.0001
cor.test(corr_input$TDN_mgL, corr_input$pct.uncond.arom.O.rich)

#throughfall volume correlations
rownames(corr[corr$Tf_depth_mm > 0.25 | corr$Tf_depth_mm < -0.25,])
#only correlates with itself, DOC, and TDN- no FTICRMS chemical parameters

#tree size (dbh) correlations
rownames(corr[corr$DBH_cm > 0.25 | corr$DBH_cm < -0.25,])
cor.test(corr_input$DBH_cm, corr_input$mean.AI.mod)
cor.test(corr_input$DBH_cm, corr_input$mean.DBE)
cor.test(corr_input$DBH_cm, corr_input$mean.N)
cor.test(corr_input$DBH_cm, corr_input$mean.S)
cor.test(corr_input$DBH_cm, corr_input$mean.non.CHO.atoms)
cor.test(corr_input$DBH_cm, corr_input$pct.CHO)
cor.test(corr_input$DBH_cm, corr_input$pct.CHON)
cor.test(corr_input$DBH_cm, corr_input$pct.CHOS)
cor.test(corr_input$DBH_cm, corr_input$pct.unsat.O.rich)
cor.test(corr_input$DBH_cm, corr_input$pct.unsat.with.N)
cor.test(corr_input$DBH_cm, corr_input$pct.saturated.O.rich)
cor.test(corr_input$DBH_cm, corr_input$pct.peptide.like)
#all are significant below p<0.05, many below p<0.0001
#this is overall what distinguishes tf from precip, testing for a couple others that agree with this
cor.test(corr_input$DBH_cm, corr_input$mean.P)
cor.test(corr_input$DBH_cm, corr_input$mean.MW)
#molecular weight significant, P content not

#building distance correlations
rownames(corr[corr$BLDG_dist > 0.25 | corr$BLDG_dist < -0.25,])
cor.test(corr_input$BLDG_dist, corr_input$pct.CHO)
cor.test(corr_input$BLDG_dist, corr_input$pct.unsat.O.rich)
#both are significant below p<0.01

#testing for road/train distance correlations
#to see if they affect chemical properties
rownames(corr[corr$ROAD_dist > 0.25 | corr$ROAD_dist < -0.25,])
rownames(corr_DOC[corr_DOC$ROAD_dist > 0.25 | corr_DOC$ROAD_dist < -0.25,])
rownames(corr[corr$TRAIN_dist > 0.25 | corr$TRAIN_dist < -0.25,])
rownames(corr_DOC[corr_DOC$TRAIN_dist > 0.25 | corr_DOC$TRAIN_dist < -0.25,])
#none have strong explanatory power
cor.test(corr_input$ROAD_dist, corr_input$pct.cond.arom)
cor.test(corr_input$ROAD_dist, corr_input$pct.cond.arom.O.rich)
cor.test(corr_input$ROAD_dist, corr_input$pct.cond.arom.O.poor)
cor.test(corr_input$TRAIN_dist, corr_input$pct.cond.arom)
cor.test(corr_input$TRAIN_dist, corr_input$pct.cond.arom.O.rich)
cor.test(corr_input$TRAIN_dist, corr_input$pct.cond.arom.O.poor)
#distance from roads/trains doesn't significantly correlate with condensed aromatics (black carbon)
sample_stats_tf_results["pct.cond.arom",]
sample_stats_precip_results["pct.cond.arom",]
min(tf_meta$ROAD_dist)
max(tf_meta$ROAD_dist)
min(tf_meta$TRAIN_dist)
max(tf_meta$TRAIN_dist)

#look for impact from other potential factors on black carbon content
rownames(corr[corr$pct.cond.arom > 0.25 | corr$pct.cond.arom < -0.25,])
#only DOC, TDN, and other chemical parameters have r>0.25
#testing building distance, tree size, species, class
cor.test(corr_input$BLDG_dist, corr_input$pct.cond.arom)
cor.test(corr_input$BLDG_dist, corr_input$pct.cond.arom.O.rich)
cor.test(corr_input$BLDG_dist, corr_input$pct.cond.arom.O.poor)
cor.test(corr_input$DBH_cm, corr_input$pct.cond.arom)
cor.test(corr_input$DBH_cm, corr_input$pct.cond.arom.O.rich)
cor.test(corr_input$DBH_cm, corr_input$pct.cond.arom.O.poor)
cor.test(corr_input$Tf_depth_mm, corr_input$pct.cond.arom)
cor.test(corr_input$Tf_depth_mm, corr_input$pct.cond.arom.O.rich)
cor.test(corr_input$Tf_depth_mm, corr_input$pct.cond.arom.O.poor)
cor.test(corr_input$DOC_mgL, corr_input$pct.cond.arom)
cor.test(corr_input$DOC_mgL, corr_input$pct.cond.arom.O.rich)
cor.test(corr_input$DOC_mgL, corr_input$pct.cond.arom.O.poor)
#none are significant
cor.test(corr_input$TDN_mgL, corr_input$pct.cond.arom)
cor.test(corr_input$TDN_mgL, corr_input$pct.cond.arom.O.rich)
cor.test(corr_input$TDN_mgL, corr_input$pct.cond.arom.O.poor)
#significant negative correlation, following the way TDN is
#negatively correlated with all aromatics
cor.test(corr_input$pct.uncond.arom, corr_input$pct.cond.arom)
cor.test(corr_input$pct.uncond.arom.O.rich, corr_input$pct.cond.arom.O.rich)
cor.test(corr_input$pct.uncond.arom.O.poor, corr_input$pct.cond.arom.O.poor)
#strong positive correlation- condensed aromatics seem to just follow other aromatics
#p<10e-15, r=0.8927868, r2=0.797
set.seed(0)
t.test(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Pinopsida"],
       means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Dicotyledoneae"])
set.seed(0)
t.test(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Leaf_type == "needles"],
       means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Leaf_type == "scales"])
set.seed(0)
t.test(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Leaf_type == "needles"],
       means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Leaf_type == "leaves"])
set.seed(0)
t.test(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Leaf_type == "leaves"],
       means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Leaf_type == "scales"])
mean(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Pinopsida"])
std.error(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Pinopsida"])
min(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Pinopsida"])
max(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Pinopsida"])
mean(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Dicotyledoneae"])
std.error(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Dicotyledoneae"])
min(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Dicotyledoneae"])
max(means_w_samp_meta$pct.cond.arom[means_w_samp_meta$Class == "Dicotyledoneae"])

print(bulk_avgs_species[order(bulk_avgs_species$pct.cond.arom_mean),
                        c("Species", "pct.cond.arom_mean", "pct.cond.arom_se")], n=37)




norm_intens_t <- norm_intens_t_all
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
#105 rows (sample IDs) by 10797 columns (mzs)- good



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
#to compare to sample FTICRMS distance matrix
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
#first make it with parents
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
#get rid of tiny exponents
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
                      "Metasequoia glyptostroboides", "Pinus resinosa", "Pinus strobus", "Abies concolor", "Platanus x. acerifolia", "Magnolia x. soulangeana")
cluster_nj_phylo2 <- nj(phylo_dist_jmp)
cluster_nj_phylo2$tip.label <- cladogram_labels
plot(cluster_nj_phylo2, type = "cladogram", cex = 0.7)
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

#alternative: making it manually to be more taxonomy focused
#and be able to control structure better
#this is final Figure 1
taxonomy_initial <- as.data.frame(tf_meta[c("Phylum", "Class", "Order", "Family", "Genus", "Latin")])
taxonomy_initial$Latin[48:50] <- "Platanus x. acerifolia"
taxonomy_initial$Latin[82:84] <- "Magnolia x. soulangeana"
taxonomy_initial$Kingdom <- "Plantae"
taxonomy_initial <- taxonomy_initial %>% mutate(Phylum_Class = paste(Phylum, Class, sep = "/"))
taxonomy <- taxonomy_initial[, c("Kingdom", "Phylum_Class", "Order", "Family", "Genus", "Latin")]
colnames(taxonomy) <- c("Kingdom", "Phylum_Class", "Order", "Family", "Genus", "Species")
taxonomy$pathString <- paste(taxonomy$Kingdom, taxonomy$Phylum_Class, taxonomy$Order,
                             taxonomy$Family, taxonomy$Genus, taxonomy$Species, sep = "|")
taxonomy_tree <- as.Node(taxonomy, pathDelimiter = "|")
plot(taxonomy_tree)
taxonomy_phylo <- as.phylo(taxonomy_tree)
plot(taxonomy_phylo, type = "cladogram", show.node.label = TRUE,
     cex = 0.7, label.offset = 0.1, no.margin = TRUE)
plot(taxonomy_phylo, type = "phylogram", show.node.label = TRUE, cex = 0.7)
#adjusting placement of names to be on top of lines instead of nodes
tree_data <- ggtree(taxonomy_phylo, ladderize = TRUE)$data
tree_data <- tree_data %>%
  mutate(parent_x = x[match(parent, node)],
         mid_x = (x + parent_x) / 2)
root_row <- tree_data %>% filter(node == parent)
root_x <- root_row$x
root_y <- root_row$y
root_label_data <- data.frame(x = -6, xend = 0, y = root_y, yend = root_y)
nonroot_tree_data <- tree_data %>% filter(node != parent)
nonroot_tree_data <- nonroot_tree_data %>%
  mutate(
    display_label = ifelse(isTip, gsub("_", " ", label), label),
    label_face = ifelse(isTip, "italic", "plain")
  )

depth <- data.frame(1:5, c(10, 30, 50, 70, 90))
colnames(depth) <- c("depth", "header_x")
header_labels <- data.frame(
  depth  = 0:5,
  rank_name = c("Kingdom", "Phylum/Class", "Order", "Family", "Genus", "Species")
)
header_data <- header_labels %>%
  left_join(depth, by = "depth")
header_data$header_x[header_data$depth == 0] <- -3
header_y <- max(tree_data$y) + 2
header_data$y <- header_y
#plot
ggtree(taxonomy_phylo, ladderize = TRUE) +
  geom_text(data = nonroot_tree_data, size = 3, vjust = -0.5,
            aes(x = mid_x, y = y, label = display_label, fontface = label_face)) +
  geom_segment(data = root_label_data, aes(x = x, xend = xend, y = y, yend = yend)) + 
  geom_text(data = root_label_data, aes(x = x + 3, y = y, label = "Plantae"),
            size = 3, vjust = -0.5) +
  geom_text(data = header_data, aes(x = header_x, y = y, label = rank_name),
            size = 4, fontface = "bold") +
  theme_tree2() +
  xlim(-8, 100) +
  coord_cartesian(clip = "off")
#Figure 1



#rda on norm intensities
#analyzes peaks, can look for peaks associated w/ env variables
#need samples as rows, peaks and parameters as columns
#peaks same input as bray, just need to remove precip from metadata
#columns wanted: tf volume, DOC and TDN conc, DBH,
#and distance, sine, and cosine for buildings, roads, and trains
rda_metadata <- as.data.frame(samp_meta[-(54:58), c(56, 16, 19, 23, 44, 46, 47, 48, 50, 51, 52, 54, 55)])
#fix that to be names instead of numbers

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
par(mfrow = c(1, 1))
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

#significance testing
dbrda_signif$CCA$tot.chi / dbrda_signif$tot.chi
dbrda_signif$CA$tot.chi / dbrda_signif$tot.chi
#27.5% of variance constrained by these 5 variables, 72.5% unconstrained
anova(dbrda_signif, permutations = 999)
anova(dbrda_signif, by = "terms", permutations = 999)
#DOC and TDN significant at P<.001 level
#tf volume and dbh significant at P<.05 level
#building distance almost significant but usually 0.5<P<0.1
anova(dbrda_signif, by = "axis", permutations = 999)
#first two axes are significant at P<.001 level, others not significant
summary(dbrda_signif)

#also want to see how chemistry plots in dbrda
chem_for_dbrda_bulk <- FTICRMS_means_t[-(54:58), c(1, 2, 3, 5, 10, 11, 12)]
arrows_chem_bulk <- envfit(dbrda_signif, chem_for_dbrda_bulk, permutations = 999, na.rm = TRUE)
chem_for_dbrda_classes <- FTICRMS_means_t[-(54:58), c(19:39)]
arrows_chem_classes <- envfit(dbrda_signif, chem_for_dbrda_classes, permutations = 999, na.rm = TRUE)
#dbrda_chem_bulk <- dbrda(peak_dist_matrix_all ~., data = chem_for_dbrda_bulk, na.action = na.exclude)
#dbrda_chem_classes <- dbrda(peak_dist_matrix_all ~., data = chem_for_dbrda_classes, na.action = na.exclude)
#plot(dbrda_chem_bulk)

#want to project the rainwater samples into the same space
peak_dist_w_precip <- vegdist(norm_intens_t_all, method = "bray")
peak_dist_matrix_w_precip <- as.matrix(peak_dist_w_precip)
tf_indices <- which(all_samples %in% tf_samples_w_DOC)
positions_all <- wcmdscale(peak_dist_w_precip, k = 2, eig = TRUE)
positions_tf <- positions_all$points[-(54:58),]
positions_precip <- positions_all$points[(54:58),]
positions_tf_w_DOC <- positions_all$points[tf_indices,]
scores_tf <- scores(dbrda_signif, display = "sites", choices = c(1, 2))
scores_tf_w_DOC <- scores_tf[tf_samples_w_DOC,]
procr_tf <- procrustes(scores_tf_w_DOC, positions_tf_w_DOC)
precip_in_dbrda <- predict(procr_tf, positions_precip)

plot(dbrda_signif, type = "n")
text(dbrda_signif, display = "bp", col = "blue", cex = 0.8)
points(dbrda_signif, display = "sites", pch = 20, col = "black", cex = 0.8)
points(precip_in_dbrda, pch = 17, col = "red", cex = 0.8)

#adding labels
plot(dbrda_signif, type = "n", choices = c(1, 2), scaling = -1, xlim = c(-3, 3), ylim = c(-2, 2))
text(dbrda_signif, display = "bp", col = "blue", cex = 0.8)
points(positions_all, pch = 17, col = "darkgray", cex = 0.8)
points(dbrda_signif, display = "sites", pch = 20, col = "black", cex = 0.8)
points(precip_in_dbrda, pch = 17, col = "red", cex = 0.8)
text(precip_in_dbrda, display = "sites", labels = rownames(precip_in_dbrda), pos = 3, cex = 0.6, offset = 0.3)
legend("topright",
       legend = c("Throughfall", "Rain"),
       col = c("black", "red"), pch = c(20, 17),
       cex = 0.8, bty = "n", inset = c(0.08, 0), y.intersp = 0.6)

#want to see if conifers/magnolias cluster uniquely
conifers <- samp_meta %>% 
  filter(Class =="Pinopsida") %>% pull(Sample_ID)
positions_conifer <- positions_all$points[conifers, ]
magnolias <- samp_meta %>% 
  filter(Order == "Magnoliales") %>% pull(Sample_ID)
positions_magnolia <- positions_all$points[magnolias, ]
other_tf <- samp_meta %>% 
  filter(Class !="Pinopsida" & Order != "Magnoliales" & Sample_type != "Precipitation") %>% pull(Sample_ID)
positions_other <- positions_all$points[other_tf, ]

order_groups <- ifelse(tf_meta$Class == "Pinopsida", "Pinopsida",
                       ifelse(tf_meta$Order == "Magnoliales", "Magnoliales",
                              "Other"))
order_colors <- c("Pinopsida" = "forestgreen", "Magnoliales" = "purple", "Other" = "gray40")
order_shapes <- c("Pinopsida" = 22, "Magnoliales" = 25, "Other" = 21)
plot(dbrda_signif, type = "n", xlim = c(-3, 3), ylim = c(-2., 1.9))
text(dbrda_signif, display = "bp", col = "blue", cex = 0.8, lwd = 1.5)
#plot(arrows_chem_bulk, col = "darkred", cex = 0.8, lwd = 1.5)
for(group in c("Other", "Pinopsida", "Magnoliales")) {
  indices <- which(order_groups == group)
  points(dbrda_signif, display = "sites", select = indices,
         col = order_colors[group], bg = order_colors[group],
         pch = order_shapes[group], cex = 0.8)
}
points(precip_in_dbrda, pch = 17, col = "red", cex = 0.8)
legend("topright",
       legend = c("Pinopsida", "Magnoliales", "Other", "Rain"),
       pch = c(22, 25, 21, 17),
       col = c("forestgreen", "purple", "gray40", "red"),
       pt.bg = c("forestgreen", "purple", "gray40", "red"),
       cex = 0.8, bty = "n", inset = c(0.15, 0), y.intersp = 0.6)
#isn't using my x and y limits, redrawing manually to make it use them

site_scores_signif <- scores(dbrda_signif, display = "sites", choices = c(1, 2))
bp_scores_signif <- scores(dbrda_signif, display = "bp", choices = c(1, 2))
label_offsets_dbrda <- data.frame(
  x = c(0.14, -0.1, 0.1, -0.27, -0.45),
  y = c(0.12, -0.1, -0.1, 0, -0.015)
)
rownames(label_offsets_dbrda) <- rownames(bp_scores_signif)

#final version
pdf("dbrda_plot.pdf", width = 8, height = 8)
plot(NA, xlim = c(-1.5, 3), ylim = c(-2.2, 1), xaxs = "i", yaxs = "i",
     xlab = "dbRDA1", ylab = "dbRDA2")
abline(h = 0, v = 0, lty = 3, col = "gray50")
for(group in c("Other", "Pinopsida", "Magnoliales")) {
  indices <- which(order_groups == group)
  points(site_scores_signif[indices, ],
         col = order_colors[group], bg = order_colors[group],
         pch = order_shapes[group], cex = 1.0)
}
arrows(0, 0, bp_scores_signif[, 1], bp_scores_signif[, 2],
       col = "blue", length = 0.1, lwd = 1.5)
text(bp_scores_signif[, 1] + label_offsets_dbrda$x, bp_scores_signif[, 2] + label_offsets_dbrda$y,
     labels = c("TF volume", "[DOC]", "[TDN]", "Tree size", "Building distance"), col = "blue", cex = 1.0)
points(precip_in_dbrda, pch = 17, col = "red", cex = 1.0)
legend("topleft",
       legend = c("Pinopsida TF", "Magnoliales TF", "Other TF", "Rainwater"),
       pch = c(22, 25, 21, 17),
       col = c("forestgreen", "purple", "gray40", "red"),
       pt.bg = c("forestgreen", "purple", "gray40", "red"),
       cex = 1.0, bty = "n", inset = c(0.02, -0.01), y.intersp = 1.4)
dev.off()
#this is final dbrda, fig 2

#writing dbrda scores to csv for inclusion in Table S4
write.csv(site_scores_signif, "site_scores_signif.csv", row.names = TRUE)
write.csv(precip_in_dbrda, "precip_in_dbrda.csv", row.names = TRUE)

#testing for correlations between explanatory variables and dbrda axes 1 and 2
cor.test(site_scores_signif[, 1], rda_meta_signif$Volume, method = "pearson")
cor.test(site_scores_signif[, 2], rda_meta_signif$Volume, method = "pearson")
cor.test(site_scores_signif[, 1], rda_meta_signif$`[DOC]`, method = "pearson")
cor.test(site_scores_signif[, 2], rda_meta_signif$`[DOC]`, method = "pearson")
cor.test(site_scores_signif[, 1], rda_meta_signif$`[TDN]`, method = "pearson")
cor.test(site_scores_signif[, 2], rda_meta_signif$`[TDN]`, method = "pearson")
cor.test(site_scores_signif[, 1], rda_meta_signif$`Tree size`, method = "pearson")
cor.test(site_scores_signif[, 2], rda_meta_signif$`Tree size`, method = "pearson")
cor.test(site_scores_signif[, 1], rda_meta_signif$`Building distance`, method = "pearson")
cor.test(site_scores_signif[, 2], rda_meta_signif$`Building distance`, method = "pearson")
vars_dbrda <- colnames(rda_meta_signif)
for (v in vars_dbrda) {
  for (axis in 1:2) {
    result <- cor.test(site_scores_signif[, axis], rda_meta_signif[[gsub("`", "", v)]])
    cat(v, "- dbRDA", axis, ": r =", round(result$estimate, 3),
        ", p =", round(result$p.value, 4), "\n")
  }
}


#adding density ellipses around the four groups
#ended up not using these as they were visually cluttering
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
#still not working, I don't think they improve the figure anyway


#get more statistics on explanatory power of the dbrda
precip_meta_rda_clean <- precip_meta[, c("Tf_volume_mL", "DOC_uM", "TDN_uM", "DBH_cm", "BLDG_dist")]
colnames(precip_meta_rda_clean) <- c("Volume", "[DOC]", "[TDN]", "Tree size", "Building distance")
precip_meta_rda_clean$`Tree size` <- as.numeric(precip_meta_rda_clean$`Tree size`)
precip_meta_rda_clean$`Building distance` <- as.numeric(precip_meta_rda_clean$`Building distance`)
precip_scores <- predict(dbrda_signif, newdata = precip_meta_rda_clean, type = "lc", na.action = na.exclude)

#test individual environmental variables for significance and variation explained
no_NA_vars <- colnames(rda_metadata[,colnames(rda_metadata) != "DOC_mgL"])
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
#this is showing only TDN being significant on its own, DBH is close (.062)
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



#test for correlations between variables
cor_matrix <- cor(rda_metadata, use = "complete.obs")
cor_matrix
#this says only variable pair with strong correlation is train sin and cos (-0.95)
#the rest are 0.43 or less
corrplot(cor_matrix, method = "color", type = "upper",
         order = "hclust", tl.cex = 0.8, addCoef.col = "black", number.cex = 0.6)




#analyzing taxonomy data for significance as another dbrda
#excluding phylum because it's no broader than class
#there are 2 different classes, 12 orders, 16 families, 24 genera, and 36 species
phylo_meta <- as.data.frame(samp_meta[-(54:58), (5:9)], as.factor)
rownames(phylo_meta) <- tf_samples
phylo_vars <- colnames(phylo_meta)
phylo_signif <- list()
for(var in phylo_vars) {
  formula_ph <- as.formula(paste("peak_dist_matrix_all ~", var))
  single_dbrda_ph <- dbrda(formula_ph, data = phylo_meta)
  sig_test_ph <- anova(single_dbrda_ph, permutations = 999)
  var_explained_ph <- single_dbrda_ph$CCA$tot.chi / single_dbrda_ph$tot.chi
  phylo_signif[[var]] <- list(
    variable = var,
    variance_explained = var_explained_ph,
    p_value = sig_test_ph$`Pr(>F)`[1],
    F_statistic = sig_test_ph$F[1],
    df = sig_test_ph$Df[1]
  )
}
phylo_signif_results <- do.call(rbind, lapply(phylo_signif, data.frame))
phylo_signif_results <- phylo_signif_results[order(phylo_signif_results$variance_explained, decreasing = TRUE), ]
print(phylo_signif_results)

norm_intens_t_meta <- cbind(phylo_meta, norm_intens_t[-(54:58),])
write.csv(norm_intens_t_meta, "norm_intens_t_meta.csv", row.names = TRUE)

#permanova on phylogeny data- this shows the same results as dbrda above
for(var in phylo_vars) {
  formula <- as.formula(paste("peak_dist_matrix_all ~", var))
  permanova_result <- adonis2(formula, data = phylo_meta, permutations = 999)
  print(paste("Variable:", var))
  print(permanova_result)
  cat("\n")
}
#there's a good linear relationship between # of categories in a taxonomic rank and % of variance explained

#now look at them sequentially instead of individually
phylo_meta[phylo_vars] <- lapply(phylo_meta[phylo_vars], as.factor)
sequential_phylo <- list()
remaining_var <- 1.0
for(i in 1:5) {
  current_var <- phylo_vars[i]
  if(i == 1) {
    formula_phylo <- as.formula(paste("peak_dist_matrix_all ~", current_var))
    dbrda_phylo <- dbrda(formula_phylo, data = phylo_meta)
  } else {
    conditioning_vars <- paste(phylo_vars[1:(i-1)], collapse = " + ")
    formula_phylo <- as.formula(paste("peak_dist_matrix_all ~", current_var, "+ Condition(", conditioning_vars, ")"))
    dbrda_phylo <- dbrda(formula_phylo, data = phylo_meta)
  }
  if(i == 1) {
    var_explained_total <- dbrda_phylo$CCA$tot.chi / dbrda_phylo$tot.chi
    var_explained_remaining <- var_explained_total
  } else {
    var_explained_total <- dbrda_phylo$CCA$tot.chi / dbrda_phylo$tot.chi
    var_explained_remaining <- var_explained_total / remaining_var
  }
  remaining_var <- remaining_var - var_explained_total
  sig_test_phylo <- anova(dbrda_phylo, permutations = 999)
  sequential_phylo[[current_var]] <- list(
    variable = current_var,
    variance_explained_total = var_explained_total,
    variance_explained_remaining = var_explained_remaining,
    remaining_variance_after = remaining_var,
    p_value = sig_test$`Pr(>F)`[1],
    F_statistic = sig_test_phylo$F[1]
  )
  cat("Level:", current_var, "\n")
  cat("Variance explained (of total):", round(var_explained_total, 4), "\n")
  cat("Variance explained (of remaining):", round(var_explained_remaining, 4), "\n")
  cat("P-value:", sig_test$`Pr(>F)`[1], "\n")
  cat("Remaining variance:", round(remaining_var, 4), "\n\n")
}
sequential_df <- do.call(rbind, lapply(sequential_phylo, data.frame))
print(sequential_df)
#this shows the same results as doing them independently
#also the p values are odd, all showing the same (0.705)

#restricted permutations to get meaningful p values, though with reduced statistical power
seq_phylo_restr <- list()
for (i in 1:5) {
  current_var_rest <- phylo_vars[i]
  if(i == 1) {
    formula_rest <- as.formula(paste("peak_dist_matrix_all ~", current_var_rest))
    dbrda_rest <- dbrda(formula_rest, data = phylo_meta)
    sig_test_rest <- anova(dbrda_rest, permutations = 999)
  } else {
    conditioning_vars_rest <- paste(phylo_vars[1:(i-1)], collapse = " + ")
    formula_rest <- as.formula(paste("peak_dist_matrix_all ~", current_var_rest, "+ Condition(", conditioning_vars_rest, ")"))
    dbrda_rest <- dbrda(formula_rest, data = phylo_meta)
    if(i == 2) {
      perm_design <- how(within = Within(type = "free"),
                         blocks = phylo_meta[[phylo_vars[1]]])
    } else {
      perm_design <- how(within = Within(type = "free"),
                         blocks = phylo_meta[[phylo_vars[i-1]]])
    }
    sig_test_rest <- anova(dbrda_rest, permutations = perm_design)
  }
  var_explained_total_rest <- dbrda_rest$CCA$tot.chi / dbrda_rest$tot.chi
  seq_phylo_restr[[current_var_rest]] <- list(
    variable = current_var_rest,
    variance_explained = var_explained_total_rest,
    p_value = sig_test_rest$`Pr(>F)`[1],
    F_statistic <- sig_test_rest$F[1]
  )
  cat("Level:", current_var_rest, "\n")
  cat("Variance explained:", round(var_explained_total_rest, 4), "\n")
  cat("P-value:", sig_test_rest$`Pr(>F)`[1], "\n\n")
}
#results look good


#analysis of env variables on residuals after phylogeny
species_formula <- as.formula(paste("peak_dist_matrix_all ~", paste(phylo_meta, collapse = " + ")))
dbrda_species <- dbrda(species_formula, data = phylo_meta)
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
#again sequential percentages are just separating the components of looking at them individually


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
#end