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

# Metadata

samp_meta <- read_csv(file = 'data/sample_metadata.csv')

# FTICRMS data

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

#transpose for NMDS
norm_transp_intensities <- t(norm_intensities)
#check if this is right
#data is good but it didn't bring row/column names with it,
#so row names (mz, formula, sample IDs) are gone and
#sequential column names have been added
colnames(norm_transp_intensities) <- norm_intensities$mz
mycolumnnames <- read.table(file = 'data/columnnames.csv', sep = ",")
mycols <- mycolumnnames[1, ]
mycolumnnames_t <- read.table(file = 'data/columnnames_t.csv', sep = ",")
mycolumnnames_t$V1
rownames(norm_transp_intensities) <- mycolumnnames_t$V1
#I think this is good
#column names are duplicate of first row (mz values),
#might be off by one with first mz over row names and second mz over first actual mz,
#but I think this is only happening when I write to a table, not in actual data

#check my data format is good
is.matrix(norm_transp_intensities)
#True, that's good
rownames(norm_transp_intensities)
#row names are mz, formula, and sample IDs- correct
colnames(norm_transp_intensities)
#mz values, it has a number paired with that number slightly rounded- seems good
norm_transp_intensities[1, ]
#same as the column names
str(norm_transp_intensities)
#107 rows by 10799 columns, values in first row are column names, values in first column are mz values- good
isSymmetric(norm_transp_intensities)
#false, that should be fine