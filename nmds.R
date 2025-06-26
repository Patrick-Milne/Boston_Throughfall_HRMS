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


# Metadata

samp_meta <- read_csv(file = 'data/sample_metadata.csv')
peak_meta <- read.csv(file = 'data/peak_metadata.csv')

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

#kevin's method
#don't want precipitation samples included
tf_samples <- samp_meta %>% 
  filter(Sample_type !="Precipitation") %>% pull(Sample_ID)
#for use later when hybrid species will be included
tf_nonhyb_samples <- samp_meta %>% 
  filter(PhyloName !="NA") %>% pull(Sample_ID)

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
#stress should be below 0.1; it's 0.06201059

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
anova(dbrda)
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
cor_matrix <- cor(rda_metadata, use = "complete.obs")
cor_matrix
#plot labeled with sample ID
plot(dbrda, type = "n")
points(dbrda, display = "sites", pch = 19, col = "blue", cex = 0.8)
text(dbrda, display = "sites", labels = rownames(rda_metadata), pos = 3, cex = 0.6, offset = 0.3)
text(dbrda, display = "bp", col = "red", cex = 0.8)

#mantel
large_phylo_dist <- as.dist(large_phylo_matrix)
#do 999 permutations for speed, 9999 for final published result
mantel <- mantel(large_phylo_dist, peak_dist,
                 method = "spearman", permutations = 999)
#r=0.01525, significance = 0.37-- not significant

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
