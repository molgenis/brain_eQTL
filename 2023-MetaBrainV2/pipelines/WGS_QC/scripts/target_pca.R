library(bigreadr)
library(bigsnpr)
library(dplyr)
library(ggplot2)
library(data.table)
library(optparse)
library(patchwork)
library(stringr)
library(rmarkdown)
library(Cairo)
library(igraph)

# Create command line argument list
option_list <- list(
    make_option(c("--target_bed"), type = "character",
    help = "Path to the target genotype file (bed/bim/fam format). Required file extension: .bed.")
    )

# Create arg parser and pars args
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Load bed file
target_bed <- bed(args$target_bed)

# Do PCA
target_pca <- bed_autoSVD(target_bed, k = 10, ncores = 4)

# Find outlier samples
prob <- bigutilsr::prob_dist(target_pca$u, ncores = 4)
S <- prob$dist.self / sqrt(prob$dist.nn)
Sthresh <- 0.4

# Plot outliers
p <- ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = 0:5 / 5, limits = c(0, NA)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)") +
  geom_vline(aes(xintercept = 0.4), colour = "red", linetype = 2)

# Save plots as png and pdf
ggsave("PC_dist_outliers_S.png", type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
ggsave("PC_dist_outliers_S.pdf", height = 7 / 2, width = 9, units = "in", dpi = 300)

# Visualise PCs and color outliers red
PCs <- predict(target_pca)

PCs <- as.data.frame(PCs)
colnames(PCs) <- paste0("PC", 1:10)
PCs$S <- S

PCs$outlier_ind <- "no"

if (any(PCs$S > Sthresh)) {
  PCs[PCs$S > Sthresh, ]$outlier_ind <- "yes"
}
PCs$sd_outlier <- "no"
sd_outlier_selection <- ((PCs$PC1 > mean(PCs$PC1) + args$SD_threshold * sd(PCs$PC1)
  | PCs$PC1 < mean(PCs$PC1) - args$SD_threshold * sd(PCs$PC1))
  | (PCs$PC2 > mean(PCs$PC2) + args$SD_threshold * sd(PCs$PC2)
  | PCs$PC2 < mean(PCs$PC2) - args$SD_threshold * sd(PCs$PC2)))
if (any(sd_outlier_selection)) {
  PCs[sd_outlier_selection, ]$sd_outlier <- "yes"
}

PCs$outlier <- "no"
if (nrow(PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "no", ]) > 0){
    PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "no", ]$outlier <- "S outlier"}
if(nrow(PCs[PCs$outlier_ind == "no" & PCs$sd_outlier == "yes", ]) > 0){
PCs[PCs$outlier_ind == "no" & PCs$sd_outlier == "yes", ]$outlier <- "SD outlier"}
if(nrow(PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "yes", ]) > 0){
PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "yes", ]$outlier <- "S and SD outlier"
}

# For first 2 PCs also remove samples which deviate from the mean
p1 <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick")) +
geom_vline(xintercept = c(mean(PCs$PC1) + 3 * sd(PCs$PC1), mean(PCs$PC1) - 3 * sd(PCs$PC1)), colour = "firebrick", linetype = 2) +
geom_hline(yintercept = c(mean(PCs$PC2) + 3 * sd(PCs$PC2), mean(PCs$PC2) - 3 * sd(PCs$PC2)), colour = "firebrick", linetype = 2)
p2 <- ggplot(PCs, aes(x = PC3, y = PC4, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p3 <- ggplot(PCs, aes(x = PC5, y = PC6, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p4 <- ggplot(PCs, aes(x = PC7, y = PC8, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p5 <- ggplot(PCs, aes(x = PC9, y = PC10, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))

p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

ggsave("PCA_outliers.png", type = "cairo", height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)
ggsave("PCA_outliers.pdf", height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)

# Filter out related samples and outlier samples, write out QCd data
message("Filter out related samples and outlier samples, write out QCd data.")
indices_of_passed_samples <- rows_along(target_bed)
indices_of_passed_samples <- indices_of_passed_samples[PCs$outlier == "no"]
samples_to_include <- data.frame(family.ID = target_bed$.fam$`family.ID`[indices_of_passed_samples], sample.IDD2 = target_bed$.fam$sample.ID[indices_of_passed_samples])

fwrite(data.table::data.table(samples_to_include), "SamplesToInclude.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
