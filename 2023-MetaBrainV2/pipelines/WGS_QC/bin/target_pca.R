#!/opt/conda/envs/eQTLGenPopAssign/bin/Rscript

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
    help = "Path to the target genotype file (bed/bim/fam format). Required file extension: .bed."),
    make_option(c("--population"), type = "character",
    help = "Super population label (EUR/AMR/AFR/SAS/EAS)"),
    make_option(c("--outlier_threshold"), type = "integer",
    help = "Z-score threshold for filtering outliers")
    )

# Create arg parser and pars args
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Load bed file
target_bed <- bed(args$target_bed)

# Do PCA
target_pca <- bed_autoSVD(target_bed, k = 10, ncores = 4)

# Visualise PCs and color outliers red
PCs <- predict(target_pca)

PCs <- as.data.frame(PCs)
colnames(PCs) <- paste0("PC", 1:10)


# Identify standard deviation outliers
PCs$sd_outlier <- "no"
sd_outlier_selection <- ((PCs$PC1 > mean(PCs$PC1) + args$outlier_threshold * sd(PCs$PC1)
  | PCs$PC1 < mean(PCs$PC1) - args$outlier_threshold * sd(PCs$PC1))
  | (PCs$PC2 > mean(PCs$PC2) + args$outlier_threshold * sd(PCs$PC2)
  | PCs$PC2 < mean(PCs$PC2) - args$outlier_threshold * sd(PCs$PC2)))

if (any(sd_outlier_selection)) {
  PCs[sd_outlier_selection, ]$sd_outlier <- "yes"
}

PCs$outlier <- "no"

if(nrow(PCs[PCs$sd_outlier == "yes", ]) > 0){
  PCs[PCs$sd_outlier == "yes", ]$outlier <- "SD outlier"
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

ggsave(paste(args$population, "PCA_outliers.png", sep="_"), type = "cairo", height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)
ggsave(paste(args$population, "PCA_outliers.pdf", sep="_"), height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)

# Filter out related samples and outlier samples, write out QCd data
message("Filter out related samples and outlier samples, write out QCd data.")
indices_of_passed_samples <- rows_along(target_bed)
indices_of_passed_samples <- indices_of_passed_samples[PCs$outlier == "no"]
samples_to_include <- data.frame(family.ID = target_bed$.fam$`family.ID`[indices_of_passed_samples], sample.IDD2 = target_bed$.fam$sample.ID[indices_of_passed_samples])

indices_of_outlier_samples <- indices_of_passed_samples[PCs$sd_outlier == "yes"]
samples_to_exclude <- data.frame(family.ID = target_bed$.fam$`family.ID`[indices_of_outlier_samples], sample.IDD2 = target_bed$.fam$sample.ID[indices_of_outlier_samples])


fwrite(data.table::data.table(samples_to_include), paste(args$population, "SamplesToInclude.txt", sep="."), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
fwrite(data.table::data.table(samples_to_exclude), paste(args$population, "OutlierSamples.txt", sep="."), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
