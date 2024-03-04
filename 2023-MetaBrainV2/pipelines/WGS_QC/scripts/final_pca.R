# Author: Urmo Võsa
# Edited by Joost Bakker

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

# Function that reads in fam file
read_fam <- function(path) {
  NAMES.FAM <- c("family.ID", "sample.ID", "paternal.ID",
                 "maternal.ID", "sex", "affection")

  pattern <- "\\.bed$"
  
  if (!grepl(pattern, path)) {
    famfile <- paste0(path, ".fam")
  } else {
    famfile <- sub(pattern, ".fam", path)
  } 

  fam <- bigreadr::fread2(famfile, col.names = NAMES.FAM, keepLeadingZeros = TRUE,
                          colClasses = list(character = c(1,2)), nThread = 1)

  return(fam)
}

# Create command line argument list
option_list <- list(
    make_option(c("--target_bed"), type = "character",
    help = "Path to the target genotype file (bed/bim/fam format). Required file extension: .bed.")
    )

# Create arg parser and pars args
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Read in bed file
bed_qc <- bed(args$target_bed)
bed_qc$.fam <- read_fam(args$target_bed)

# Do PCA
target_pca_qcd <- bed_autoSVD(bed_qc, k = 10, ncores = 4)

# Visualise loadings
plot(target_pca_qcd, type = "loadings", loadings = 1:10, coeff = 0.6)
ggsave("Target_PCs_postQC_Loadings.png", type = "cairo", height = (5 * 7) * 0.7, width = (5 * 7) * 0.7, units = "in", dpi = 300)

PCsQ <- predict(target_pca_qcd)
PCsQ <- as.data.frame(PCsQ)

colnames(PCsQ) <- paste0("PC", 1:10)
rownames(PCsQ) <- bed_qc$fam$sample.ID

# Create plots of the loadings
p1 <- ggplot(PCsQ, aes(x = PC1, y = PC2)) + theme_bw() + geom_point(alpha = 0.5)
p2 <- ggplot(PCsQ, aes(x = PC3, y = PC4)) + theme_bw() + geom_point(alpha = 0.5)
p3 <- ggplot(PCsQ, aes(x = PC5, y = PC6)) + theme_bw() + geom_point(alpha = 0.5)
p4 <- ggplot(PCsQ, aes(x = PC7, y = PC8)) + theme_bw() + geom_point(alpha = 0.5)
p5 <- ggplot(PCsQ, aes(x = PC9, y = PC10)) + theme_bw() + geom_point(alpha = 0.5)
p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

# Save plots
ggsave("Target_PCs_postQC.png", type = "cairo", height = 10 * 1.5, width = 9 * 1.3, units = "in", dpi = 300)
ggsave("Target_PCs_postQC.pdf", height = 10 * 1.5, width = 9 * 1.3, units = "in", dpi = 300)

# Write PCs to text file
fwrite(PCsQ, "GenotypePCs.txt", row.names = TRUE, sep = "\t", quote = FALSE)

# Create scree plots
singlar_value <- data.frame(
  PC = paste0("PC", 1:10), 
  sv = target_pca_qcd$d
  )
singlar_value$PC <- factor(singlar_value$PC, levels = as.character(singlar_value$PC))

p <- ggplot(singlar_value, aes(x = PC, y = sv)) + 
geom_bar(stat = "identity") + 
theme_bw() + 
ylab("Singular value")

# Save scree plots
ggsave("Target_PCs_scree_postQC.png", type = "cairo", height = 5, width = 9, units = "in", dpi = 300)
ggsave("Target_PCs_scree_postQC.pdf", height = 5, width = 9, units = "in", dpi = 300)
