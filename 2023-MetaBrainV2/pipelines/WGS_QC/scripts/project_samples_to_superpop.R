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
    make_option(c("--ref_bed"), type = "character",
    help = "Path to the reference genotype file (bed/bim/fam format). Required file extension: .bed."),
    make_option(c("--target_bed"), type = "character",
    help = "Path to the target genotype file (bed/bim/fam format). Required file extension: .bed."),
    make_option(c("--ref_pops"), type = "character",
    help = "Path to the file indicating the population for each sample in reference data.")
    )

# Create arg parser and pars args
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Load reference and target bed files
ref_bed <- bed(args$ref_bed)
target_bed <- bed(args$target_bed)

# Do PCA on combined reference and target data
proj_PCA <- bed_projectPCA(
    obj.bed.ref = ref_bed,
    obj.bed.new = target_bed,
    k = 10,
    strand_flip = TRUE,
    join_by_pos = TRUE,
    match.min.prop = 0.01,
    build.new = "b38",
    build.ref = "b38",
    verbose = TRUE,
    ncores = 4
  )

# Create plots for first 10 PCs
abi <- as.data.frame(proj_PCA$OADP_proj)
colnames(abi) <- paste0("PC", 1:10)

PCs_ref <- predict(proj_PCA$obj.svd.ref)
abi2 <- as.data.frame(PCs_ref)
colnames(abi2) <- paste0("PC", 1:10)

abi2$sample <- ref_bed$fam$`sample.ID`
abi2 <- abi2[, c(11, 1:10)]

pops <- fread(args$ref_pops, keepLeadingZeros = TRUE, colClasses = list(character = c(2, 6, 7)))
abi2 <- merge(abi2, pops, by.x = "sample", by.y = "SampleID")
abi2 <- abi2[, c(1, 12, 13, 2:11)]

abi <- data.frame(sample = target_bed$fam$`sample.ID`,
                  Population = "Target", Superpopulation = "Target", abi)

abi$type <- "Target"
abi2$type <- "1000G"

combined <- rbind(abi, abi2)

combined$Superpopulation <- factor(combined$Superpopulation, levels = c("Target", "EUR", "EAS", "AMR", "SAS", "AFR"))

p00 <- ggplot(combined, aes(x = PC1, y = PC2, alpha = type)) +
geom_point() +
theme_bw() +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0)) +
ggtitle("Target sample projections\nin 1000G PC space")

combined_h <- combined[combined$Superpopulation == "Target", ]

p0 <- ggplot(combined_h, aes(x = PC1, y = PC2)) +
geom_point() +
theme_bw() +
ggtitle("Target sample projections\nzoomed in")

p1 <- ggplot(combined, aes(x = PC1, y = PC2, colour = Superpopulation, alpha = type)) +
geom_point() +
theme_bw() +
scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p2 <- ggplot(combined, aes(x = PC3, y = PC4, colour = Superpopulation, alpha = type)) +
geom_point() + theme_bw() +
scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p3 <- ggplot(combined, aes(x = PC5, y = PC6, colour = Superpopulation, alpha = type)) +
geom_point() + theme_bw() +
scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p4 <- ggplot(combined, aes(x = PC7, y = PC8, colour = Superpopulation, alpha = type)) +
  geom_point() + theme_bw() +
  scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
  "EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
  scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p5 <- ggplot(combined, aes(x = PC9, y = PC10, colour = Superpopulation, alpha = type)) +
  geom_point() + theme_bw() +
  scale_color_manual(values = c("Target" = "black", "EUR" = "blue",
  "EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
  scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p <- p00 + p0 + p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 4)

# Save plots as png and pdf
ggsave("SamplesPCsProjectedTo1000G.png", type = "cairo", height = 20, width = 9.5 * 1.6, units = "in", dpi = 300)
ggsave("SamplesPCsProjectedTo1000G.pdf", height = 20, width = 9.5 * 1.6, units = "in", dpi = 300)
fwrite(abi[, -c(2, 3, ncol(abi))], "1000G_PC_projections.txt", sep = "\t", quote = FALSE )

## Assign each sample to the superpopulation
message("Assign each sample to 1000G superpopulation.")
### Calculate distance of each sample to all samples per each population
target_samples <- abi[, -c(2, 3, ncol(abi))]

#### Use 3 PCs
target_samples <- target_samples[, c(1:4)]
rownames(target_samples) <- target_samples$sample
target_samples <- target_samples[, -1]

population_assign_res <- data.frame(sample = rownames(target_samples), abi = rownames(target_samples))

#### EUR
    for(population in c("EUR", "EAS", "AMR", "SAS", "AFR")){
    abi_e <- abi2[abi2$Superpopulation == population, ]
    head(abi_e)

    sup_pop_samples <- abi_e[, -c(2, 3, ncol(abi_e))]

    sup_pop_samples <- sup_pop_samples[, c(1:4)]
    rownames(sup_pop_samples) <- sup_pop_samples$sample
    sup_pop_samples <- sup_pop_samples[, -1]

    head(sup_pop_samples)

    comb <- rbind(target_samples, sup_pop_samples)
    head(comb)
    distance <- as.matrix(dist(comb, method = "euclidean"))
    head(distance)
    distance <- distance[c(1:nrow(target_samples)), -c(1:nrow(target_samples))]
    head(distance)

    head(rowMeans(distance))

    distance <- data.frame(sample = rownames(target_samples), MeanDistance = rowMeans(distance))
    colnames(distance)[2] <- population

    population_assign_res <- cbind(population_assign_res, distance[, -1])


}

colnames(population_assign_res)[3:ncol(population_assign_res)] <- c("EUR", "EAS", "AMR", "SAS", "AFR")
fwrite(population_assign_res[, -1], "PopAssignResults.txt", sep = "\t", quote = FALSE )

