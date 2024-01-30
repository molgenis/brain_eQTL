library(bigreadr)
library(ggplot2)
library(data.table)
library(dplyr)
library(optparse)

# Create command line argument list
option_list <- list(
    make_option(c("--het_check_file"), type = "character",
    help = "Path to the plink heterozygosity check output file")
    )

# Create arg parser and pars args
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Read data from plink het check and determine failed samples
het <- fread(args$het_check_file, header = TRUE, keepLeadingZeros = TRUE, colClasses = list(character = c(1,2)))
het$het_rate <- (het$OBS_CT - het$`O(HOM)`) / het$OBS_CT
het_fail_samples <- het[het$het_rate < mean(het$het_rate) - 3 * sd(het$het_rate) | het$het_rate > mean(het$het_rate) + 3 * sd(het$het_rate), ]

# Create file with failed samples
fwrite(het_fail_samples, "HeterozygosityFailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Create plot of the heterozygosity check
p <- ggplot(het, aes(x = het_rate)) + geom_histogram(color = "#000000", fill = "#000000", alpha = 0.5) +
xlab("Heterozygosity rate") +
geom_vline(xintercept = c(mean(het$het_rate), mean(het$het_rate) + 3 * sd(het$het_rate), mean(het$het_rate) - 3 * sd(het$het_rate)), linetype = 2, colour = "red") +
theme_bw()

# Save plots as png and pdf
ggsave("HetCheck.png", type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
ggsave("HetCheck.pdf", height = 7 / 2, width = 9, units = "in", dpi = 300)
