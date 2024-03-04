# Author: Urmo Võsa
# Edited by Joost Bakker

library(ggplot2)
library(data.table)
library(optparse)

# Create command line argument list
option_list <- list(
    make_option(c("--ref_afs"), type = "character",
    help = "Path to the reference allele frequencies file (.gz)"),
    make_option(c("--target_afs"), type = "character",
    help = "Path to the target allele frequencies file (.gz)")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Read in frequency files
reference_freq <- fread(args$ref_afs, header = TRUE)
target_freq <- fread(args$target_afs, header = TRUE)

# Merge target and reference 
af_merged <- merge(reference_freq,target_freq, by="ID")

# Create a scatter plot
ggplot(af_merged, aes(x = af_merged$ALT_FREQS.x, y = af_merged$ALT_FREQS.y)) +
  geom_point() +
  labs(title = "Reference vs Target Allele Frequencies", x = "Reference Frequency", y = "Target Frequency")

# Save plot as png
ggsave("allele_frequencies.png")