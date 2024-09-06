#!/opt/conda/envs/eQTLGenPopAssign/bin/Rscript
# Author: Urmo Võsa
# Edited by Joost Bakker


library(bigreadr)
library(ggplot2)
library(data.table)
library(dplyr)
library(optparse)

# Create command line argument list
option_list <- list(
    make_option(c("--sex_check_file"), type = "character",
    help = "Sex check file from plink")
    )

# Create arg parser and pars args
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

sexcheck <- fread(args$sex_check_file, keepLeadingZeros = TRUE,
                colClasses = list(character = c(1,2)))

# Annotate samples who have clear sex
sexcheck$F_PASS <- !(sexcheck$F > 0.4 & sexcheck$F < 0.6)
sexcheck$MATCH_PASS <- case_when(sexcheck$PEDSEX == 0 ~ T,
                                sexcheck$STATUS == "PROBLEM" ~ F,
                                TRUE ~ T)
sexcheck$PASS <- sexcheck$MATCH_PASS & sexcheck$F_PASS

# Create sex check plot
sex_cols <- c("0" = "black", "1" = "orange", "2" = "blue")
p <- ggplot(sexcheck, aes(x = F, fill = factor(PEDSEX))) +
geom_histogram(position="stack", color = "black", alpha = 0.5) +
scale_fill_manual(values = sex_cols, breaks = c("0", "1", "2"), labels = c("Unknown", "Male", "Female"), name = "Reported sex") +
geom_vline(xintercept = c(0.4, 0.6), colour = "red", linetype = 2) + theme_bw()

# Save plots as png and pdf
ggsave("SexCheck.png", p, type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
ggsave("SexCheck.pdf", p, height = 7 / 2, width = 9, units = "in", dpi = 300)

# Write out sex check results and failed samples
fwrite(sexcheck, "SexCheck.txt", sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(sexcheck[!sexcheck$PASS,], "SexCheckFailed.txt", sep = "\t", quote = FALSE, row.names = FALSE)
