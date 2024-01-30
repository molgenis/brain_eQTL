library(optparse)
library(data.table)

# Create command line argument list
option_list <- list(
    make_option(c("--sample_file"), type = "character",
    help = "Path to the sample file that should be shuffled.")
    )

# Create arg parser and pars args
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)


# Read in sample file
samples <- read.delim(args$sample_file, header=FALSE)

length(samples)

# Shuffle sample order
rows <- sample(nrow(samples))
shuffled_samples <- samples[rows, ]



# Write shuffled samples to file
fileConn<-file("ShuffledSampleOrder.txt")
writeLines(shuffled_samples, fileConn)
close(fileConn)
