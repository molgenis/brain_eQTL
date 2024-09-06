#!/opt/conda/envs/eQTLGenPopAssign/bin/Rscript
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

# Create command line argument list
option_list <- list(
    make_option(c("--kin_file"), type = "character",
    help = "Path to plink kinship file."),
    make_option(c("--target_bed"), type = "character",
    help = "Path to the target genotype file (bed/bim/fam format). Required file extension: .bed.")
    )

# Create arg parser and pars args
parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Load target bed file
target_bed <- bed(args$target_bed)

# Read kin file
related <- read.delim(args$kin_file)

# Remove samples that are related to each other
related$IID1 <- as.character(related$IID1)
related$IID2 <- as.character(related$IID2)

# Get list of related samples
related_individuals <- unique(c(related$IID1, related$IID2))


# Create an empty vector with samples to remove
samples_to_remove_due_to_relatedness <- c()

# If there are related individuals, remove these in the following step.
if (length(related_individuals) > 0) {

  # Define a graph wherein each relation depicts an edge between vertices (samples)
  relatedness_graph <- graph_from_edgelist(
    as.matrix(related[,c("IID1", "IID2")]),
    directed = F)

  relatedness_graph <- simplify(
    relatedness_graph,
    remove.multiple = TRUE,
    remove.loops = FALSE,
    edge.attr.comb = igraph_opt("edge.attr.comb")
  )

  # Now, get a list of samples that should be removed due to relatedness
  while (length(V(relatedness_graph)) > 1) {

    # Get the degrees (how many edges does each vertex have)
    degrees_named <- degree(relatedness_graph)

    # Get the vertex with the least amount of degrees (edges)
    least_vertex_samples <- names(degrees_named)[min(degrees_named) == degrees_named]

    curr_vertex <- least_vertex_samples[1]

    # Get all vertices that have an edge with curr_vertex
    related_vertices <- names(relatedness_graph[curr_vertex][relatedness_graph[curr_vertex] > 0])

    # Add these vertexes to the list of vertices to remove
    samples_to_remove_due_to_relatedness <- c(samples_to_remove_due_to_relatedness, related_vertices)

    # Remove the vertices to remove
    relatedness_graph <- delete_vertices(relatedness_graph, c(curr_vertex, related_vertices))
  }

  # Remove these indices from the indices that remained after the previous check
  passed_samples <- target_bed$fam$`sample.ID`[
    (!target_bed$fam$`sample.ID` %in% samples_to_remove_due_to_relatedness)]
} else{
    passed_samples <- target_bed$fam$`sample.ID`
}

# Write passed samples to file
fileConn<-file("RelatednessPassedSamples.txt")
writeLines(as.character(passed_samples), fileConn)
close(fileConn)

