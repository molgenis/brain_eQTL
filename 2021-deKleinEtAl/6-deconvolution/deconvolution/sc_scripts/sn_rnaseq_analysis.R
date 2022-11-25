# Title     : SN RNA-seq analysis
# Objective : -
# Created by: mvochteloo
# Created on: 22/02/2021

# install.packages('Seurat')
library(Seurat)
library(Matrix)
library(ggplot2)


work_dir <- '/groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/2020-10-22-MVochteloo-Copy/'
list.files(work_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = work_dir)

projid_table <- read.table('/groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/filtered_column_metadata.txt.gz', header = T, stringsAsFactors = F)
projid_table$TAG <- gsub('\\.', '-', projid_table$TAG)

rownames(projid_table) <- projid_table$TAG
projid_table$TAG <- NULL

seurat_object <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                            min.cells = 0,
                                            min.features = 0,
                                            project = "scRNAseq",
                                            meta.data = projid_table)

filtered_gene_row_file <- read.csv("/groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/filtered_gene_row_names.txt.gz", sep = "\t")
filtered_gene_row_tags <- filtered_gene_row_file[, 1]
length(filtered_gene_row_tags)
# [1] 17925

seurat_object <- seurat_object[filtered_gene_row_tags, !is.na(seurat_object@meta.data$projid)]
seurat_object
# An object of class Seurat
# 17352 features across 70634 samples within 1 assay
# Active assay: RNA (17352 features, 0 variable features)

seurat_object <- SCTransform(seurat_object)

length(rownames(seurat_object))
# [1] 16866 = genes
length(colnames(seurat_object))
# [1] 70634 = cells

Idents(object=seurat_object) <- "broad.cell.type"

mean_expression_per_cell_type <- function(seurat, sample.id.column.name = "projid", gene_trans_df = NULL) {
  individuals <- unique(seurat@meta.data[, sample.id.column.name])
  idents <- unique(seurat@active.ident)

  cell_counts <- matrix(nrow = length(idents), ncol = length(individuals),
                        dimnames = list(idents, individuals))

  for (ident in idents) {

    cells_cell_type <- seurat[, seurat@active.ident == ident]

    mean_expression_matrix <- matrix(nrow = nrow(cells_cell_type), ncol = length(individuals),
                                     dimnames = list(rownames(cells_cell_type), individuals))

    for (individual in individuals) {
      if (sum(cells_cell_type@meta.data[, sample.id.column.name] == individual) == 0) {
        mean_expression_matrix[, toString(individual)] <- 0
        cell_counts[ident, toString(individual)] <- 0
      } else {
        cells_cell_type_individual <- cells_cell_type[, cells_cell_type@meta.data[, sample.id.column.name] == individual]

        mean_expression_matrix[, toString(individual)] <- rowMeans(cells_cell_type_individual)
        cell_counts[ident, toString(individual)] <- ncol(cells_cell_type_individual)
      }
      print(cell_counts[ident, toString(individual)])
    }

    if (!is.null(gene_trans_df)) {
      rownames(mean_expression_matrix) <- gene_trans_df[match(rownames(mean_expression_matrix), gene_trans_df[,2]), 1]
      mean_expression_matrix <- mean_expression_matrix[!is.na(rownames(mean_expression_matrix)), ]
    }

    write.table(mean_expression_matrix,
                file = paste0(work_dir, ident, "_expression", ".tsv"),
                quote = F, sep = "\t", col.names = NA)

  }

  write.table(cell_counts,
              file = paste0(work_dir, "cell_counts.txt"),
              quote = F, sep = "\t", col.names = NA)
}

mean_expression_per_cell_type(seurat_object, gene_trans_df = gene_trans_df)