# Title     : create_seurat_object.R
# Objective :
# Created by: mvochteloo
# Created on: 22/10/2020

# https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
# https://satijalab.org/seurat/v3.2/sctransform_vignette.html

# install.packages('Seurat')
library(Seurat)
library(Matrix)
library(ggplot2)

work_dir <- '/groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/2020-10-22-MVochteloo-Copy/'
list.files(work_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = work_dir)

projid_table <- read.table('/groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/filtered_column_metadata.txt.gz', header = T, stringsAsFactors = F)
projid_table$TAG <- gsub('\\.', '-', projid_table$TAG)

# projid_rosmapid_translate_table <- read.csv("/groups/umcg-biogen/tmp01/input/ROSMAP-scRNAseq/meta/ROSMAP_IDkey.csv", sep = ",")
# phenotype_table <- read.csv("/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt", sep="\t")
# meta_table <- merge(projid_table, merge(projid_rosmapid_translate_table, phenotype_table, by.x="rnaseq_id", by.y="rnaseq_id"), by.x="projid", by.y="projid", all.x = TRUE)
meta_table <- projid_table

rownames(meta_table) <- meta_table$TAG
meta_table$TAG <- NULL

seurat_object <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                            min.cells = 0,
                                            min.features = 0,
                                            project = "scRNAseq",
                                            meta.data = meta_table)

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

#################################################################################

seurat_object <- RunPCA(seurat_object)
seurat_object <- RunTSNE(object = seurat_object, reduction = "pca",  dims = 1:20)
seurat_object <- RunUMAP(object = seurat_object, reduction = "pca",  dims = 1:20)
seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:50)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

for (name in c("pre.cluster", "broad.cell.type", "Subcluster", "seurat_clusters")) {
  DimPlot(seurat_object, reduction = 'tsne', group.by = name)
  ggsave(paste0("ROSMAP_scRNAseq_tsne_groupedBy", name ,".pdf"), width = 10, height = 10)
}

visualisation_table <- merge(merge(merge(seurat_object@meta.data, seurat_object[['pca']]@cell.embeddings, by="row.names", all=TRUE), seurat_object[['tsne']]@cell.embeddings, by="row.names", all=TRUE), seurat_object[['umap']]@cell.embeddings, by="row.names", all=TRUE)
visualisation_table$Row.names <- NULL
visualisation_table$Row.names <- NULL
visualisation_table$Row.names <- NULL
write.table(visualisation_table,
            file = paste0(work_dir, "2020-12-01-metadata-withDimReduction.txt"),
            quote = F, sep = "\t", col.names = NA)

#################################################################################

ds_genes <- scan('/groups/umcg-biogen/tmp01/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/data/multiple_sclerosis_patsopoulos_harm_jan_enrichtments_exHla.txt', character(), quote = "")
gene_trans_df <- read.csv("/groups/umcg-biogen/tmp01/annotation/gencode.v32.primary_assembly.annotation.collapsedGenes.ProbeAnnotation.TSS.txt.gz", sep = "\t")[, c("ArrayAddress", "Symbol")]
gene_trans_df <- gene_trans_df[complete.cases(gene_trans_df),]
gene_trans_df <- within(gene_trans_df, gene <- data.frame(do.call('rbind', strsplit(as.character(gene_trans_df$ArrayAddress),'.',fixed=TRUE))))
ds_genes_hgnc_smybols <- gene_trans_df[gene_trans_df$gene$X1 %in% ds_genes, "Symbol"]
ds_genes_matrix <- as.matrix(GetAssayData(object = seurat_object[ds_genes_hgnc_smybols, ], slot = "counts", assay = "SCT"))

ds_genes_max <- apply(ds_genes_matrix, 2, max)
names(ds_genes_max) <- colnames(x = seurat_object)

seurat_object <- AddMetaData(
  object = seurat_object,
  metadata = ds_genes_max,
  col.name = 'ds_max_expr'
)

for (dim_red in c("umap", "tsne")) {
  FeaturePlot(object = seurat_object, reduction = dim_red, features = "ds_max_expr")
  ggsave(paste0("ROSMAP_scRNAseq_", dim_red, "_groupedByds_max_expr.pdf"), width = 10, height = 10)
}

ds_genes_mean <- colMeans(ds_genes_matrix)
seurat_object <- AddMetaData(
  object = seurat_object,
  metadata = ds_genes_mean,
  col.name = 'ds_mean_expr'
)

for (dim_red in c("umap", "tsne")) {
  FeaturePlot(object = seurat_object, reduction = dim_red, features = "ds_mean_expr")
  ggsave(paste0("ROSMAP_scRNAseq_", dim_red, "_groupedByds_mean_expr.pdf"), width = 10, height = 10)
}

ds_top_gene_expr <- t(as.matrix(GetAssayData(object = seurat_object["UBALD1", ], slot = "counts", assay = "SCT")))
seurat_object <- AddMetaData(
  object = seurat_object,
  metadata = ds_top_gene_expr,
  col.name = 'ds_top_geneUBALD1'
)

for (dim_red in c("umap", "tsne")) {
  FeaturePlot(object = seurat_object, reduction = dim_red, features = "ds_top_geneUBALD1")
  ggsave(paste0("ROSMAP_scRNAseq_", dim_red, "_groupedByds_top_geneUBALD1.pdf"), width = 10, height = 10)
}

trans_genes <- c("UBE2L3", "PSD3", "MPP1", "KLHDC3", "FAM111A", "SPRYD3", "TPRG1L", "CLSTN1", "WIPI2", "THY1", "TMEM246", "HPCAL4", "B4GALNT1")
trans_genes_matrix <- as.matrix(GetAssayData(object = seurat_object[trans_genes, ], slot = "counts", assay = "SCT"))
trans_genes_mean <- colMeans(trans_genes_matrix)
seurat_object <- AddMetaData(
  object = seurat_object,
  metadata = trans_genes_mean,
  col.name = 'trans_mean_expr'
)

for (dim_red in c("umap", "tsne")) {
  FeaturePlot(object = seurat_object, reduction = dim_red, features = "trans_mean_expr")
  ggsave(paste0("ROSMAP_scRNAseq_", dim_red, "_groupedBytrans_mean_expr.pdf"), width = 10, height = 10)
}

cells_cell_type <- seurat_object[, seurat_object@active.ident == "Ex"]
for (dim_red in c("umap", "tsne")) {
  for (name in c("broad.cell.type", "Subcluster")) {
    DimPlot(cells_cell_type, reduction = dim_red, group.by = name)
    ggsave(paste0("ROSMAP_scRNAseq_EX_", dim_red, "_groupedBy", name ,".pdf"), width = 10, height = 10)
  }
  FeaturePlot(object = cells_cell_type, reduction = dim_red, features = "trans_mean_expr")
  ggsave(paste0("ROSMAP_scRNAseq_EX_", dim_red, "_groupedBytrans_mean_expr.pdf"), width = 10, height = 10)
}

#################################################################################

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

###############################################################################

link1_df <- data.frame(colnames(cells_cell_type),
                       sprintf("cell%d", 1:length(colnames(cell_type_matrix))))
colnames(link1_df) <- c("cellid", "cellMaskedID")

link2_df <- cells_cell_type@meta.data
link2_df["cellid"] <- rownames(link2_df)
rownames(link2_df) <- NULL
link2_df <- link2_df[, c("cellid", "projid")]

count <- 1
for (sample in unique(link2_df["projid"])) {
  link2_df[link2_df["projid"] == sample, "sampleMaskedID"] <- paste0("sample", count)
  count <- count + 1
}

link_df <- merge(link1_df, link2_df, by="cellid")
link_df <- mtcars[order(mpg),]
write.table(link_df[, c("maskedid")],
            file = paste0(work_dir, "cell_links.txt"),
            quote = F, sep = "\t", col.names = NA)

###############################################################################


cells_cell_type <- seurat_object[, seurat_object@active.ident == "Ex"]
cell_type_matrix <- as.matrix(GetAssayData(object = cells_cell_type, slot = "counts", assay = "SCT"))
rownames(cell_type_matrix) <- gene_trans_df[match(rownames(cell_type_matrix), gene_trans_df[,2]), 1]
cell_type_matrix <- cell_type_matrix[!is.na(rownames(cell_type_matrix)), ]
colnames(cell_type_matrix) <- sprintf("cell%d", 1:length(colnames(cell_type_matrix)))
write.table(cell_type_matrix, file = paste0(work_dir, "Ex_allCells_expression", ".tsv"), quote = F, sep = "\t", col.names = NA)