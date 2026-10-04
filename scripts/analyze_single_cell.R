suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(jsonlite)
  library(Matrix)
})

set.seed(123)
input_dir <- "data/source/single_cell"
output_dir <- "results/tables/single_cell"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
configuration <- jsonlite::fromJSON("config/analysis.json")
samples <- data.frame(
  sample_id = c("Tcells_noIL2", "Tcells_IL2", "HeLa_coculture_R1", "HeLa_coculture_R2"),
  condition = c("Tcells_noIL2", "Tcells_IL2", "HeLa_coculture", "HeLa_coculture"),
  filename = c("Supplementary_Table_S4_Tcells_noIL2.xlsx",
               "Supplementary_Table_S1_Tcells_IL2.xlsx",
               "Supplementary_Table_S2_Tcells_HeLa_coculture_rep1.xlsx",
               "Supplementary_Table_S3_Tcells_HeLa_coculture_rep2.xlsx"),
  expected_cells = c(4051L, 1146L, 2240L, 3695L)
)

read_sample <- function(i) {
  values <- as.data.frame(readxl::read_excel(file.path(input_dir, samples$filename[i]), sheet = "CD3+ cells"))
  stopifnot(nrow(values) == samples$expected_cells[i], ncol(values) == 260L)
  stopifnot(names(values)[1] == "Cell_Index", !anyNA(values), !anyDuplicated(values$Cell_Index))
  cell_id <- paste(samples$sample_id[i], values$Cell_Index, sep = "_")
  counts <- t(as.matrix(values[, -1, drop = FALSE]))
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)))
  stopifnot(all(colSums(counts[c("CD3D", "CD3E", "CD3G"), , drop = FALSE]) > 0))
  colnames(counts) <- cell_id
  metadata <- data.frame(sample_id = samples$sample_id[i], condition = samples$condition[i],
                         Cell_Index = values$Cell_Index, row.names = cell_id)
  list(counts = counts, metadata = metadata)
}

input <- lapply(seq_len(nrow(samples)), read_sample)
stopifnot(all(vapply(input, function(x) identical(rownames(x$counts), rownames(input[[1]]$counts)), logical(1))))
counts <- do.call(cbind, lapply(input, `[[`, "counts"))
metadata <- do.call(rbind, lapply(input, `[[`, "metadata"))
stopifnot(ncol(counts) == 11132L, !anyDuplicated(colnames(counts)))
seu <- CreateSeuratObject(counts = counts, meta.data = metadata, project = "HeLa_Tcell_coculture")
write.csv(data.frame(source_gene = rownames(counts), analysis_gene = rownames(seu)),
          file.path(output_dir, "gene_name_mapping.csv"), row.names = FALSE)
seu$condition <- factor(seu$condition, levels = c("Tcells_noIL2", "Tcells_IL2", "HeLa_coculture"))
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
features <- rownames(seu)
seu <- ScaleData(seu, features = features, verbose = FALSE)
seu <- RunPCA(seu, features = features, npcs = 20, seed.use = 42, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:15, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.4, algorithm = 1, random.seed = 0, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:15, seed.use = 42, verbose = FALSE)

normalized <- GetAssayData(seu, assay = "RNA", layer = "data")
for (module_name in names(configuration$single_cell_modules)) {
  genes <- configuration$single_cell_modules[[module_name]]
  stopifnot(all(genes %in% rownames(normalized)))
  seu[[module_name]] <- Matrix::colMeans(normalized[genes, , drop = FALSE])
}

coordinates <- as.data.frame(Embeddings(seu, "umap"))
names(coordinates) <- c("UMAP_1", "UMAP_2")
cell_table <- cbind(data.frame(cell_id = rownames(seu@meta.data)), seu@meta.data, coordinates)
write.csv(cell_table, file.path(output_dir, "cell_metadata_and_coordinates.csv"), row.names = FALSE)
write.csv(samples, file.path(output_dir, "sample_manifest.csv"), row.names = FALSE)
cluster_counts <- as.data.frame(table(sample_id = seu$sample_id, cluster = seu$seurat_clusters))
names(cluster_counts)[3] <- "cells"
write.csv(cluster_counts, file.path(output_dir, "cluster_counts_by_sample.csv"), row.names = FALSE)

clusters <- sort(unique(as.character(seu$seurat_clusters)))
mean_expression <- sapply(clusters, function(cluster) Matrix::rowMeans(normalized[, seu$seurat_clusters == cluster, drop = FALSE]))
fraction_detected <- sapply(clusters, function(cluster) Matrix::rowMeans(counts[, seu$seurat_clusters == cluster, drop = FALSE] > 0))
colnames(mean_expression) <- colnames(fraction_detected) <- paste0("cluster_", clusters)
write.csv(data.frame(gene = rownames(mean_expression), mean_expression), file.path(output_dir, "cluster_mean_log_expression.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(fraction_detected), fraction_detected), file.path(output_dir, "cluster_fraction_detected.csv"), row.names = FALSE)
object_path <- file.path(output_dir, "CD3_Tcell_Seurat_object.rds")
temporary_object_path <- paste0(object_path, ".tmp")
saveRDS(seu, temporary_object_path)
restored <- readRDS(temporary_object_path)
stopifnot(validObject(restored),
          identical(seu@meta.data, restored@meta.data),
          identical(GetAssayData(seu, layer = "counts"), GetAssayData(restored, layer = "counts")),
          identical(GetAssayData(seu, layer = "data"), GetAssayData(restored, layer = "data")),
          identical(Embeddings(seu, "umap"), Embeddings(restored, "umap")))
stopifnot(file.rename(temporary_object_path, object_path))
rm(restored)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_session_info.txt"))
jsonlite::write_json(list(cells = ncol(seu), genes = nrow(seu), clusters = length(clusters),
                         sample_cells = as.list(table(seu$sample_id)),
                         normalization = "LogNormalize", scale_factor = 10000,
                         pca_components = 20, neighbor_dimensions = 1:15, resolution = 0.4,
                         clustering_algorithm = "Louvain", clustering_seed = 0, pca_seed = 42, umap_seed = 42),
                    file.path(output_dir, "analysis_summary.json"), pretty = TRUE, auto_unbox = TRUE)
message("Single-cell analysis completed: ", ncol(seu), " cells; ", length(clusters), " clusters.")
