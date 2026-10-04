suppressPackageStartupMessages({
  library(fgsea)
  library(BiocParallel)
  library(jsonlite)
})

set.seed(1234)
output_dir <- "results/tables/enrichment"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
gene_sets <- jsonlite::fromJSON("data/annotations/hallmark_selected.json")
pathways <- lapply(gene_sets, `[[`, "geneSymbols")
pathways <- pathways[c("HALLMARK_APOPTOSIS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE")]
results <- list()
curves <- list()
for (condition in c("T_cells", "TNF_IFNg")) {
  table <- read.csv(file.path("results/tables/bulk", paste0(condition, "_gsea_ranks.csv")))
  stopifnot(!anyDuplicated(table$gene), all(is.finite(table$rank_score)))
  ranking <- setNames(table$rank_score, table$gene)
  stopifnot(all(diff(ranking) <= 0))
  result <- as.data.frame(fgsea::fgseaSimple(pathways = pathways, stats = ranking,
                                           nperm = 10000, minSize = 10, maxSize = 500,
                                           gseaParam = 1, scoreType = "std",
                                           BPPARAM = BiocParallel::SerialParam()))
  stopifnot(nrow(result) == 3L)
  result$condition <- condition
  result$permutations <- 10000L
  result$tested_pathways <- 3L
  result$leading_edge <- vapply(result$leadingEdge, paste, collapse = ";", character(1))
  result$leadingEdge <- NULL
  results[[condition]] <- result
  for (pathway in names(pathways)) {
    hit <- names(ranking) %in% pathways[[pathway]]
    weights <- abs(ranking)
    increments <- ifelse(hit, weights / sum(weights[hit]), -1 / sum(!hit))
    running <- cumsum(increments)
    independent_es <- if (max(running) > -min(running)) max(running) else min(running)
    package_es <- result$ES[result$pathway == pathway]
    stopifnot(abs(independent_es - package_es) < 1e-12)
    curves[[paste(condition, pathway)]] <- data.frame(condition = condition, pathway = pathway,
                                                     rank = seq_along(ranking), gene = names(ranking),
                                                     rank_score = as.numeric(ranking), member = hit,
                                                     running_enrichment_score = as.numeric(running))
  }
}
write.csv(do.call(rbind, results), file.path(output_dir, "hallmark_gsea_results.csv"), row.names = FALSE)
connection <- gzfile(file.path(output_dir, "hallmark_gsea_running_scores.csv.gz"), open = "wt")
write.csv(do.call(rbind, curves), connection, row.names = FALSE)
close(connection)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "gsea_R_session_info.txt"))
print(do.call(rbind, results)[, c("condition", "pathway", "size", "ES", "NES", "pval", "padj")])
