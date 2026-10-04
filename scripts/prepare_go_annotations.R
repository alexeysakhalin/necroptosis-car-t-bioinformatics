suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

output_dir <- "data/annotations"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
files <- file.path("results/tables/bulk", paste0(c("Common", "T_cells_only", "TNF_IFNg_only"), "_enrichment_universe.csv"))
symbols <- sort(unique(unlist(lapply(files, function(path) read.csv(path)$gene))))
available <- intersect(symbols, keys(org.Hs.eg.db, keytype = "SYMBOL"))
mapping <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = available, columns = "ENTREZID", keytype = "SYMBOL"))
mapping <- mapping[!is.na(mapping$ENTREZID), ]
write.csv(mapping, file.path(output_dir, "gene_symbol_to_entrez.csv"), row.names = FALSE)
write.csv(data.frame(gene = setdiff(symbols, mapping$SYMBOL)), file.path(output_dir, "unmapped_gene_symbols.csv"), row.names = FALSE)
go <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = unique(mapping$ENTREZID),
                                columns = c("GOALL", "ONTOLOGYALL"), keytype = "ENTREZID"))
go <- go[!is.na(go$GOALL), ]
terms <- AnnotationDbi::select(GO.db, keys = unique(go$GOALL), columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
stopifnot(!anyDuplicated(terms$GOID))
connection <- gzfile(file.path(output_dir, "go_memberships.csv.gz"), open = "wt")
write.csv(go, connection, row.names = FALSE)
close(connection)
write.csv(terms, file.path(output_dir, "go_terms.csv"), row.names = FALSE)
write.csv(AnnotationDbi::metadata(org.Hs.eg.db), file.path(output_dir, "org_hs_eg_db_metadata.csv"), row.names = FALSE)
write.csv(AnnotationDbi::metadata(GO.db), file.path(output_dir, "go_db_metadata.csv"), row.names = FALSE)
writeLines(c(paste("org.Hs.eg.db", packageVersion("org.Hs.eg.db")), paste("GO.db", packageVersion("GO.db"))),
           file.path(output_dir, "annotation_versions.txt"))
message("GO memberships exported: ", nrow(go), "; mapped symbols: ", length(unique(mapping$SYMBOL)))
