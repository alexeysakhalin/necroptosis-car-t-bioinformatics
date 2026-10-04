suppressPackageStartupMessages(library(readxl))
out <- "results/tables/validation"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
record <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = passed, detail = detail)
  if (!isTRUE(passed)) stop("Verification failed: ", name)
}

corrections <- read.csv("config/gene_identifier_corrections.csv", check.names = FALSE)
mapping <- setNames(corrections$corresponding_cytokine_id, corrections$supplied_Tcell_id)
up <- list()
for (condition in c("T_cells", "TNF_IFNg")) {
  source <- if (condition == "T_cells") {
    read.csv("data/source/figure4/RNAseq_filtered_without_tcell_genes.csv", sep = ";", check.names = FALSE)
  } else {
    as.data.frame(read_excel("data/source/figure4/Cytokines.xlsx", sheet = "WT vs TNF-IFN",
                            col_types = c("text", "numeric", "text", "numeric", "numeric"), na = c("", "NA")))
  }
  effect <- suppressWarnings(as.numeric(gsub("\u2013", "-", source$log2FoldChange, fixed = TRUE)))
  gene <- trimws(source$id)
  matched <- gene %in% names(mapping)
  gene[matched] <- mapping[gene[matched]]
  chosen <- !is.na(source$padj) & is.finite(effect) & !is.na(source$baseMean) & source$baseMean >= 30 & source$padj < 0.05 & effect > 1
  up[[condition]] <- sort(gene[chosen])
  python <- read.csv(file.path("results/tables/bulk", paste0(condition, "_statistics.csv")), check.names = FALSE)
  record(paste(condition, "upregulated gene identities"), identical(up[[condition]], sort(python$gene[python$upregulated == "True"])), length(up[[condition]]))
  ranking <- read.csv(file.path("results/tables/bulk", paste0(condition, "_gsea_ranks.csv")), check.names = FALSE)
  rankable <- !is.na(effect) & !is.na(source$pval) & source$pval > 0 & source$pval <= 1 & !is.na(source$baseMean) & source$baseMean >= 10
  expected <- setNames(sign(effect[rankable]) * -log10(source$pval[rankable]), gene[rankable])
  record(paste(condition, "ranking membership"), setequal(names(expected), ranking$gene), nrow(ranking))
  difference <- max(abs(expected[ranking$gene] - ranking$rank_score))
  record(paste(condition, "ranking scores"), difference < 1e-12, difference)
}
record("Shared upregulated genes", length(intersect(up$T_cells, up$TNF_IFNg)) == 93L, "93")

expression <- read.csv("results/tables/public_data/depmap_selected_expression.csv", check.names = FALSE)
metadata <- read.csv("data/source/depmap/Model.csv", check.names = FALSE)
expression <- expression[expression$IsDefaultEntryForModel == "Yes", ]
joined <- merge(expression, metadata, by = "ModelID", sort = FALSE)
keep <- !is.na(joined$OncotreeLineage) & !is.na(joined$StrippedCellLineName)
keep <- keep & !joined$OncotreeLineage %in% c("Lymphoid", "Myeloid", "Fibroblast", "Embryonal", "Other")
keep <- keep & !grepl("Non-Cancerous|Unknown", joined$OncotreePrimaryDisease, ignore.case = TRUE)
cohort <- joined[keep, ]
gene_columns <- grepl(" \\([0-9]+\\)$", names(cohort))
values <- cohort[, gene_columns]
colnames(values) <- sub(" \\([0-9]+\\)$", "", colnames(values))
scaled <- scale(as.matrix(values), center = TRUE, scale = TRUE)
rownames(scaled) <- cohort$ModelID
python <- read.csv("results/tables/public_data/depmap_gene_z_scores.csv", check.names = FALSE)
record("DepMap analysis cohort", setequal(cohort$ModelID, python$ModelID) && nrow(cohort) == 1343L, nrow(cohort))
difference <- max(abs(scaled[python$ModelID, colnames(values)] - as.matrix(python[, colnames(values)])))
record("DepMap gene-wise normalization", difference < 1e-12, difference)

timer <- read.csv("data/source/figure2/Immunotherapy_outcome_PANoptosis.csv", check.names = FALSE)
z <- suppressWarnings(as.numeric(timer$heatTable.p))
p <- suppressWarnings(as.numeric(timer$heatTable.rho))
valid <- is.finite(z) & is.finite(p) & p > 0 & p <= 1
timer <- data.frame(cohort = timer$heatTable.cancer[valid], gene = timer$heatTable.variable[valid], z = z[valid], p = p[valid])
record("TIMER field interpretation", max(abs(timer$p - 2 * pnorm(-abs(timer$z)))) < 1e-12, nrow(timer))
groups <- list(Apoptosis = c("CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "APAF1"),
               Pyroptosis = c("GSDMD", "CASP1", "CASP4", "CASP5"), Necroptosis = c("RIPK3", "MLKL"))
python <- read.csv("results/tables/public_data/timer_pathway_summaries.csv")
for (name in names(groups)) {
  chosen <- timer$gene %in% groups[[name]] & !grepl("skcm|melano", timer$cohort, ignore.case = TRUE)
  probability <- pchisq(-2 * sum(log(timer$p[chosen])), df = 2 * sum(chosen), lower.tail = FALSE)
  expected <- python$fisher_p_nominal[python$cohort_group == "Non-melanoma" & python$pathway == name]
  record(paste("TIMER Fisher probability", name), length(expected) == 1L && abs(probability - expected) < 1e-12, probability)
}
pairs <- read.csv("results/tables/public_data/figure_s1g_paired_auc.csv", check.names = FALSE)
paired <- t.test(pairs$HeLa, pairs$`NCI-H460`, paired = TRUE)
record("CTRP paired probability", abs(paired$p.value - 0.02553484633540972) < 1e-12, paired$p.value)
ora <- read.csv("results/tables/enrichment/over_representation_all_tests.csv", check.names = FALSE)
expected <- phyper(ora$count - 1, ora$background_term_genes,
                   ora$background_annotated_genes - ora$background_term_genes,
                   ora$query_annotated_genes, lower.tail = FALSE)
difference <- max(abs(expected - ora$pvalue))
record("Over-representation probabilities", difference < 1e-12, difference)
for (group in unique(ora$group)) {
  for (collection in unique(ora$collection)) {
    chosen <- ora$group == group & ora$collection == collection
    corrected <- p.adjust(ora$pvalue[chosen], method = "BH")
    record(paste("Over-representation FDR", group, collection),
           max(abs(corrected - ora$padj[chosen])) < 1e-12, sum(chosen))
  }
}
write.csv(do.call(rbind, checks), file.path(out, "independent_R_checks.csv"), row.names = FALSE)
message("Independent R checks passed: ", length(checks))
