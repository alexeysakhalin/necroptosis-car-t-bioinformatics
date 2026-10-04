# Bioinformatics methods

## Bulk differential expression and gene selection

Downstream analysis used the supplied differential-expression summaries for the T-cell co-culture versus control and TNF+IFNγ versus control comparisons. The input fields were gene identifier, overall baseMean, log2 fold change, nominal P value and adjusted P value. The count-based model fitting was performed upstream. Reported effects and probabilities were retained without refitting or reconstructing a count matrix. The T-cell input had already undergone the upstream exclusion of selected T-cell marker genes.

Gene lists required overall baseMean ≥ 30, a finite log2 fold change, supplied adjusted P < 0.05 and absolute log2 fold change > 1. Missing or infinite effect estimates were retained in the full result tables and excluded from these lists and volcano plots. Thirty-one comma-containing identifiers in the T-cell input were reconciled with their dot-containing counterparts in the cytokine table using an explicit, collision-free mapping. Shared upregulated genes were the intersection of the two positive gene lists. The remaining genes were assigned to the corresponding comparison-specific set. The 30 shared genes with the largest mean log2 fold change across the two comparisons were plotted on an uncentered sequential scale.

## Over representation analysis

GO resource attribution follows the Gene Ontology Consortium (2023). Human symbols were mapped to unique Entrez gene identifiers using org.Hs.eg.db 3.18.0. All recorded symbol-to-identifier mappings were retained, with each Entrez identifier counted once. GO Biological Process, Molecular Function and Cellular Component annotations were obtained from GO.db 3.18.0 and the propagated GOALL mappings. Human KEGG pathway annotations were retrieved from the official REST service on 4 October 2026; retrieval hashes and release information accompany the results.

Each comparison-specific set was evaluated against the genes eligible for list selection in that comparison. The shared set was evaluated against the intersection of the two eligible backgrounds. These backgrounds contained 11,922, 11,668 and 11,481 gene symbols for the T-cell, TNF+IFNγ and shared analyses, respectively. For each annotation collection, the background and query were restricted to mapped genes present in that collection. Terms containing 10–500 background genes were tested using the upper tail of the hypergeometric distribution. Benjamini–Hochberg adjustment included all eligible terms, including terms with no query overlap, separately for each gene set and annotation collection. Up to ten terms with adjusted P < 0.05 were shown per collection, ordered by adjusted P, nominal P and term identifier. Gene ratio denotes overlapping query genes divided by the number of annotated query genes in that collection.

## Preranked enrichment

The Hallmark collection (Liberzon et al., 2015) was analyzed using the preranked GSEA framework (Subramanian et al., 2005; Korotkevich et al., 2021). The Hallmark apoptosis, TNF signaling via NF-κB and interferon-γ response gene sets were taken from MSigDB 2026.1.Hs. Genes with overall baseMean ≥ 10, an available effect direction and a nominal P value in (0, 1] were ranked by sign(log2 fold change) multiplied by −log10(P). Infinite supplied effects contributed only their sign when the ranking statistic was finite; missing effects were excluded. Equal scores were ordered alphabetically by gene identifier. The resulting rankings contained 14,128 and 13,777 genes for the T-cell and TNF+IFNγ comparisons, respectively.

Enrichment was calculated with fgsea 1.32.2 using fgseaSimple, 10,000 permutations, standard signed enrichment scoring, a weight exponent of 1, gene-set sizes of 10–500 and serial execution. The random seed was 1234. Benjamini–Hochberg correction was applied to the three selected pathways separately within each comparison. These adjusted probabilities do not represent a screen of the complete Hallmark collection. Running enrichment scores were independently recomputed from the ranked inputs.

## Targeted single cell expression

The CD3-positive worksheets of Supplementary Tables S1–S4 contained nonnegative integer molecular counts for 259 genes. The inputs comprised 4,051 cells cultured without IL-2, 1,146 cells cultured with IL-2 and two HeLa co-culture samples containing 2,240 and 3,695 cells. All 11,132 supplied CD3-positive cells were retained. No additional whole-transcriptome quality-control thresholds were inferred for this targeted panel. Sample identifiers were prepended to cell identifiers before combining matrices.

The count matrices were analyzed with Seurat (Hao et al., 2024), using Seurat 5.0.1 and SeuratObject 5.0.1 in R 4.3.3. Seurat's standard conversion of underscores to hyphens in feature names was recorded. Counts were normalized with LogNormalize and a scale factor of 10,000. All 259 genes were scaled, and 20 principal components were calculated. Components 1–15 were used for neighbor-graph construction, Louvain clustering at resolution 0.4 and UMAP visualization. PCA and UMAP seeds were 42; the clustering seed was 0. Eight clusters were annotated using their marker-expression patterns, with complete cluster means and detection fractions retained.

Per-cell module scores were the arithmetic mean of log-normalized expression for the specified signature genes, without a control-gene subtraction. The signatures covered TCR activation, FAS/FASL signaling, cytotoxicity, proliferation and naive/memory states. All signatures and state labels are specified in the analysis configuration. For the condition heatmap, per-cell scores were averaged within condition, pooling the two co-culture samples, and each module was standardized across the three condition means using the sample standard deviation. Cluster proportions used captured cells as the denominator. These figures describe the captured samples and do not use individual cells as independent biological replicates.

## DepMap expression

CCLE expression data were accessed through DepMap (Ghandi et al., 2019).

DepMap Public 25Q3 (released 25 September 2025) was identified by exact MD5 matches of both supplied files to the official download catalogue. The supplied DepMap protein-coding log2(TPM+1) matrix contained 1,754 expression profiles. Only profiles marked IsDefaultEntryForModel = Yes were selected, leaving 1,699 models. Models without lineage or cell-line labels were removed, along with the Lymphoid, Myeloid, Fibroblast, Embryonal and Other lineages and records labeled Non-Cancerous or Unknown in the primary-disease field. The analysis cohort comprised 1,343 models. Each gene was standardized across this cohort using its mean and sample standard deviation.

Figure 2E displays 41 genes from the original five-group heatmap panel across A431, H460, HeLa, HepG2, SW837 and T-47D. Colors saturate at z = −2 and z = 2; no measurements are replaced by missing-value tiles. S1C uses the original core apoptosis, pyroptosis and necroptosis signatures. S1E uses the broader sensing and execution signatures; the execution score is the equally weighted mean of the three execution-module scores. The different panel gene lists are recorded separately. These expression scores are descriptive and do not directly measure pathway activation or cell-death competence.

## Clinical cohort summaries

The supplied TIMER3 export (manuscript reference 24) was interpreted using its numeric fields: heatTable.p contained signed Cox Z statistics and heatTable.rho contained their two-sided normal-approximation P values. This interpretation was confirmed for all 234 available numeric gene/cohort records. Sixteen cohorts were present, including eight non-melanoma cohorts. The panels summarize the non-melanoma subset using the core gene signatures.

Nominal Fisher probabilities were calculated from available gene/cohort P values. Signed Stouffer scores were calculated as the sum of available Z statistics divided by the square root of their number. Leave-one-cohort-out values used the same calculation after omitting each non-melanoma cohort. Missing gene/cohort results were not imputed. Gene tests within a cohort may be correlated; the combined statistics are descriptive nominal summaries, without a dependence correction. They are not estimates of pathway strength or analyses of the six experimental cell lines.

## Compound sensitivity

The drug-response analysis used CTRPv2 v22 AUC and metadata tables (Basu et al., 2013; Seashore-Ludlow et al., 2015; Rees et al., 2016). Twenty-five DNMT/HDAC-targeting compounds were selected by the explicit source-name rules recorded in the configuration. HeLa was identified by index_ccl 123; the two NCI-H460 records, 359 and 362, were averaged for each compound. Source compound and cell-line identifiers are retained in the output tables. The source metadata label the HeLa record as SNP-unconfirmed.

Complete compound pairs were compared with a two-sided paired t test. The reported effect is the mean within-compound AUC difference, NCI-H460 minus HeLa, with a 95% t interval. The unit of pairing is the compound. Compounds can share mechanisms, so this comparison describes the selected panel and does not estimate variation across independently sampled patients or cell-line replicates. Lower CTRPv2 AUC denotes greater sensitivity.

## Software and verification

Python 3.12.14, NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0, Matplotlib 3.10.8 and openpyxl 3.1.5 were used for table processing, selected statistical calculations and plotting. R sessions and installed package versions accompany the results. Source checksums, independent calculations in R and Python, and consistency checks of the figure source tables were completed before rendering the revised panels.
