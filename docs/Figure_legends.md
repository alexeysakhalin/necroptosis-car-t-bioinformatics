# Bioinformatics figure legends

**Figure 2D.** Nominal Fisher combination of available gene-level clinical association P values for core apoptosis, pyroptosis and necroptosis signatures in the non-melanoma TIMER cohorts. Bars show −log10(P). There were 44, 21 and 8 gene/cohort tests, respectively, drawn from 8, 8 and 5 cohorts with data. Within-cohort gene dependence is not corrected; these are descriptive clinical summaries.

**Figure 2E.** Baseline expression of 41 cell-death-related genes in six experimental cell-line models. Gene-wise z scores were calculated across 1,343 selected default DepMap model profiles. Colors saturate at ±2; saturated cells retain their measured values in the source table. Gene-group labels identify sensing, apoptosis, necroptosis, pyroptosis/inflammasome and survival-associated genes.

**Figure 3D and E.** UMAP of 11,132 captured CD3-positive cells from the 259-gene targeted panel, colored by culture condition (D) or eight annotated clusters (E). Coordinates were calculated jointly and are identical in both panels. The two HeLa co-culture samples are pooled for condition coloring.

**Figure 3F.** Distribution of the eight annotated states among captured CD3-positive cells in each condition. Each bar sums to 100%; the denominator is 4,051 cells without IL-2, 1,146 cells with IL-2 or 5,935 pooled HeLa co-culture cells. Colors match E. Proportions describe captured cells and are not donor-level estimates.

**Figure 3G.** Condition means of five T-cell expression modules, standardized separately for each module across the three conditions using the sample standard deviation. Scores are means of log-normalized signature-gene expression per cell. Values in the tiles are z scores; the co-culture mean is weighted by captured cell numbers.

**Figure 4A.** Volcano plots from the full supplied differential-expression summaries. Genes with overall baseMean ≥ 30, finite log2 fold change and a positive adjusted P value are shown. Red and blue indicate supplied adjusted P < 0.05 and log2 fold change > 1 or < −1, respectively; gray denotes other eligible genes. Dashed lines mark these thresholds. Selected death-response genes are labeled whether or not they meet the significance criteria.

**Figure 4B.** Area-proportional overlap of upregulated genes selected with overall baseMean ≥ 30, finite log2 fold change > 1 and supplied adjusted P < 0.05. The sets contain 93 shared genes, 501 genes meeting the criteria only in the TNF+IFNγ comparison and 134 only in the T-cell comparison. Set-specific membership does not imply a significant difference between treatments.

**Figure 4C.** The 30 shared upregulated genes with the largest mean log2 fold change across the two comparisons. Colors show uncentered log2 fold changes on one sequential scale. Rows are ordered by the mean effect, with alphabetical ordering for ties.

**Figure 4D.** Preranked enrichment for three selected MSigDB Hallmark gene sets. Rankings use sign(log2 fold change) × −log10(nominal P) for genes with overall baseMean ≥ 10 and a finite ranking score. Curves show running enrichment scores; ticks mark gene-set members. NES denotes normalized enrichment score. Adjusted P values use Benjamini–Hochberg correction across these three pathways separately for each comparison, with 10,000 fgseaSimple permutations.

**Figure S1C.** Core death-pathway expression scores across all 1,343 selected default DepMap model profiles. Each score is the mean of the relevant gene-wise z scores. The six experimental models are highlighted. All profiles are shown using common axis limits across the three comparisons.

**Figure S1D.** Signed Stouffer summaries of the non-melanoma TIMER records. Red diamonds use all available records for each pathway. Black circles omit one of the eight non-melanoma cohorts at a time; omission of a cohort without data for a pathway leaves its score unchanged. These nominal summaries do not correct dependence among gene tests.

**Figure S1E.** TNF/IFNγ sensing versus death-execution expression scores in 1,343 DepMap models. Sensing is the mean gene-wise z score of the specified sensing signature. Execution is the equally weighted mean of the broader apoptosis, necroptosis and pyroptosis/inflammasome module scores. Highlighted models comprise the experimental panel; expression alone does not establish functional competence.

**Figure S1G.** CTRPv2 response-curve AUC for 25 selected DNMT/HDAC-targeting compounds in HeLa and NCI-H460. Each line joins the same compound. The two source NCI-H460 records are averaged per compound. Violin outlines summarize distributions; horizontal marks show medians. The annotated P value is from a two-sided paired t test across compounds. Lower AUC denotes greater sensitivity.

**Figure S2I.** Per-cell distributions of five expression-module scores in the three culture conditions. Boxes indicate medians and interquartile ranges; whiskers extend to the most extreme observations within 1.5 interquartile ranges; points outside are shown individually. Scores are mean log-normalized signature-gene expression. Sample sizes refer to captured cells, not independent donors.

**Figure S2J.** The joint UMAP split by culture condition, with identical coordinates, axis limits and state colors. The two HeLa co-culture samples are pooled. State labels and marker summaries are provided in the analysis tables.

**Figure S3A to C.** GO and KEGG over-representation of shared upregulated genes (A), genes meeting upregulation criteria only in the T-cell comparison (B) and only in the TNF+IFNγ comparison (C). Each set uses its corresponding eligible tested-gene background, restricted to annotated Entrez identifiers. The upper-tail hypergeometric P values were adjusted by Benjamini–Hochberg across all eligible terms in each collection and gene set, including zero-overlap terms. Up to ten terms with adjusted P < 0.05 are displayed per collection. Position shows the fraction of annotated query genes in a term, bubble size the number of overlapping genes and color −log10(adjusted P). Terms were tested only when they contained 10–500 background genes. Complete tested results and background counts accompany the figures.
