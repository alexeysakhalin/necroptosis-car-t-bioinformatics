# Necroptosis and CAR T cell resistance

Analysis code accompanying **Tumor cell death competence shapes responses to CAR T cells and inflammatory cytokines**, by Ksenia Levchuk, Alexey Petukhov, Shynggys Tursymbek, Bayansulu Ilyassova and Nikolai A. Barlev. The archived version 1.0.0 retains its deposited title, **Loss of Necroptosis Drives Tumor Resistance to CAR-T-mediated Immunologic Cell Death: Analysis Code and Figure Data**.

The analysis covers Figures 2D–E, 3D–G, 4A–D, S1C–E, S1G, S2I–J and S3A–C. Figure 6 is a separate analysis: its [software repository](https://github.com/alexeysakhalin/tumor-model-comparison-patient-3d-2d) accompanies the [version 1.1.1 source-data archive](https://doi.org/10.5281/zenodo.22863363). That Zenodo archive contains derived data, methods, provenance, validation records and figure exports; the companion scripts are distributed through GitHub. Neither archive contains the separate spatial simulations in Supplementary Figures S5 and S6. The six-line correlation matrix underlying S1F and the primary experimental measurements underlying the remaining panels are not inputs to this analysis.

The archived analysis and separate bulk and targeted single-cell source-table bundles are deposited under DOI [10.5281/zenodo.23136820](https://doi.org/10.5281/zenodo.23136820). The software repository is [alexeysakhalin/necroptosis-car-t-bioinformatics](https://github.com/alexeysakhalin/necroptosis-car-t-bioinformatics).

## Contents

| Directory | Contents |
| --- | --- |
| `scripts/` | Analysis, independent numerical verification and plotting |
| `config/` | Gene sets, thresholds, identifier corrections and source checksums |
| `data/source/` | Supplied differential-expression, targeted single-cell and public-resource tables |
| `data/annotations/` | Frozen GO annotations and three selected MSigDB Hallmark gene sets |
| `results/tables/` | Complete analysis results, panel values and validation records |
| `figures/` | Vector PDF and 600 dpi PNG panels |
| `docs/` | Methods, figure legends, results and bibliography |
| `environment/` | Software versions used for the reported results |

The complete analysis bundle and source data are distributed separately from the GitHub code repository. Third-party resources retain their own terms; see `docs/Data_sources.md`. Downloaded input files are identified by SHA-256 in `config/source_data_manifest.json`.

Download `Bioinformatics_Publication_Package.zip` from the repository release and unpack it into a new directory. It contains the code, annotations, saved results, figures and shareable source inputs in one project folder. The GitHub file tree contains the code and documentation; the release archive restores the complete project. The two full DepMap inputs and the full KEGG database must be obtained separately under their applicable terms.

## Reproduce the analysis

Run commands from the project root using Python 3.12 and R 4.3.3. The Python dependencies are pinned in `requirements.txt`; the recorded R package versions are in `environment/R_packages.tsv`. Required R packages include Seurat 5.0.1, SeuratObject 5.0.1, Matrix 1.6-5, readxl 1.4.3, jsonlite 1.8.8, AnnotationDbi 1.64.1, org.Hs.eg.db 3.18.0, GO.db 3.18.0, fgsea 1.32.2 and BiocParallel 1.36.0.

Place the input files in the paths recorded in the source manifest. The DepMap Public 25Q3 matrix and model annotations must match the recorded checksums; see `config/depmap_release.json`. The pipeline reads the supplied bulk differential-expression tables; it does not fit a new count-based model.

```bash
python -m pip install -r requirements.txt
python scripts/check_inputs.py
python scripts/run_pipeline.py --kegg-cache /path/to/reference/kegg
```

The driver checks source checksums and the required R package versions before analysis. It uses the distributed GO annotations; `--regenerate-go` recreates them with the recorded packages. Use `--rscript /path/to/Rscript` if R is not on the executable search path. The reference run and its scope are described in `docs/Reproducibility.md`.

To verify the saved result files and regenerate panels without downloading the full DepMap matrix or KEGG database:

```bash
python scripts/verify_cached_results.py
python scripts/plot_figures.py
```

`analyze_enrichment.py` obtains KEGG annotations from the official REST service if they are absent from its local cache. The published results use the retrieval recorded in `results/tables/enrichment/kegg_resource_manifest.json`. An exact repeat requires the same annotation snapshot. A newer KEGG download can change enrichment results. The script rejects changed resources by default; `--allow-new-kegg-snapshot` explicitly starts a new enrichment analysis. The complete KEGG database is not redistributed. Existing figure data can be inspected without downloading KEGG.

The plotting step verifies the numerical validation record and the hashes of the inputs, scripts and result tables. Numerical verification comprises source-integrity checks, independent R/Python calculations and consistency checks between panel data and analysis definitions. PDF panels retain vector text and axes; dense scatter points are rasterized within the PDF.

## Interpretation

Using overall baseMean ≥ 30, finite log2 fold change > 1 and supplied adjusted P < 0.05 gives 93 shared upregulated genes, 134 genes specific to the T-cell comparison and 501 genes specific to the TNF+IFNγ comparison. These are threshold-defined sets; membership in only one set does not establish a significant difference between the two treatment effects.

The single-cell analysis contains 11,132 captured CD3-positive cells across four input samples and 259 targeted genes. The two co-culture samples are pooled for the condition-level figures. Cell-level distributions are descriptive and are not treated as independent donor replicates.

TIMER gene/cohort results are summarized using nominal Fisher probabilities and signed Stouffer scores. Genes within a cohort can be dependent, so the combined values are descriptive summaries rather than dependence-adjusted pathway significance tests. They describe clinical cohorts and are separate from the six experimental cell lines.

CTRPv2 supplies the drug-response measurements in S1G. The 25 paired observations are compounds, not biological replicates. Lower AUC denotes greater compound sensitivity within this dataset.

## Citation and attribution

Use `CITATION.cff` to cite the archived analysis package. `docs/References_verified.bib` records the package bibliography; its reference identifiers are independent of the revised manuscript numbering. The DOI-based correspondence is given in `docs/Reference_number_mapping.csv`. Dataset and method references are listed in `docs/Additional_method_references.txt` and `docs/Data_sources.md`. The original data analysis and visualization contribution of Bayansulu Ilyassova is retained in the manuscript attribution.
