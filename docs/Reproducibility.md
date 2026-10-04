# Reproducing the manuscript bioinformatics

The release concerns **Loss of Necroptosis Drives Tumor Resistance to CAR-T-mediated Immunologic Cell Death**, by Levchuk, Petukhov, Tursymbek, Ilyassova and Barlev. It contains 19 analyzed panels in 18 PDF and 18 PNG files, with Figures 3D and 3E combined. Figure 6 is outside the analysis scope.

## Reference computation

The downstream pipeline was rerun from the recorded source files in a separate working directory on 4 October 2026. Bulk selection, public-data analysis, targeted single-cell analysis, GO annotation preparation, over-representation analysis, preranked enrichment, independent R calculations and final validation completed successfully. All numerical result tables agreed with the saved reference results to a tolerance of 10⁻¹². Recreated serialized Seurat objects can differ in internal session metadata while retaining identical reported tables and coordinates.

The reference validation contains 77 source, numerical and panel-consistency checks, plus 27 independent R checks. The counts describe this release and are not a certification of experimental measurements outside its scope. Full check details and file hashes are included in `results/tables/validation/`.

## Software

The reference environment uses Python 3.12.14 and R 4.3.3. Install the Python packages pinned in `requirements.txt`. `environment/R_packages.tsv` records the complete R package inventory, and the single-cell and enrichment results include session information. Required direct R dependencies are checked by `scripts/check_environment.R`. In particular, the reference uses fgsea 1.32.2 with R 4.3.3; installing a current Bioconductor release without pinning versions does not recreate this environment. Archive package sources or an equivalent environment with the recorded versions are required. A fresh operating-system installation has not been independently validated by this release.

Floating-point results and stochastic embeddings may vary with different numerical libraries, package versions or hardware. Seeds, thread settings and all numerical parameters used here are recorded in the scripts and configuration. Compare numeric tables with an explicit tolerance rather than requiring identical bytes for serialized R objects.

## Full downstream rerun

Unpack the analysis bundle and obtain the two DepMap Public 25Q3 inputs under their applicable terms. Place them under `data/source/depmap/`. The source manifest identifies all 16 input files by path, byte size and SHA-256. Before calculation, run `python scripts/check_inputs.py`; a truncated, altered or substituted source file stops the pipeline.

Use the exact KEGG resources recorded in `config/kegg_reference_manifest.json`. They were retrieved on 4 October 2026 and are not redistributed as a database. The enrichment script can download current resources from the official service, but it rejects a checksum mismatch by default. A different resource snapshot is a new analysis, not a reproduction of this release. Other gene annotations and the three selected Hallmark gene sets are frozen in `data/annotations/`.

Run `python scripts/run_pipeline.py --kegg-cache /path/to/reference/kegg`. The command checks the inputs and required R versions, recalculates results, validates them and renders the figures. Add `--regenerate-go` to regenerate the distributed GO annotations from org.Hs.eg.db and GO.db 3.18.0. It does not refit the upstream bulk RNA-seq model; the supplied differential-expression summaries are the starting point.

## Rebuilding figures from saved results

Run `python scripts/verify_cached_results.py` followed by `python scripts/plot_figures.py`. These commands use the distributed results and do not require the restricted full DepMap matrix or KEGG database. The verifier checks the protected files and independently recalculates every over-representation probability and multiple-testing adjustment from the exported contingency counts. It does not re-establish the original KEGG annotation membership.

The plotting script refuses result files that changed after numerical validation. All panel values, cell coordinates, module definitions, full enrichment tests and displayed-term selections are available for inspection.

## Boundary of the release

S1F cannot be reconstructed from the available files because the paired six-cell-line source matrix is absent. Raw qPCR values and the transformation underlying S2C–D were not supplied, and the meaning of blank tiles in S2G is unresolved. The sequence-generation description in the submitted Methods also requires reconciliation with the accession metadata. These issues are recorded in the author review document and are not resolved by successful execution of the downstream pipeline.
