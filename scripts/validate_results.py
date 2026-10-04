"""Verify source integrity, independent numerical checks and panel source tables."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/tables/validation"
OUT.mkdir(parents=True, exist_ok=True)
config = json.loads((ROOT / "config/analysis.json").read_text())
checks = []

def check(stage, name, passed, detail):
    checks.append({"stage": stage, "check": name, "passed": bool(passed), "detail": detail})
    if not passed:
        raise ValueError(name + ": " + str(detail))

manifest = json.loads((ROOT / "config/source_data_manifest.json").read_text())
for row in manifest:
    path = ROOT / row["path"]
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024*1024):
            digest.update(block)
    check(1, "Source integrity: " + path.name, path.stat().st_size == row["bytes"] and digest.hexdigest() == row["sha256"], row["sha256"])

r_checks = pd.read_csv(OUT / "independent_R_checks.csv")
check(2, "Independent R calculations", r_checks.passed.all() and len(r_checks) >= 27, int(len(r_checks)))
panel_genes = sum(config["figure_2e_gene_groups"].values(), [])
panel_data = pd.read_csv(ROOT / "results/tables/public_data/figure_2e_expression_z_scores.csv")
check(3, "Figure 2E gene panel", len(panel_genes) == len(set(panel_genes)) == 41 and
      set(panel_genes) <= set(panel_data.columns) and len(panel_data) == 6,
      "41 source-panel genes; six models; separate gene lists for supplementary module scores")
sc_dir = ROOT / "results/tables/single_cell"
metadata = pd.read_csv(sc_dir / "cell_metadata_and_coordinates.csv")
samples = pd.read_csv(sc_dir / "sample_manifest.csv")
modules = config["single_cell_modules"]
normalized_parts = []
for row in samples.itertuples(index=False):
    counts = pd.read_excel(ROOT / "data/source/single_cell" / row.filename, sheet_name="CD3+ cells", index_col=0)
    check(1, "Cell count and count scale: " + row.sample_id,
          len(counts) == row.expected_cells and counts.shape[1] == 259 and np.isfinite(counts.to_numpy()).all() and
          (counts.to_numpy() >= 0).all() and (counts.to_numpy() == np.floor(counts.to_numpy())).all(), int(len(counts)))
    normalized = np.log1p(counts.div(counts.sum(axis=1), axis=0) * 10000)
    cell_ids = [row.sample_id + "_" + str(int(index)) for index in counts.index]
    normalized.index = cell_ids
    reference = metadata.set_index("cell_id").loc[cell_ids]
    for module, genes in modules.items():
        score = normalized[genes].mean(axis=1)
        error = float(np.max(np.abs(score.to_numpy() - reference[module].to_numpy())))
        check(2, "Single-cell module: " + row.sample_id + "/" + module, error < 1e-12, error)
    normalized.columns = normalized.columns.str.replace("_", "-", regex=False)
    normalized["cluster"] = reference.seurat_clusters.to_numpy()
    normalized_parts.append(normalized)
normalized = pd.concat(normalized_parts)
means = normalized.groupby("cluster").mean().T
means.columns = ["cluster_" + str(cluster) for cluster in means.columns]
reference = pd.read_csv(sc_dir / "cluster_mean_log_expression.csv", index_col=0)
error = float(np.max(np.abs(means.loc[reference.index, reference.columns].to_numpy() - reference.to_numpy())))
check(2, "Cluster marker means", error < 1e-12, error)
check(3, "State labels cover every cluster", set(map(str, metadata.seurat_clusters.unique())) == set(config["cell_states"]), int(metadata.seurat_clusters.nunique()))
marker_maxima = {"CCR7": "cluster_5", "TCF7": "cluster_5", "IL7R": "cluster_0", "GNLY": "cluster_3",
                 "GZMB": "cluster_3", "CTLA4": "cluster_2", "MKI67": "cluster_6", "TOP2A": "cluster_6",
                 "CD69": "cluster_7", "FASLG": "cluster_7"}
for gene, expected in marker_maxima.items():
    check(3, "Annotation marker: " + gene, reference.loc[gene].idxmax() == expected, expected)
check(3, "UMAP coordinates", np.isfinite(metadata[["UMAP_1", "UMAP_2"]].to_numpy()).all(), len(metadata))

bulk_dir = ROOT / "results/tables/bulk"
membership = pd.read_csv(bulk_dir / "upregulated_gene_membership.csv")
top = pd.read_csv(bulk_dir / "figure_4c_top_30_shared_genes.csv")
check(3, "Shared-gene heatmap selection", len(top) == 30 and set(top.gene) <= set(membership.loc[membership.group.eq("Common"), "gene"]) and
      (top[["T_cells_log2FoldChange", "TNF_IFNg_log2FoldChange"]] > 1).all().all(), len(top))
for group in membership.group.unique():
    background = set(pd.read_csv(bulk_dir / (group + "_enrichment_universe.csv")).gene)
    check(3, "Enrichment input and background: " + group, set(membership.loc[membership.group.eq(group), "gene"]) <= background, len(background))

enrichment_dir = ROOT / "results/tables/enrichment"
ora = pd.read_csv(enrichment_dir / "over_representation_all_tests.csv")
check(3, "Over-representation contingency tables", ((ora["count"] <= ora.query_annotated_genes) &
       (ora["count"] <= ora.background_term_genes) & (ora.background_term_genes <= ora.background_annotated_genes) &
       ora.pvalue.between(0, 1) & ora.padj.between(0, 1)).all(), len(ora))
shown = pd.read_csv(enrichment_dir / "figure_s3_displayed_terms.csv")
check(3, "Displayed enrichment terms", shown.padj.lt(.05).all() and shown.groupby(["group", "collection"]).size().le(10).all(), len(shown))

gsea = pd.read_csv(enrichment_dir / "hallmark_gsea_results.csv")
curves = pd.read_csv(enrichment_dir / "hallmark_gsea_running_scores.csv.gz")
pathways = json.loads((ROOT / "data/annotations/hallmark_selected.json").read_text())
for condition, frame in gsea.groupby("condition"):
    check(3, "GSEA multiplicity: " + condition,
          len(frame) == 3 and np.allclose(stats.false_discovery_control(frame.pval.to_numpy()), frame.padj, atol=1e-12, rtol=0), 3)
    ranks = pd.read_csv(bulk_dir / (condition + "_gsea_ranks.csv"))
    check(3, "GSEA rank order: " + condition, ranks.rank_score.diff().dropna().le(0).all() and not ranks.gene.duplicated().any(), len(ranks))
    for row in frame.itertuples(index=False):
        hits = ranks.gene.isin(pathways[row.pathway]["geneSymbols"]).to_numpy()
        weights = abs(ranks.rank_score.to_numpy())
        increments = np.where(hits, weights/weights[hits].sum(), -1/(~hits).sum())
        running = np.cumsum(increments)
        es = running.max() if running.max() > -running.min() else running.min()
        saved = curves[curves.condition.eq(condition) & curves.pathway.eq(row.pathway)]
        error = max(float(abs(es-row.ES)), float(np.max(abs(running-saved.running_enrichment_score.to_numpy()))))
        check(2, "GSEA running statistic: " + condition + "/" + row.pathway, error < 1e-11 and hits.sum() == row.size, error)
        check(3, "GSEA terminal running score: " + condition + "/" + row.pathway, abs(running[-1]) < 1e-11, float(running[-1]))

pd.DataFrame(checks).to_csv(OUT / "validation_checks.csv", index=False)
protected = [path for path in (ROOT / "config").glob("*") if path.is_file()]
protected += [path for path in (ROOT / "scripts").glob("*") if path.is_file() and not path.name.startswith("plot_")]
protected += [path for path in (ROOT / "data/annotations").rglob("*") if path.is_file()]
protected += [path for path in (ROOT / "results/tables").rglob("*") if path.is_file() and "validation" not in path.parts]
protected = [path for path in protected if not path.name.endswith((".tmp", ".pyc")) and not any(part.startswith(".") or part == "__pycache__" for part in path.relative_to(ROOT).parts)]
hashes = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(protected)}
report = {"stages_passed": [1, 2, 3], "checks_passed": len(checks),
          "scope": "Generated bioinformatics panels; bulk differential-expression fits are supplied inputs.",
          "panels": ["2D", "2E", "3D", "3E", "3F", "3G", "4A", "4B", "4C", "4D", "S1C", "S1D", "S1E", "S1G", "S2I", "S2J", "S3A", "S3B", "S3C"],
          "input_and_result_hashes": hashes}
(OUT / "validation_summary.json").write_text(json.dumps(report, indent=2) + "\n")
print("Validation stages passed: 3; checks passed:", len(checks))
