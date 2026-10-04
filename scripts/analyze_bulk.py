"""Prepare bulk RNA-seq gene lists from supplied differential-expression tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/source/figure4"
OUT = ROOT / "results/tables/bulk"
OUT.mkdir(parents=True, exist_ok=True)
CONFIG = json.loads((ROOT / "config/analysis.json").read_text())["bulk"]
corrections = pd.read_csv(ROOT / "config/gene_identifier_corrections.csv")
identifier_map = dict(zip(corrections.supplied_Tcell_id, corrections.corresponding_cytokine_id))

tables, summary = {}, {}
for condition, filename in [("T_cells", "RNAseq_filtered_without_tcell_genes.csv"),
                            ("TNF_IFNg", "Cytokines.xlsx")]:
    path = DATA / filename
    table = pd.read_csv(path, sep=";") if path.suffix == ".csv" else pd.read_excel(path, sheet_name="WT vs TNF-IFN")
    if list(table) != ["id", "baseMean", "log2FoldChange", "pval", "padj"]:
        raise ValueError("Unexpected columns in " + filename)
    table = table.rename(columns={"id": "source_gene_id"})
    table["gene"] = table.source_gene_id.astype(str).str.strip().replace(identifier_map)
    if table.gene.isna().any() or table.gene.duplicated().any():
        raise ValueError("Gene identifiers must be present and unique: " + condition)
    effects = table.log2FoldChange.astype("string").str.replace("\u2013", "-", regex=False)
    table["log2FoldChange"] = pd.to_numeric(effects, errors="raise").astype(float)
    for field in ("baseMean", "pval", "padj"):
        table[field] = pd.to_numeric(table[field], errors="raise")
    if (table.baseMean.dropna() < 0).any():
        raise ValueError("Negative mean expression: " + condition)
    for field in ("pval", "padj"):
        if not table[field].dropna().between(0, 1).all():
            raise ValueError("Invalid probability in " + field)
    finite = np.isfinite(table.log2FoldChange)
    table["effect_status"] = np.select(
        [finite, np.isposinf(table.log2FoldChange), np.isneginf(table.log2FoldChange)],
        ["finite", "positive_infinity", "negative_infinity"], default="missing")
    table["eligible_for_gene_lists"] = finite & table.padj.notna() & (table.baseMean >= CONFIG["base_mean_min"])
    significant = table.eligible_for_gene_lists & (table.padj < CONFIG["adjusted_p_max"])
    table["upregulated"] = significant & (table.log2FoldChange > CONFIG["absolute_log2_fold_change_min"])
    table["downregulated"] = significant & (table.log2FoldChange < -CONFIG["absolute_log2_fold_change_min"])
    table.to_csv(OUT / (condition + "_statistics.csv"), index=False)
    table.loc[table.upregulated | table.downregulated].to_csv(OUT / (condition + "_differential_genes.csv"), index=False)
    table.loc[~finite].to_csv(OUT / (condition + "_nonfinite_effects.csv"), index=False)
    rankable = (table.baseMean >= CONFIG["gsea_base_mean_min"]) & table.log2FoldChange.notna() & table.pval.gt(0) & table.pval.le(1)
    ranking = table.loc[rankable, ["gene", "log2FoldChange", "pval", "baseMean"]].copy()
    ranking["rank_score"] = np.sign(ranking.log2FoldChange) * -np.log10(ranking.pval)
    ranking.sort_values(["rank_score", "gene"], ascending=[False, True], kind="stable").to_csv(OUT / (condition + "_gsea_ranks.csv"), index=False)
    tables[condition] = table.set_index("gene", drop=False)
    summary[condition] = {
        "supplied_rows": len(table), "finite_effects": int(finite.sum()),
        "nonfinite_effects": int((~finite).sum()), "eligible_for_gene_lists": int(table.eligible_for_gene_lists.sum()),
        "upregulated": int(table.upregulated.sum()), "downregulated": int(table.downregulated.sum()),
        "gsea_ranked_genes": len(ranking), "gene_identifier_corrections": int(table.gene.ne(table.source_gene_id).sum())}

up = {condition: set(table.index[table.upregulated]) for condition, table in tables.items()}
groups = {"Common": up["T_cells"] & up["TNF_IFNg"],
          "T_cells_only": up["T_cells"] - up["TNF_IFNg"],
          "TNF_IFNg_only": up["TNF_IFNg"] - up["T_cells"]}
membership = []
for group, genes in groups.items():
    for gene in sorted(genes):
        row = {"gene": gene, "group": group}
        for condition, table in tables.items():
            for field in ("baseMean", "log2FoldChange", "pval", "padj"):
                row[condition + "_" + field] = table.at[gene, field] if gene in table.index else np.nan
        membership.append(row)
membership = pd.DataFrame(membership)
membership.to_csv(OUT / "upregulated_gene_membership.csv", index=False)
common = membership[membership.group.eq("Common")].copy()
common["mean_log2FoldChange"] = common[["T_cells_log2FoldChange", "TNF_IFNg_log2FoldChange"]].mean(axis=1)
top = common.sort_values(["mean_log2FoldChange", "gene"], ascending=[False, True], kind="stable").head(30)
top.to_csv(OUT / "figure_4c_top_30_shared_genes.csv", index=False)

eligible = {condition: set(table.index[table.eligible_for_gene_lists]) for condition, table in tables.items()}
universes = {"Common": eligible["T_cells"] & eligible["TNF_IFNg"],
             "T_cells_only": eligible["T_cells"], "TNF_IFNg_only": eligible["TNF_IFNg"]}
for group, universe in universes.items():
    if not groups[group] <= universe:
        raise ValueError("Gene list is outside its tested background: " + group)
    pd.DataFrame({"gene": sorted(universe)}).to_csv(OUT / (group + "_enrichment_universe.csv"), index=False)
summary["upregulated_gene_sets"] = {name: len(genes) for name, genes in groups.items()}
summary["enrichment_universe_genes"] = {name: len(genes) for name, genes in universes.items()}
summary["parameters"] = CONFIG
(OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
