"""Calculate DepMap expression scores, TIMER summaries and paired CTRPv2 results."""
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/source"
OUT = ROOT / "results/tables/public_data"
OUT.mkdir(parents=True, exist_ok=True)
config = json.loads((ROOT / "config/analysis.json").read_text())
modules, core = config["depmap_modules"], config["core_death_modules"]
heatmap_groups = config["figure_2e_gene_groups"]
genes = sorted(set(sum(modules.values(), []) + sum(core.values(), []) + sum(heatmap_groups.values(), [])))
expression_path = DATA / "depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
headers = pd.read_csv(expression_path, nrows=0).columns
columns = [column for column in headers if column in ["ModelID", "SequencingID", "IsDefaultEntryForModel"] or column.split(" (")[0] in genes]
expression = pd.read_csv(expression_path, usecols=columns)
expression.to_csv(OUT / "depmap_selected_expression.csv", index=False)
gene_columns = [column for column in columns if re.search(r" \(\d+\)$", column)]
if set(genes) != {column.split(" (")[0] for column in gene_columns}:
    raise ValueError("The expression file does not contain every required gene exactly once")
metadata = pd.read_csv(DATA / "depmap/Model.csv")
defaults = expression[expression.IsDefaultEntryForModel.eq("Yes")].copy()
joined = defaults.merge(metadata, on="ModelID", how="left", validate="one_to_one", indicator=True)
if not joined._merge.eq("both").all():
    raise ValueError("Unmatched default expression profiles")
keep = joined.OncotreeLineage.notna() & joined.StrippedCellLineName.notna()
keep &= ~joined.OncotreeLineage.isin(["Lymphoid", "Myeloid", "Fibroblast", "Embryonal", "Other"])
keep &= ~joined.OncotreePrimaryDisease.str.contains("Non-Cancerous|Unknown", case=False, na=False)
solid = joined.loc[keep].copy()
matrix = solid[gene_columns].copy()
matrix.columns = [column.split(" (")[0] for column in gene_columns]
if not np.isfinite(matrix.to_numpy()).all() or (matrix.std(ddof=1) == 0).any():
    raise ValueError("Invalid or constant expression feature")
z = (matrix - matrix.mean()) / matrix.std(ddof=1)
identity = solid[["ModelID", "SequencingID", "CellLineName", "StrippedCellLineName", "OncotreeLineage", "OncotreePrimaryDisease"]]
identity.to_csv(OUT / "depmap_analysis_cohort.csv", index=False)
pd.concat([identity, z], axis=1).to_csv(OUT / "depmap_gene_z_scores.csv", index=False)
panel = solid.ModelID.isin(config["panel_models"].values())
if panel.sum() != 6:
    raise ValueError("Expected six distinct experimental models")
pd.concat([identity.loc[panel], z.loc[panel]], axis=1).to_csv(OUT / "figure_2e_expression_z_scores.csv", index=False)
for signatures, filename in [(core, "figure_s1c_pathway_scores.csv"), (modules, "figure_s1e_pathway_scores.csv")]:
    scores = identity.copy()
    for label, members in signatures.items():
        scores[label] = z[members].mean(axis=1)
    if signatures is modules:
        scores["Death_execution"] = scores[["Apoptosis", "Necroptosis", "Pyroptosis_inflammasome"]].mean(axis=1)
    scores.to_csv(OUT / filename, index=False)

export = pd.read_csv(DATA / "figure2/Immunotherapy_outcome_PANoptosis.csv")
numeric = pd.DataFrame({"cohort": export["heatTable.cancer"], "gene": export["heatTable.variable"],
                        "cox_z": pd.to_numeric(export["heatTable.p"], errors="coerce"),
                        "pvalue": pd.to_numeric(export["heatTable.rho"], errors="coerce")}).dropna(subset=["cox_z", "pvalue"])
if numeric.duplicated(["cohort", "gene"]).any() or not numeric.pvalue.between(0, 1, inclusive="right").all():
    raise ValueError("Invalid TIMER summary records")
if not np.allclose(numeric.pvalue, 2*stats.norm.sf(numeric.cox_z.abs()), atol=1e-12, rtol=0):
    raise ValueError("The TIMER score and probability columns are inconsistent")
numeric["cohort_group"] = np.where(numeric.cohort.str.contains("skcm|melano", case=False), "Melanoma", "Non-melanoma")
numeric.to_csv(OUT / "timer_gene_cohort_statistics.csv", index=False)
pathway_rows, loco_rows = [], []
selected = numeric[numeric.gene.isin(sum(core.values(), []))]
for label, group in [("All cohorts", selected)] + list(selected.groupby("cohort_group")):
    for pathway, members in core.items():
        frame = group[group.gene.isin(members)]
        pathway_rows.append({"cohort_group": label, "pathway": pathway, "valid_gene_cohort_tests": len(frame),
                             "cohorts_with_data": frame.cohort.nunique(),
                             "fisher_p_nominal": float(stats.combine_pvalues(frame.pvalue).pvalue),
                             "signed_z_nominal": float(frame.cox_z.sum()/np.sqrt(len(frame)))})
        for omitted in sorted(group.cohort.unique()):
            subset = frame[frame.cohort.ne(omitted)]
            loco_rows.append({"cohort_group": label, "pathway": pathway, "omitted_cohort": omitted,
                              "valid_gene_cohort_tests": len(subset),
                              "signed_z_nominal": float(subset.cox_z.sum()/np.sqrt(len(subset)))})
pd.DataFrame(pathway_rows).to_csv(OUT / "timer_pathway_summaries.csv", index=False)
pd.DataFrame(loco_rows).to_csv(OUT / "timer_leave_one_cohort_out.csv", index=False)

compounds = pd.read_csv(DATA / "figureS1G/v22.meta.per_compound.txt", sep="\t")
compounds["Drug"] = [next((rule["label"] for rule in config["drug_name_patterns"] if re.search(rule["pattern"], str(name).lower())), None)
                     for name in compounds.cpd_name]
compound_map = compounds.dropna(subset=["Drug"])
auc = pd.read_csv(DATA / "figureS1G/v22.data.auc_sensitivities.txt", sep="\t")
cells = pd.read_csv(DATA / "figureS1G/v22.meta.per_cell_line.txt", sep="\t")
cell_subset = cells[cells.index_ccl.isin([123, 359, 362])]
cell_subset.to_csv(OUT / "ctrp_selected_cell_lines.csv", index=False)
selected = auc[auc.index_ccl.isin([123, 359, 362])].merge(compound_map, on="index_cpd", validate="many_to_one")
selected["line"] = np.where(selected.index_ccl.eq(123), "HeLa", "NCI-H460")
selected.to_csv(OUT / "ctrp_selected_source_records.csv", index=False)
means = selected.groupby(["Drug", "cpd_name", "index_cpd", "line"], as_index=False).area_under_curve.mean()
pairs = means.groupby(["Drug", "line"]).area_under_curve.mean().unstack().dropna()
pairs["difference_H460_minus_HeLa"] = pairs["NCI-H460"] - pairs.HeLa
pairs.to_csv(OUT / "figure_s1g_paired_auc.csv")
test = stats.ttest_rel(pairs.HeLa, pairs["NCI-H460"])
delta = pairs.difference_H460_minus_HeLa
summary = {"depmap_profiles": len(expression), "depmap_default_profiles": len(defaults), "depmap_analysis_models": len(solid),
           "timer_numeric_records": len(numeric), "timer_available_cohorts": int(numeric.cohort.nunique()),
           "ctrp_compound_pairs": len(pairs), "paired_t_HeLa_minus_H460": float(test.statistic), "paired_df": float(test.df),
           "paired_p": float(test.pvalue), "mean_auc_HeLa": float(pairs.HeLa.mean()), "mean_auc_H460": float(pairs["NCI-H460"].mean()),
           "mean_difference_H460_minus_HeLa": float(delta.mean()),
           "difference_95ci": list(stats.t.interval(.95, len(pairs)-1, loc=delta.mean(), scale=stats.sem(delta)))}
(OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
