"""Render manuscript panels from validated analysis tables."""
from pathlib import Path
import hashlib
import json
import math
import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle
import numpy as np
import pandas as pd
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results/tables"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
validation = json.loads((TABLES / "validation/validation_summary.json").read_text())
if validation["stages_passed"] != [1, 2, 3]:
    raise ValueError("Complete numerical validation before rendering figures")
for name, digest in validation["input_and_result_hashes"].items():
    if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
        raise ValueError("An analysis file changed after validation: " + name)
config = json.loads((ROOT / "config/analysis.json").read_text())
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 13,
                     "axes.titleweight": "bold", "axes.labelsize": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "svg.fonttype": "none", "savefig.facecolor": "white"})
CMAP = LinearSegmentedColormap.from_list("expression", ["#3B4CC0", "#FFFFFF", "#B40426"])
COLORS = {"Tcells_noIL2": "#F8766D", "Tcells_IL2": "#00BA38", "HeLa_coculture": "#619CFF"}
CONDITIONS = list(COLORS)
LABELS = {"Tcells_noIL2": "T cells, no IL-2", "Tcells_IL2": "T cells, IL-2", "HeLa_coculture": "HeLa co-culture"}
STATE_COLORS = {0: "#CD9600", 1: "#F8766D", 2: "#00BFC4", 3: "#00B0F6", 4: "#C77CFF", 5: "#7CAE00", 6: "#F564E3", 7: "#00BE67"}
DISPLAY_MODELS = {"A431": "A431", "H460": "H460", "HeLa": "HeLa", "HepG2": "HepG2", "SW837": "SW837", "T47D": "T-47D"}
saved = []

def save(fig, name, description):
    pdf_temp = OUT / (name + ".pdf.tmp")
    png_temp = OUT / (name + ".png.tmp")
    fig.savefig(pdf_temp, format="pdf", bbox_inches="tight", metadata={"Title": description, "Author": "", "Subject": "Manuscript figure", "Keywords": "", "Creator": "Matplotlib", "CreationDate": None, "ModDate": None})
    if not pdf_temp.read_bytes().rstrip().endswith(b"%%EOF"):
        raise ValueError("Incomplete PDF export: " + name)
    pdf_temp.replace(OUT / (name + ".pdf"))
    fig.savefig(png_temp, format="png", dpi=600, bbox_inches="tight", metadata={"Title": description, "Software": "Matplotlib"})
    png_temp.replace(OUT / (name + ".png"))
    saved.append({"figure": name, "description": description, "formats": "PDF; PNG", "png_dpi": 600})
    plt.close(fig)

def panel_label(ax, label):
    ax.text(-.12, 1.06, label, transform=ax.transAxes, fontsize=20, weight="bold", va="bottom")

def mark_models(ax, data, x, y):
    positions = {"A431": (.79, .82), "H460": (.17, .24), "HeLa": (.76, .43),
                 "HepG2": (.74, .18), "SW837": (.18, .81), "T47D": (.15, .58)}
    for label, identifier in config["panel_models"].items():
        row = data.loc[data.ModelID.eq(identifier)].iloc[0]
        ax.scatter(row[x], row[y], s=38, color="#E67E22", edgecolor="white", linewidth=.6, zorder=4)
        ax.annotate(DISPLAY_MODELS[label], (row[x], row[y]), xytext=positions[label], textcoords="axes fraction", fontsize=9,
                    weight="bold", bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "#999999", "lw": .5},
                    arrowprops={"arrowstyle": "-", "lw": .5, "color": "#666666"}, zorder=5)

public = TABLES / "public_data"
summary = json.loads((public / "analysis_summary.json").read_text())
timer = pd.read_csv(public / "timer_pathway_summaries.csv")
nonmel = timer[timer.cohort_group.eq("Non-melanoma")].set_index("pathway").loc[["Apoptosis", "Pyroptosis", "Necroptosis"]]
fig, ax = plt.subplots(figsize=(4.6, 3.7), layout="constrained")
height = -np.log10(nonmel.fisher_p_nominal)
ax.bar(nonmel.index, height, color="#555555", width=.64)
for i, row in enumerate(nonmel.itertuples()):
    ax.text(i, height.iloc[i]+.10, "P = " + format(row.fisher_p_nominal, ".3g"), ha="center", fontsize=10)
ax.set(ylim=(0, height.max()+.65), ylabel="−log₁₀ nominal Fisher P", title="Clinical cohort associations")
ax.text(.02, .97, "Non-melanoma cohorts", transform=ax.transAxes, va="top", fontsize=9)
save(fig, "Figure2D", "TIMER clinical cohort associations")

z = pd.read_csv(public / "figure_2e_expression_z_scores.csv").set_index("ModelID")
gene_groups = config["figure_2e_gene_groups"]
gene_order = sum(gene_groups.values(), [])
model_order = list(config["panel_models"])
values = z.loc[[config["panel_models"][model] for model in model_order], gene_order].to_numpy()
fig, ax = plt.subplots(figsize=(7.2, 3.0))
fig.subplots_adjust(left=.07, right=.95, bottom=.30, top=.79)
im = ax.imshow(values, cmap=CMAP, vmin=-2, vmax=2, aspect="auto")
ax.set_xticks(range(len(gene_order)), gene_order, rotation=90, fontsize=6.5)
ax.set_yticks(range(6), [DISPLAY_MODELS[x] for x in model_order], fontsize=8)
ax.tick_params(length=0)
titles = {"TNF_IFN_input": "TNF/IFNγ\nsensing", "Apoptosis": "Death receptor /\napoptosis", "Necroptosis": "Necroptosis",
          "Pyroptosis_inflammasome": "Pyroptosis /\ninflammasome", "Survival_brakes": "Survival\nbrakes"}
start = 0
for name, genes in gene_groups.items():
    ax.add_patch(Rectangle((start-.5, -1.25), len(genes), .74, facecolor="#E5E5E5", edgecolor="#555555", lw=.6, clip_on=False))
    ax.text(start+(len(genes)-1)/2, -.88, titles[name], ha="center", va="center", fontsize=6.0 if name == "Necroptosis" else 6.6, weight="bold")
    ax.axvline(start-.5, color="#555555", lw=.8)
    start += len(genes)
ax.set_xticks(np.arange(-.5, len(gene_order), 1), minor=True)
ax.set_yticks(np.arange(-.5, 6, 1), minor=True)
ax.grid(which="minor", color="white", linewidth=.6)
ax.tick_params(which="minor", length=0)
fig.suptitle("Baseline cell-death pathway expression", y=1.02, fontsize=12, weight="bold")
cax = fig.add_axes([.962, .32, .013, .42])
fig.colorbar(im, cax=cax, extend="both", ticks=[-2, -1, 0, 1, 2]).set_label("Gene-wise z-score", fontsize=8)
cax.tick_params(labelsize=7)
save(fig, "Figure2E", "DepMap expression across the experimental cell-line panel")

sc = pd.read_csv(TABLES / "single_cell/cell_metadata_and_coordinates.csv")
order = np.random.default_rng(123).permutation(len(sc))
shuffled = sc.iloc[order]
fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.8))
fig.subplots_adjust(left=.06, right=.98, top=.89, bottom=.32, wspace=.19)
for ax, field, colors, title, letter in [(axes[0], "condition", COLORS, "CD3+ T cells by condition", "D"),
                                       (axes[1], "seurat_clusters", STATE_COLORS, "CD3+ T-cell states", "E")]:
    ax.scatter(shuffled.UMAP_1, shuffled.UMAP_2, c=shuffled[field].map(colors), s=2.3, alpha=.75, linewidths=0, rasterized=True)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", title=title)
    ax.set_xticks([]); ax.set_yticks([])
    panel_label(ax, letter)
condition_handles = [Line2D([], [], marker="o", linestyle="", color=COLORS[k], label=LABELS[k], markersize=5) for k in CONDITIONS]
axes[0].legend(handles=condition_handles, fontsize=9, frameon=False, loc="upper left")
state_handles = [Line2D([], [], marker="o", linestyle="", color=STATE_COLORS[int(key)],
                        label=textwrap.fill(value, 39), markersize=5) for key, value in config["cell_states"].items()]
fig.legend(handles=state_handles, ncol=2, fontsize=9, frameon=False, loc="lower center", bbox_to_anchor=(.53, -.015),
           columnspacing=2.2, labelspacing=.75, handletextpad=.45)
save(fig, "Figure3DE", "Single-cell UMAP by condition and annotated state")

composition = pd.crosstab(sc.condition, sc.seurat_clusters, normalize="index").loc[CONDITIONS] * 100
fig, ax = plt.subplots(figsize=(7.6, 4.5))
fig.subplots_adjust(left=.1, right=.55, bottom=.25, top=.9)
bottom = np.zeros(3)
for cluster in range(8):
    height = composition[cluster].to_numpy()
    ax.bar(range(3), height, bottom=bottom, color=STATE_COLORS[cluster], width=.75, edgecolor="white", linewidth=.3)
    bottom += height
ax.set_xticks(range(3), [LABELS[c] for c in CONDITIONS], rotation=30, ha="right")
ax.set(ylabel="Cells (%)", ylim=(0, 100), title="CD3+ T-cell state composition")
ax.legend(handles=state_handles, bbox_to_anchor=(1.03, 1.03), loc="upper left", frameon=False, fontsize=8.8,
          labelspacing=.75, handletextpad=.4)
save(fig, "Figure3F", "Single-cell state composition by condition")

module_order = ["Cytotoxicity", "FAS_FASL_axis", "Naive_memory", "Proliferation", "TCR_activation"]
module_labels = ["Cytotoxicity", "FAS/FASL axis", "Naive/memory", "Proliferation", "TCR activation"]
condition_means = sc.groupby("condition")[module_order].mean().loc[CONDITIONS]
module_z = ((condition_means-condition_means.mean())/condition_means.std(ddof=1)).T
module_z.to_csv(OUT / "Figure3G_source_values.csv")
fig, ax = plt.subplots(figsize=(5.1, 4.4), layout="constrained")
im = ax.imshow(module_z, cmap=CMAP, vmin=-1.2, vmax=1.2, aspect="auto")
ax.set_xticks(range(3), [LABELS[c] for c in CONDITIONS], rotation=35, ha="right", fontsize=10)
ax.set_yticks(range(5), module_labels)
ax.set_title("T-cell module scores")
for y in range(5):
    for x in range(3):
        v = module_z.iloc[y, x]
        ax.text(x, y, f"{v:.1f}", ha="center", va="center", color="white" if abs(v) > .85 else "black", fontsize=11)
fig.colorbar(im, ax=ax, shrink=.68, pad=.04).set_label("Z-score across conditions")
save(fig, "Figure3G", "T-cell module scores standardized across conditions")

fig, axes = plt.subplots(2, 1, figsize=(4.3, 6.7), layout="constrained")
for ax, condition, title in zip(axes, ["TNF_IFNg", "T_cells"], ["TNF+IFNγ", "T-cell co-culture"]):
    table = pd.read_csv(TABLES / "bulk" / (condition + "_statistics.csv"))
    table = table[table.eligible_for_gene_lists & table.padj.gt(0)].copy()
    table["y"] = -np.log10(table.padj)
    for name, mask, color in [("Other tested genes", ~(table.upregulated | table.downregulated), "#B9B9B9"),
                              ("Downregulated", table.downregulated, "#287BB5"), ("Upregulated", table.upregulated, "#D64B3D")]:
        ax.scatter(table.loc[mask, "log2FoldChange"], table.loc[mask, "y"], color=color, s=4, alpha=.8, linewidths=0,
                   label=f"{name} (n={int(mask.sum())})", rasterized=True)
    ax.axvline(-1, color="#777777", linestyle="--", linewidth=.6); ax.axvline(1, color="#777777", linestyle="--", linewidth=.6)
    ax.axhline(-np.log10(.05), color="#777777", linestyle="--", linewidth=.6)
    annotations = table[table.gene.isin(["AIM2", "MLKL", "CASP1", "CASP7", "CASP8", "STAT1"])].sort_values("y")
    label_y = 0
    for row in annotations.itertuples(index=False):
        label_y = max(row.y + .8, label_y + 1.7)
        ax.annotate(row.gene, (row.log2FoldChange, row.y), xytext=(10.6, label_y),
                    textcoords="data", fontsize=8, va="center",
                    bbox={"fc": "white", "ec": "none", "pad": .15},
                    arrowprops={"arrowstyle": "-", "lw": .45, "color": "#555555"})
    ax.set(xlabel="log₂ fold change", ylabel="−log₁₀ adjusted P", title=title, ylim=(0, table.y.max()*1.12))
    ax.legend(frameon=False, loc="upper left", fontsize=6.7, handlelength=1)
save(fig, "Figure4A", "Volcano plots from full supplied differential-expression tables")

membership = pd.read_csv(TABLES / "bulk/upregulated_gene_membership.csv")
counts = membership.group.value_counts()
left, shared, right = int(counts.TNF_IFNg_only), int(counts.Common), int(counts.T_cells_only)
r1, r2 = math.sqrt((left+shared)/math.pi), math.sqrt((right+shared)/math.pi)
def overlap(d):
    a = np.clip((d*d+r1*r1-r2*r2)/(2*d*r1), -1, 1)
    b = np.clip((d*d+r2*r2-r1*r1)/(2*d*r2), -1, 1)
    return r1*r1*np.arccos(a) + r2*r2*np.arccos(b) - .5*np.sqrt(max(0, (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2)))
distance = brentq(lambda d: overlap(d)-shared, abs(r1-r2)+1e-8, r1+r2-1e-8)
fig, ax = plt.subplots(figsize=(5.0, 3.8), layout="constrained")
ax.add_patch(Circle((0, 0), r1, color="#FF705C", alpha=.75, ec="black", lw=1))
ax.add_patch(Circle((distance, 0), r2, color="#FFD92F", alpha=.75, ec="black", lw=1))
for x, value in [(-r1*.36, left), ((r1+distance-r2)/2, shared), (distance+r2*.4, right)]:
    ax.text(x, 0, str(value), ha="center", va="center", fontsize=18, weight="bold")
ax.text(-r1*.14, r1+1.4, "TNF+IFNγ", ha="center", fontsize=12, weight="bold")
ax.text(distance, r2+1.4, "T cells", ha="center", fontsize=12, weight="bold")
ax.set(xlim=(-r1-1, distance+r2+1), ylim=(-r1-2, r1+4), aspect="equal")
ax.axis("off")
save(fig, "Figure4B", "Shared and condition-specific upregulated genes")

top = pd.read_csv(TABLES / "bulk/figure_4c_top_30_shared_genes.csv")
heat = top[["TNF_IFNg_log2FoldChange", "T_cells_log2FoldChange"]].to_numpy()
fig, ax = plt.subplots(figsize=(4.1, 7.2), layout="constrained")
im = ax.imshow(heat, cmap="Reds", vmin=0, vmax=np.ceil(heat.max()), aspect="auto")
ax.set_xticks([0, 1], ["TNF+IFNγ", "T cells"], fontsize=11)
ax.set_yticks(range(len(top)), top.gene, fontsize=9)
ax.tick_params(length=0)
ax.set_title("Top 30 shared upregulated genes", fontsize=12)
fig.colorbar(im, ax=ax, pad=.04, shrink=.45).set_label("log₂ fold change")
save(fig, "Figure4C", "Top shared upregulated genes on an absolute fold-change scale")

gsea = pd.read_csv(TABLES / "enrichment/hallmark_gsea_results.csv")
curves = pd.read_csv(TABLES / "enrichment/hallmark_gsea_running_scores.csv.gz")
pathways = ["HALLMARK_APOPTOSIS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE"]
titles = ["Apoptosis", "TNF signaling via NF-κB", "IFNγ response"]
fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.5), layout="constrained")
for y, condition in enumerate(["TNF_IFNg", "T_cells"]):
    for x, (pathway, title) in enumerate(zip(pathways, titles)):
        ax = axes[y, x]
        data = curves[curves.condition.eq(condition) & curves.pathway.eq(pathway)]
        result = gsea[gsea.condition.eq(condition) & gsea.pathway.eq(pathway)].iloc[0]
        ax.plot(data["rank"], data.running_enrichment_score, color="#168B37", linewidth=1.4)
        ax.axhline(0, color="#888888", linewidth=.55)
        hits = data.loc[data.member, "rank"]
        ax.vlines(hits, -.09, -.035, color="#555555", linewidth=.28)
        ax.set(xlabel="Gene rank", ylim=(-.12, .95), title=title)
        ax.text(.03, .97, f"NES = {result.NES:.2f}\nAdjusted P = {result.padj:.3g}", transform=ax.transAxes, va="top", fontsize=9)
        if x == 0:
            ax.set_ylabel(("TNF+IFNγ" if y == 0 else "T cells") + "\nEnrichment score")
        ax.tick_params(labelsize=8)
save(fig, "Figure4D", "Preranked enrichment of three selected Hallmark pathways")

scores = pd.read_csv(public / "figure_s1c_pathway_scores.csv")
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9), layout="constrained")
lo = min(scores[x].min() for x in ["Apoptosis", "Pyroptosis", "Necroptosis"])-.15
hi = max(scores[x].max() for x in ["Apoptosis", "Pyroptosis", "Necroptosis"])+.15
for ax, x, y in zip(axes, ["Pyroptosis", "Necroptosis", "Pyroptosis"], ["Apoptosis", "Apoptosis", "Necroptosis"]):
    ax.scatter(scores[x], scores[y], s=6, color="#7898B6", alpha=.36, linewidths=0, rasterized=True)
    mark_models(ax, scores, x, y)
    ax.axhline(0, color="#555555", linestyle="--", lw=.6); ax.axvline(0, color="#555555", linestyle="--", lw=.6)
    ax.set(xlabel=x+" score (mean z)", ylabel=y+" score (mean z)", xlim=(lo, hi), ylim=(lo, hi))
    ax.set_title(x + " vs " + y, fontsize=11)
save(fig, "FigureS1C", "Cell-death pathway scores across all 1343 selected DepMap models")

loco = pd.read_csv(public / "timer_leave_one_cohort_out.csv")
loco = loco[loco.cohort_group.eq("Non-melanoma")]
fig, ax = plt.subplots(figsize=(5.3, 3.6), layout="constrained")
for y, pathway in enumerate(["Apoptosis", "Pyroptosis", "Necroptosis"]):
    points = loco[loco.pathway.eq(pathway)].signed_z_nominal
    ax.scatter(points, y + np.linspace(-.065, .065, len(points)), color="#444444", s=25)
    ax.scatter(nonmel.loc[pathway, "signed_z_nominal"], y, color="#C83E36", marker="D", s=28, zorder=4)
ax.axvline(0, color="black", linestyle="--", lw=.8)
ax.set_yticks(range(3), ["Apoptosis", "Pyroptosis", "Necroptosis"])
ax.invert_yaxis()
ax.set(xlabel="Combined signed Z (nominal)", title="Leave-one-cohort-out analysis")
ax.legend(handles=[Line2D([], [], color="#444444", marker="o", ls="", label="One cohort omitted"),
                   Line2D([], [], color="#C83E36", marker="D", ls="", label="All available cohorts")],
          fontsize=8, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(.5, -.22))
save(fig, "FigureS1D", "TIMER non-melanoma leave-one-cohort-out summaries")

scores = pd.read_csv(public / "figure_s1e_pathway_scores.csv")
fig, ax = plt.subplots(figsize=(5.7, 4.5), layout="constrained")
ax.scatter(scores.TNF_IFN_input, scores.Death_execution, s=8, color="#7898B6", alpha=.4, linewidths=0, rasterized=True)
mark_models(ax, scores, "TNF_IFN_input", "Death_execution")
ax.axvline(0, color="#555555", linestyle="--", lw=.7); ax.axhline(0, color="#555555", linestyle="--", lw=.7)
ax.set(xlabel="TNF/IFNγ input score (mean z)", ylabel="Death execution score (mean z)", title="Inflammatory death-pathway landscape")
save(fig, "FigureS1E", "DepMap inflammatory sensing and death-execution scores")

pairs = pd.read_csv(public / "figure_s1g_paired_auc.csv")
fig, ax = plt.subplots(figsize=(5.3, 4.4), layout="constrained")
violins = ax.violinplot([pairs.HeLa, pairs["NCI-H460"]], positions=[0, 1], showmedians=True, showextrema=False, widths=.75)
for body in violins["bodies"]:
    body.set_facecolor("#CCCCCC"); body.set_edgecolor("black"); body.set_alpha(.7)
violins["cmedians"].set_color("black")
for row in pairs.itertuples(index=False, name=None):
    ax.plot([0, 1], [row[1], row[2]], color="#888888", linewidth=.45, alpha=.6)
ax.scatter(np.zeros(len(pairs)), pairs.HeLa, s=12, color="#386CB0", zorder=3)
ax.scatter(np.ones(len(pairs)), pairs["NCI-H460"], s=12, color="#E67E22", zorder=3)
ax.set_xticks([0, 1], ["HeLa", "NCI-H460"])
ax.set(ylabel="CTRPv2 response-curve AUC", title="DNMT/HDAC inhibitor responses")
ax.text(.04, .96, f"25 paired compounds\nΔAUC (H460 − HeLa) = {summary['mean_difference_H460_minus_HeLa']:.2f}\nPaired P = {summary['paired_p']:.4f}", transform=ax.transAxes, va="top", fontsize=9)
ax.set_ylim(min(pairs.HeLa.min(), pairs["NCI-H460"].min())-1, max(pairs.HeLa.max(), pairs["NCI-H460"].max())+3.7)
save(fig, "FigureS1G", "Paired drug-response comparison between HeLa and H460")

fig, axes = plt.subplots(2, 3, figsize=(9.3, 6.3), layout="constrained")
for ax, module, title in zip(axes.flat, module_order, module_labels):
    vectors = [sc.loc[sc.condition.eq(condition), module] for condition in CONDITIONS]
    box = ax.boxplot(vectors, tick_labels=["No IL-2", "IL-2", "HeLa\nco-culture"], patch_artist=True,
                     flierprops={"markersize": 1, "markeredgecolor": "#888888", "alpha": .2}, widths=.6)
    for patch, condition in zip(box["boxes"], CONDITIONS):
        patch.set_facecolor(COLORS[condition]); patch.set_alpha(.75)
    for median in box["medians"]: median.set_color("black")
    ax.set(title=title, ylabel="Mean log-normalized expression")
    ax.tick_params(labelsize=9)
axes.flat[-1].axis("off")
axes.flat[-1].text(.05, .80, "CD3+ cells\n\nNo IL-2: 4,051\nIL-2: 1,146\nHeLa co-culture: 5,935", transform=axes.flat[-1].transAxes, va="top", fontsize=11)
save(fig, "FigureS2I", "Single-cell distributions of T-cell module scores")

fig, axes = plt.subplots(1, 3, figsize=(11.3, 5.1))
fig.subplots_adjust(left=.055, right=.99, bottom=.37, top=.89, wspace=.13)
for ax, condition in zip(axes, CONDITIONS):
    frame = shuffled[shuffled.condition.eq(condition)]
    ax.scatter(frame.UMAP_1, frame.UMAP_2, c=frame.seurat_clusters.map(STATE_COLORS), s=3.2, linewidths=0, alpha=.8, rasterized=True)
    ax.set(xlabel="UMAP 1", title=LABELS[condition], xlim=(sc.UMAP_1.min()-.4, sc.UMAP_1.max()+.4), ylim=(sc.UMAP_2.min()-.4, sc.UMAP_2.max()+.4))
    ax.set_xticks([]); ax.set_yticks([])
axes[0].set_ylabel("UMAP 2")
fig.legend(handles=state_handles, ncol=2, fontsize=9, loc="lower center", frameon=False, bbox_to_anchor=(.52, -.02), labelspacing=.75)
save(fig, "FigureS2J", "Single-cell states across culture conditions")

enrichment = pd.read_csv(TABLES / "enrichment/figure_s3_displayed_terms.csv")
for group, name, title in [("Common", "FigureS3A", "Shared upregulated genes"),
                          ("T_cells_only", "FigureS3B", "Genes upregulated only in the T-cell comparison"),
                          ("TNF_IFNg_only", "FigureS3C", "Genes upregulated only in the TNF+IFNγ comparison")]:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 10.7), layout="constrained")
    for ax, collection, heading in zip(axes.flat, ["KEGG", "GO_BP", "GO_MF", "GO_CC"], ["KEGG pathways", "GO biological process", "GO molecular function", "GO cellular component"]):
        frame = enrichment[enrichment.group.eq(group) & enrichment.collection.eq(collection)].iloc[::-1]
        if frame.empty:
            ax.text(.5, .5, "No terms with adjusted P < 0.05", ha="center", transform=ax.transAxes); ax.set_axis_off(); continue
        colors = -np.log10(frame.padj)
        scatter = ax.scatter(frame.gene_ratio, range(len(frame)), s=25 + frame["count"]*4, c=colors, cmap="viridis", edgecolors="#444444", linewidth=.35)
        ax.set_yticks(range(len(frame)), [textwrap.fill(term, 31) for term in frame.term_name], fontsize=10)
        ax.set(xlabel="Gene ratio", title=heading, ylim=(-.8, len(frame)-.2))
        ax.grid(axis="x", color="#DDDDDD", linewidth=.5); ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
        colorbar = fig.colorbar(scatter, ax=ax, shrink=.50, pad=.02)
        colorbar.set_label("−log₁₀ adjusted P", fontsize=9)
        colorbar.ax.tick_params(labelsize=8)
    count_handles = [Line2D([], [], marker="o", linestyle="", markersize=np.sqrt(25+4*n),
                            markerfacecolor="#999999", markeredgecolor="#444444", markeredgewidth=.35,
                            label=str(n)) for n in [5, 20, 50, 100]]
    fig.legend(handles=count_handles, title="Genes per term", ncol=4, loc="outside lower center", frameon=False, fontsize=9)
    fig.suptitle(title + "\nbaseMean ≥ 30; log₂ fold change > 1; adjusted P < 0.05", fontsize=13, weight="bold")
    save(fig, name, title + ": GO and KEGG over-representation")

pd.DataFrame(saved).to_csv(OUT / "figure_manifest.csv", index=False)
print("Rendered figure files:", len(saved), "PDF and", len(saved), "PNG")
