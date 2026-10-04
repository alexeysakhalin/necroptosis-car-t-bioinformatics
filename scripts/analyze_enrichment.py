"""Run over-representation analyses with explicitly defined tested-gene backgrounds."""
from pathlib import Path
import argparse
import hashlib
import json
import urllib.request
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--kegg-cache", type=Path, default=ROOT / "data/cache/kegg")
parser.add_argument("--allow-new-kegg-snapshot", action="store_true",
                    help="Explicitly start a new analysis with different KEGG resources")
args = parser.parse_args()
args.kegg_cache.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "results/tables/enrichment"
OUT.mkdir(parents=True, exist_ok=True)
ANNOTATIONS = ROOT / "data/annotations"
mapping = pd.read_csv(ANNOTATIONS / "gene_symbol_to_entrez.csv", dtype=str)
symbol_to_ids = mapping.groupby("SYMBOL").ENTREZID.agg(set).to_dict()
id_to_symbols = mapping.groupby("ENTREZID").SYMBOL.agg(lambda x: ";".join(sorted(set(x)))).to_dict()
memberships = pd.read_csv(ANNOTATIONS / "go_memberships.csv.gz", dtype=str)
terms = pd.read_csv(ANNOTATIONS / "go_terms.csv", dtype=str)
term_names = dict(zip(terms.GOID, terms.TERM))
collections = {}
for ontology in ("BP", "MF", "CC"):
    frame = memberships[memberships.ONTOLOGYALL.eq(ontology)]
    collections["GO_" + ontology] = frame.groupby("GOALL").ENTREZID.agg(set).to_dict()

resources = []
expected = {row["filename"]: row["sha256"] for row in
            json.loads((ROOT / "config/kegg_reference_manifest.json").read_text())}
for filename, url in [("kegg_hsa_links.tsv", "https://rest.kegg.jp/link/pathway/hsa"),
                       ("kegg_hsa_pathways.tsv", "https://rest.kegg.jp/list/pathway/hsa"),
                       ("kegg_release.txt", "https://rest.kegg.jp/info/kegg")]:
    path = args.kegg_cache / filename
    if not path.exists():
        with urllib.request.urlopen(url, timeout=60) as response:
            path.write_bytes(response.read())
    resources.append({"filename": filename, "url": url, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    if resources[-1]["sha256"] != expected[filename] and not args.allow_new_kegg_snapshot:
        raise ValueError("KEGG resource differs from the reference snapshot: " + filename +
                         ". Supply the recorded snapshot or explicitly use --allow-new-kegg-snapshot "
                         "for a new analysis. Published results have not been overwritten.")
kegg = pd.read_csv(args.kegg_cache / "kegg_hsa_links.tsv", sep="\t", names=["gene", "pathway"], dtype=str)
kegg["gene"] = kegg.gene.str.removeprefix("hsa:")
kegg["pathway"] = kegg.pathway.str.removeprefix("path:")
collections["KEGG"] = kegg.groupby("pathway").gene.agg(set).to_dict()
kegg_names = pd.read_csv(args.kegg_cache / "kegg_hsa_pathways.tsv", sep="\t", names=["pathway", "name"], dtype=str)
term_names.update(dict(zip(kegg_names.pathway.str.removeprefix("path:"), kegg_names.name.str.replace(" - Homo sapiens (human)", "", regex=False))))

membership = pd.read_csv(ROOT / "results/tables/bulk/upregulated_gene_membership.csv")
all_results, diagnostics = [], []
for group in ("Common", "T_cells_only", "TNF_IFNg_only"):
    universe_symbols = set(pd.read_csv(ROOT / "results/tables/bulk" / (group + "_enrichment_universe.csv")).gene)
    query_symbols = set(membership.loc[membership.group.eq(group), "gene"])
    universe_ids = set().union(*(symbol_to_ids.get(gene, set()) for gene in universe_symbols))
    query_ids = set().union(*(symbol_to_ids.get(gene, set()) for gene in query_symbols))
    if not query_ids <= universe_ids:
        raise ValueError("Gene list is outside its mapped background")
    for collection, term_sets in collections.items():
        annotated = set().union(*term_sets.values())
        universe = universe_ids & annotated
        query = query_ids & universe
        N, n = len(universe), len(query)
        results = []
        for term_id, members in sorted(term_sets.items()):
            relevant = members & universe
            M = len(relevant)
            if not 10 <= M <= 500:
                continue
            overlap = query & relevant
            k = len(overlap)
            pvalue = float(stats.hypergeom.sf(k-1, N, M, n))
            results.append({"group": group, "collection": collection, "term_id": term_id,
                            "term_name": term_names[term_id], "count": k,
                            "query_annotated_genes": n, "background_term_genes": M, "background_annotated_genes": N,
                            "gene_ratio": k/n if n else np.nan, "pvalue": pvalue,
                            "overlap_entrez_ids": ";".join(sorted(overlap, key=int)),
                            "overlap_symbols": ";".join(id_to_symbols.get(gene, gene) for gene in sorted(overlap, key=int))})
        frame = pd.DataFrame(results)
        if not frame.empty:
            frame["padj"] = stats.false_discovery_control(frame.pvalue.to_numpy(), method="bh")
            frame = frame.sort_values(["padj", "pvalue", "term_id"], kind="stable")
            all_results.append(frame)
        diagnostics.append({"group": group, "collection": collection, "input_symbols": len(query_symbols),
                            "mapped_entrez_ids": len(query_ids), "annotated_query_ids": n,
                            "annotated_background_ids": N, "tested_terms": len(results),
                            "significant_terms": int((frame.padj < .05).sum()) if not frame.empty else 0})
combined = pd.concat(all_results, ignore_index=True)
combined.to_csv(OUT / "over_representation_all_tests.csv", index=False)
combined[combined.padj < .05].to_csv(OUT / "over_representation_significant.csv", index=False)
shown = combined[combined.padj < .05].groupby(["group", "collection"], sort=False).head(10)
shown.to_csv(OUT / "figure_s3_displayed_terms.csv", index=False)
pd.DataFrame(diagnostics).to_csv(OUT / "background_and_mapping_summary.csv", index=False)
(OUT / "kegg_resource_manifest.json").write_text(json.dumps(resources, indent=2) + "\n")
(OUT / "kegg_release_information.txt").write_bytes((args.kegg_cache / "kegg_release.txt").read_bytes())
print(pd.DataFrame(diagnostics).to_string(index=False))
