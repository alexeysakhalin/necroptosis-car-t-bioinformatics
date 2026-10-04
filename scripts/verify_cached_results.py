"""Verify distributed result integrity and ORA arithmetic without restricted inputs.

This is an audit of saved results, not a new fit from the original source data.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
summary = json.loads((ROOT / "results/tables/validation/validation_summary.json").read_text())
for name, expected in summary["input_and_result_hashes"].items():
    path = ROOT / name
    with path.open("rb") as stream:
        actual = hashlib.file_digest(stream, "sha256").hexdigest()
    if actual != expected:
        raise ValueError("Saved analysis file changed: " + name)

ora = pd.read_csv(ROOT / "results/tables/enrichment/over_representation_all_tests.csv")
recomputed = stats.hypergeom.sf(ora["count"] - 1, ora.background_annotated_genes,
                               ora.background_term_genes, ora.query_annotated_genes)
np.testing.assert_allclose(ora.pvalue, recomputed, rtol=1e-10, atol=1e-300)
for _, frame in ora.groupby(["group", "collection"]):
    adjusted = stats.false_discovery_control(frame.pvalue.to_numpy(), method="bh")
    np.testing.assert_allclose(frame.padj, adjusted, rtol=1e-10, atol=1e-300)
print(f"Verified {len(summary['input_and_result_hashes'])} saved file hashes and "
      f"{len(ora)} over-representation tests, including multiple-testing correction.")
print("Source-data recomputation requires the exact inputs described in docs/Data_sources.md.")
