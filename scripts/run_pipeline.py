"""Run the downstream analysis with checked inputs and fixed reference resources."""
from pathlib import Path
import argparse
import os
import subprocess
import sys

from check_inputs import ROOT, verify_inputs

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rscript", default="Rscript", help="Path to the Rscript executable")
parser.add_argument("--kegg-cache", type=Path, default=ROOT / "data/cache/kegg")
parser.add_argument("--regenerate-go", action="store_true",
                    help="Regenerate the supplied GO snapshot with the recorded annotation packages")
parser.add_argument("--skip-figures", action="store_true")
args = parser.parse_args()
verify_inputs()
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", TZ="UTC")
stages = [
    [args.rscript, "scripts/check_environment.R"],
    [sys.executable, "scripts/analyze_bulk.py"],
    [sys.executable, "scripts/analyze_public_data.py"],
    [args.rscript, "scripts/analyze_single_cell.R"],
]
if args.regenerate_go:
    stages.append([args.rscript, "scripts/prepare_go_annotations.R"])
stages += [
    [sys.executable, "scripts/analyze_enrichment.py", "--kegg-cache", str(args.kegg_cache.resolve())],
    [args.rscript, "scripts/analyze_gsea.R"],
    [args.rscript, "scripts/verify_numerics.R"],
    [sys.executable, "scripts/validate_results.py"],
]
if not args.skip_figures:
    stages.append([sys.executable, "scripts/plot_figures.py"])
for command in stages:
    print("Running " + Path(command[1]).name, flush=True)
    subprocess.run(command, cwd=ROOT, env=env, check=True)
print("Pipeline completed.")
