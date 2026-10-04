"""Check exact source bytes before starting an analysis."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def verify_inputs():
    failures = []
    manifest = json.loads((ROOT / "config/source_data_manifest.json").read_text())
    for row in manifest:
        path = ROOT / row["path"]
        if not path.is_file():
            failures.append("Missing: " + row["path"])
            continue
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if path.stat().st_size != row["bytes"] or digest != row["sha256"]:
            failures.append("Different source bytes: " + row["path"])
    if failures:
        raise SystemExit("Input verification failed. No analysis started.\n" +
                         "\n".join(failures) + "\nSee docs/Data_sources.md.")
    print(f"All {len(manifest)} source files match the recorded checksums.")


if __name__ == "__main__":
    verify_inputs()
