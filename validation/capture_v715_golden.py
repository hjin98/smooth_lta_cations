"""Capture v7.15.0 golden outputs for the v7.16 legacy/dry byte-identity tests.

Run once on the pre-v7.16 code (commit 769d999). Writes tests/golden/v715/ files
plus manifest.json (command line -> SHA-256 of every output file).
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "tests" / "golden" / "v715"
CIF = "tests/real_lta_primitive.cif"
POP = "lta_cation_population_grouped.py"
SMOOTH = "generate_smooth_lta_cations_grouped.py"
DRY = ["--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6"]
SPARSE = ["--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "2", "--beta-cl", "2",
          "--occlusion-n4", "Na:4"]

CASES = []
for seed in (0, 1, 2):
    for ext in ("cif", "POSCAR", "data"):
        CASES.append((f"pop_dry_s{seed}.{ext}", POP, DRY, seed))
for seed in (0, 1):
    for ext in ("cif", "POSCAR", "data"):
        CASES.append((f"pop_sparse_s{seed}.{ext}", POP, SPARSE, seed))
    CASES.append((f"smooth_dry_s{seed}.cif", SMOOTH, DRY + ["--max-attempts", "5"], seed))
    CASES.append((f"smooth_sparse_s{seed}.cif", SMOOTH, SPARSE + ["--max-attempts", "5"], seed))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for name, script, args, seed in CASES:
        target = OUT / name
        cmd = [sys.executable, script, CIF, str(target.relative_to(REPO)), *args, "--seed", str(seed)]
        subprocess.run(cmd, cwd=REPO, check=True, capture_output=True)
        manifest[name] = {"argv": cmd[1:], "sha256": hashlib.sha256(target.read_bytes()).hexdigest()}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"captured {len(manifest)} outputs")


if __name__ == "__main__":
    main()
