"""Capture v7.15.0 golden outputs for the v7.16 legacy/dry byte-identity tests.

Evidence capture on unchanged v7.15.0 code (S0a at 949f4fe; S0b extension). Every case
runs in its own directory under tests/golden/v715/<case>/; all produced files, stdout and
the --report text (where supported) are hashed into manifest.json together with the
numerical stack. Run twice and compare to confirm determinism.
"""
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "tests" / "golden" / "v715"
CIF = "tests/real_lta_primitive.cif"
POP = "lta_cation_population_grouped.py"
SMOOTH = "generate_smooth_lta_cations_grouped.py"
DRY = ["--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6"]
SPARSE = DRY + ["--alpha-cl", "2", "--beta-cl", "2", "--occlusion-n4", "Na:4"]

CASES = []  # (case, script, args, output filename, report?)
for seed in (0, 1, 2):
    for ext in ("cif", "POSCAR", "data"):
        CASES.append((f"pop_dry_s{seed}_{ext}", POP, DRY + ["--seed", str(seed)], f"out.{ext}", True))
for seed in (0, 1):
    for ext in ("cif", "POSCAR", "data"):
        CASES.append((f"pop_sparse_s{seed}_{ext}", POP, SPARSE + ["--seed", str(seed)], f"out.{ext}", True))
        CASES.append((f"smooth_dry_s{seed}_{ext}", SMOOTH, DRY + ["--max-attempts", "5", "--seed", str(seed)], f"out.{ext}", False))
        CASES.append((f"smooth_sparse_s{seed}_{ext}", SMOOTH, SPARSE + ["--max-attempts", "5", "--seed", str(seed)], f"out.{ext}", False))
    for mode in ("regenerate", "move"):
        CASES.append((f"smooth_dry_{mode}_s{seed}", SMOOTH, DRY + ["--mode", mode, "--max-attempts", "5", "--seed", str(seed)], "out.cif", False))
    CASES.append((f"pop_dry_vf098_s{seed}", POP, DRY + ["--volume-factor", "0.98", "--seed", str(seed)], "out.cif", True))
    CASES.append((f"smooth_dry_vf098_s{seed}", SMOOTH, DRY + ["--volume-factor", "0.98", "--max-attempts", "5", "--seed", str(seed)], "out.cif", False))
    CASES.append((f"pop_sparse_vf098_s{seed}", POP, SPARSE + ["--volume-factor", "0.98", "--seed", str(seed)], "out.cif", True))
CASES.append(("pop_unique_batch2_s0", POP, DRY + ["--unique-batch-size", "2", "--seed", "0"], "out.cif", True))
CASES.append(("pop_maxsym", POP, DRY + ["--population-mode", "maximum-symmetry"], "out.cif", True))
CASES.append(("smooth_maxsym", SMOOTH, DRY + ["--population-mode", "maximum-symmetry"], "out.cif", False))
CASES.append(("pop_dry_candidate_s0", POP, DRY + ["--write-candidate", "tests/golden/v715/pop_dry_candidate_s0/candidate.cif", "--seed", "0"], "out.cif", True))


def stack() -> dict:
    import numpy, pymatgen.core
    try:
        import scipy
        sv = scipy.__version__
    except ImportError:
        sv = None
    return {"python": platform.python_version(), "numpy": numpy.__version__, "scipy": sv,
            "pymatgen": pymatgen.core.__version__, "platform": platform.platform()}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {"_stack": stack(), "_base_commit": "769d999 (v7.15.0 code)"}
    for case, script, args, outname, report in CASES:
        d = OUT / case
        d.mkdir(parents=True, exist_ok=True)
        for old in d.iterdir():
            old.unlink()
        rel = d.relative_to(REPO)
        cmd = [sys.executable, script, CIF, str(rel / outname), *args]
        if report:
            cmd += ["--report", str(rel / "report.txt")]
        proc = subprocess.run(cmd, cwd=REPO, check=True, capture_output=True, text=True)
        (d / "stdout.txt").write_text(proc.stdout)
        manifest[case] = {
            "argv": cmd[1:],
            "files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(d.iterdir())},
        }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"captured {len(CASES)} cases")


if __name__ == "__main__":
    main()
