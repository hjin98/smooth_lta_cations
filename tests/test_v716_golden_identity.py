"""v7.16 S4: dry and legacy outputs stay byte-identical to the captured v7.15.0 goldens.

Each case of ``tests/golden/v715/manifest.json`` is re-run with its recorded argv in
a scratch copy of the repository layout (so relative paths, and hence stdout, are
unchanged), and every captured file is compared by SHA-256.

- Dry cases are compared with no selector and under ``--occlusion-model sparse-v715``.
- Sparse cases are compared under the selector (their only allowed migration).
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "tests" / "golden" / "v715"
MANIFEST = json.loads((GOLDEN / "manifest.json").read_text())
CASES = sorted(k for k in MANIFEST if not k.startswith("_"))
SELECTOR = ["--occlusion-model", "sparse-v715"]


def _variants():
    for case in CASES:
        if "sparse" in case:
            yield case, True
        else:
            yield case, False
            yield case, True


def _run_case(case: str, with_selector: bool, tmp_path: Path) -> dict:
    entry = MANIFEST[case]
    argv = list(entry["argv"])
    work = tmp_path / "repo"
    (work / "tests").mkdir(parents=True)
    shutil.copy(ROOT / "tests" / "real_lta_primitive.cif", work / "tests" / "real_lta_primitive.cif")
    out_dir = work / "tests" / "golden" / "v715" / case
    out_dir.mkdir(parents=True)
    script = ROOT / argv[0]
    cmd = [sys.executable, str(script), *argv[1:]]
    if with_selector:
        cmd += SELECTOR
    proc = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    (out_dir / "stdout.txt").write_text(proc.stdout)
    # The maximum-symmetry report prints the absolute input path; the goldens were
    # captured in the repository root, so map the scratch root back to it.
    for name in ("stdout.txt", "report.txt"):
        path = out_dir / name
        if path.exists():
            path.write_text(path.read_text().replace(str(work), str(ROOT)))
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out_dir.iterdir())}


@pytest.mark.parametrize("case,with_selector", list(_variants()))
def test_v715_golden_byte_identity(case: str, with_selector: bool, tmp_path: Path):
    assert _run_case(case, with_selector, tmp_path) == MANIFEST[case]["files"]
