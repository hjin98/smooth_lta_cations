"""E1: exploratory CLI campaign for v7.16 cage-salt (agent probe; not acceptance evidence).

Asserter: AI agent (Claude Code session https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi).
Runs the public CLIs on requests outside the S-12 grid, then checks every written
structure with the test suite's independent verifier plus physical sanity checks
that V1-V6 do not make (salt-T distances, salt position inside its cage, CN versus
the density-reduced target, Cl-Cl spacing, beta content).

    python workplans/v716-probes/explore_e1_cli_campaign.py OUTDIR
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))
CIF = str(REPO / "tests" / "real_lta_primitive.cif")
POP = str(REPO / "lta_cation_population_grouped.py")
SMOOTH = str(REPO / "generate_smooth_lta_cations_grouped.py")
NA = ["--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6"]
NA8 = ["--n4", "Na:8", "--n6", "Na:16"]                 # 8R empty, 16 free 4R
LI = ["--n4", "Li:8", "--n6", "Li:16"]
K = ["--n4", "K:2", "--n6", "K:16", "--n8", "K:6"]
MIX = ["--n4", "Na:2", "--n6", "Na:8,K:8", "--n8", "K:6"]
SR = ["--n4", "Na:2", "--n6", "Sr:8", "--n8", "Na:6"]
CA = ["--n4", "Na:2", "--n6", "Ca:8", "--n8", "Na:6"]


def split(n):  # eutectic Li/K split per alpha-cage, as the census does
    li = round(0.592 * n)
    return f"Cl:{n},Li:{li},K:{n - li}"


CASES = []  # (id, script, args, expected exit, output name)
for n in (1, 3, 5, 7, 9, 11):
    CASES.append((f"A/odd_n{n}", POP, NA + ["--alpha-salt", split(n), "--beta-salt", "Cl:1,Li:1", "--seed", "1"], 0, "out.cif"))
CASES += [
    ("B/asym_12_0", POP, NA + ["--alpha0-salt", "Cl:12,Li:7,K:5", "--alpha1-salt", "Cl:0"], 0, "out.cif"),
    ("B/asym_12_2", POP, NA + ["--alpha0-salt", "Cl:12,Li:7,K:5", "--alpha1-salt", "Cl:2,Li:1,K:1", "--seed", "2"], 0, "out.cif"),
    ("B/asym_beta", POP, NA + ["--alpha-salt", "Cl:6,Li:4,K:2", "--beta0-salt", "Cl:1,Li:1", "--seed", "3"], 0, "out.cif"),
]
for name, spec in (("LiCl", "Cl:12,Li:12"), ("NaCl", "Cl:12,Na:12"), ("KCl", "Cl:12,K:12"), ("RbCl", "Cl:12,Rb:12"),
                   ("CsCl", "Cl:12,Cs:12"), ("SrCl2", "Cl:12,Sr:6"), ("BaCl2", "Cl:12,Ba:6"),
                   ("LiKBa", "Cl:12,Li:4,K:4,Ba:2"), ("CsLi", "Cl:12,Cs:6,Li:6"), ("Sr_dilute", "Cl:8,Li:4,K:2,Sr:1")):
    CASES.append((f"C/{name}", POP, NA + ["--alpha-salt", spec, "--seed", "4"], 0, "out.cif"))
CASES += [
    ("D/beta_only", POP, NA + ["--beta-salt", "Cl:1,Li:1", "--seed", "5"], 0, "out.cif"),
    ("D/beta_cl_counterion", POP, NA8 + ["--beta-salt", "Cl:1", "--salt-counterion-n4", "Li:2", "--seed", "5"], 0, "out.cif"),
    ("E/cl_only_4R", POP, NA + ["--alpha-salt", "Cl:4", "--salt-counterion-n4", "Li:8", "--seed", "6"], 0, "out.cif"),
    ("E/cl_only_4R_dense", POP, NA + ["--alpha-salt", "Cl:10", "--salt-counterion-n4", "Li:20", "--seed", "6"], 0, "out.cif"),
    ("E/mixed_counter", POP, NA + ["--alpha-salt", "Cl:8,Li:4,K:2", "--salt-counterion-n4", "Li:2,Sr:1", "--seed", "7"], 0, "out.cif"),
    ("F/8R_counter", POP, NA8 + ["--alpha-salt", "Cl:8,Li:4,K:2", "--salt-counterion-n8", "K:4", "--seed", "8"], 0, "out.cif"),
    ("F/4R8R_counter", POP, NA8 + ["--alpha-salt", "Cl:10,Li:4,K:2", "--salt-counterion-n4", "Li:4",
                                   "--salt-counterion-n8", "K:4", "--seed", "8"], 0, "out.cif"),
    ("G/Li_ring", POP, LI + ["--alpha-salt", "Cl:12,Li:7,K:5", "--seed", "9"], 0, "out.cif"),
    ("G/K_ring", POP, K + ["--alpha-salt", "Cl:12,Li:7,K:5", "--seed", "9"], 0, "out.cif"),
    ("G/mixed_ring", POP, MIX + ["--alpha-salt", "Cl:10,Li:6,K:4", "--seed", "9"], 0, "out.cif"),
    ("G/Sr_ring", POP, SR + ["--alpha-salt", "Cl:10,Li:6,K:4", "--seed", "9"], 0, "out.cif"),
    ("G/Ca_ring_rejected", POP, CA + ["--alpha-salt", "Cl:4,Li:2,K:2"], 2, "out.cif"),
    ("H/poscar", POP, NA + ["--alpha-salt", "Cl:8,Li:5,K:3", "--seed", "10"], 0, "out.POSCAR"),
    ("H/lammps", POP, NA + ["--alpha-salt", "Cl:8,Li:5,K:3", "--seed", "10"], 0, "out.data"),
    ("H/vf_override_095", POP, NA + ["--alpha-salt", "Cl:12,Li:7,K:5", "--volume-factor", "0.95", "--seed", "11"], 0, "out.cif"),
    ("H/vf_override_100", POP, NA + ["--alpha-salt", "Cl:12,Li:7,K:5", "--volume-factor", "1.0", "--seed", "11"], 0, "out.cif"),
    ("H/vf_tight_090", POP, NA + ["--alpha-salt", "Cl:12,Li:7,K:5", "--volume-factor", "0.90", "--seed", "11"], None, "out.cif"),
    ("H/seed_none", POP, NA + ["--alpha-salt", "Cl:6,Li:4,K:2"], 0, "out.cif"),
    ("H/max_starts_1", POP, NA + ["--alpha-salt", "Cl:2,Li:1,K:1", "--salt-max-starts", "1", "--seed", "0"], None, "out.cif"),
    ("I/smooth_regen", SMOOTH, NA + ["--alpha-salt", "Cl:10,Li:6,K:4", "--beta-salt", "Cl:1,Li:1", "--seed", "12",
                                     "--max-attempts", "3", "--mode", "regenerate"], 0, "out.cif"),
    ("I/smooth_move", SMOOTH, MIX + ["--alpha-salt", "Cl:10,Li:6,K:4", "--seed", "13", "--max-attempts", "3"], 0, "out.cif"),
    ("I/smooth_counter", SMOOTH, NA8 + ["--alpha-salt", "Cl:10,Li:4,K:2", "--salt-counterion-n4", "Li:8", "--seed", "14",
                                        "--max-attempts", "3"], 0, "out.cif"),
    ("J/cap13", POP, NA + ["--alpha-salt", "Cl:13,Li:13"], 2, "out.cif"),
    ("J/charge_mismatch", POP, NA + ["--alpha-salt", "Cl:5,Li:2,K:2"], 2, "out.cif"),
    ("J/counter_no_salt", POP, NA + ["--salt-counterion-n4", "Li:2"], 2, "out.cif"),
    ("J/beta2", POP, NA + ["--beta-salt", "Cl:2,Li:2"], 2, "out.cif"),
    ("J/cation_no_cl", POP, NA + ["--alpha0-salt", "Li:1", "--alpha1-salt", "Cl:1"], 2, "out.cif"),
    ("J/Mg_salt", POP, NA + ["--alpha-salt", "Cl:2,Mg:1"], 2, "out.cif"),
    ("J/legacy_no_selector", POP, NA + ["--alpha-cl", "1"], 2, "out.cif"),
    ("J/salt_with_legacy", POP, NA + ["--alpha-salt", "Cl:2,Li:2", "--occlusion-model", "sparse-v715"], 2, "out.cif"),
    ("J/ratio", POP, ["--ratio4", "Na:0.1", "--n6", "Na:16", "--n8", "Na:6", "--alpha-salt", "Cl:2,Li:2"], 2, "out.cif"),
    ("J/bad_spec", POP, NA + ["--alpha-salt", "Cl=2"], 2, "out.cif"),
    ("J/zero_salt_dry", POP, NA + ["--alpha-salt", "Cl:0", "--seed", "0"], 0, "out.cif"),
    ("J/capacity", POP, NA + ["--alpha-salt", "Cl:12", "--salt-counterion-n8", "Li:24"], 2, "out.cif"),
    ("K/sparse_legacy", POP, NA + ["--occlusion-model", "sparse-v715", "--alpha-cl", "1", "--occlusion-n4", "Cs:1",
                                   "--seed", "1"], 0, "out.cif"),
]


def run_case(case, outdir):
    cid, script, args, expect, outname = case
    d = outdir / cid.replace("/", "__")
    d.mkdir(parents=True, exist_ok=True)
    out = d / outname
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    t = time.time()
    proc = subprocess.run([sys.executable, script, CIF, str(out), *args, "--report", str(d / "report.txt")],
                          cwd=d, capture_output=True, text=True, env=env, timeout=3600)
    row = {"id": cid, "exit": proc.returncode, "expected_exit": expect, "seconds": round(time.time() - t, 1),
           "stderr_tail": proc.stderr.strip().splitlines()[-1:] if proc.stderr.strip() else []}
    (d / "stdout.txt").write_text(proc.stdout)
    (d / "stderr.txt").write_text(proc.stderr)
    side = Path(str(out) + ".salt.json")
    if side.exists():
        rec = json.loads(side.read_text())
        row["status"] = rec["status"]
        row["message"] = rec.get("message")
        row["start_used"] = rec.get("verification", {}).get("start_used")
        row["n_starts"] = len(rec.get("starts", []))
        row["relocations"] = len(rec.get("relocations", []))
        row["rebalance_moves"] = len(rec.get("rebalance", {}).get("moves", []))
    if proc.returncode == 0 and side.exists() and out.exists():
        row.update(analyse(out, rec))
    return row


def analyse(path, rec):
    from test_v716_cage_salt import independent_verify, _read_structure, _align_to_sidecar

    out = {"independent": independent_verify(path)}
    st = _read_structure(path)
    rows = rec["atoms"]
    F = _align_to_sidecar(st, rows)
    M = np.asarray(st.lattice.matrix, float)
    centres = {"alpha-0": (.25, .25, .25), "alpha-1": (.75, .75, .75), "beta-0": (0, 0, 0), "beta-1": (.5, .5, .5)}

    def d(a, b):
        x = np.asarray(b, float) - np.asarray(a, float)
        x -= np.round(x)
        return float(np.linalg.norm(x @ M))

    salt = [i for i, r in enumerate(rows) if r["role"] == "cage_salt"]
    T = [i for i, r in enumerate(rows) if r["element"] in ("Si", "Al")]
    cl = [i for i in salt if rows[i]["element"] == "Cl"]
    out["min_salt_T"] = min((d(F[i], F[j]) for i in salt for j in T), default=None)
    out["min_cl_cl"] = min((d(F[i], F[j]) for a, i in enumerate(cl) for j in cl[a + 1:]), default=None)
    radial = {}
    com = {}
    for c, cen in centres.items():
        mine = [i for i in salt if rows[i]["cage_id"] == c]
        if not mine:
            continue
        rs = [d(cen, F[i]) for i in mine]
        radial[c] = [round(min(rs), 2), round(float(np.mean(rs)), 2), round(max(rs), 2)]
        vec = np.array([((F[i] - np.asarray(cen)) + 0.5) % 1.0 - 0.5 for i in mine]) @ M
        com[c] = round(float(np.linalg.norm(vec.mean(axis=0))), 2)
    out["radial_min_mean_max"] = radial
    out["salt_com_offset"] = com
    # CN of salt cations versus the density-reduced target C (record floors)
    per = rec["verification"]["detail"]["V3"]["per_ion"]
    cn = {}
    for ion_id, v in per.items():
        _, a, b, el, _ = ion_id.split("-")
        cn.setdefault(f"{a}-{b}/{el}", []).append(v["cn"])
    out["cn_mean"] = {k: round(float(np.mean(v)), 2) for k, v in sorted(cn.items())}
    out["C_target"] = {c: {m: round(x, 2) for m, x in rec["floors"][c]["c_by_species"].items()}
                       for c in ("alpha-0", "alpha-1") if c in rec["floors"]}
    out["floors"] = {c: rec["floors"][c]["floor_by_species"] for c in ("alpha-0", "alpha-1") if c in rec["floors"]}
    out["volume_factor"] = rec["volume"]["applied_factor"]
    out["q_cage_e"] = {c: round(rec["verification"]["detail"]["V4"]["q_exc"][c]
                                + rec["verification"]["detail"]["V4"]["q_fw"][c], 2) for c in centres}
    return out


def main():
    outdir = Path(sys.argv[1]).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    log = REPO / "workplans" / "v716-probes" / "results" / "e1_cli_campaign.jsonl"
    log.write_text("")
    with ThreadPoolExecutor(max_workers=4) as pool:
        for row in pool.map(lambda c: run_case(c, outdir), CASES):
            ok = row["expected_exit"] is None or row["exit"] == row["expected_exit"]
            indep = row.get("independent")
            if indep is not None and not all(indep.values()):
                ok = False
            row["as_expected"] = ok
            with open(log, "a") as fh:
                fh.write(json.dumps(row, sort_keys=True) + "\n")
            print(f"{'OK ' if ok else 'BAD'} {row['id']}: exit {row['exit']} ({row['seconds']} s)"
                  f" {row.get('status', '')} starts={row.get('start_used')}", flush=True)


if __name__ == "__main__":
    main()
