"""v7.16 S-12 qualification census (D2 §9; D4 workplan S6).

Two steps, so the plan is frozen before anything runs:

    python validation/v716_s12_census.py freeze --rule B   # writes qualification/v716/census-plan.json
    python validation/v716_s12_census.py run    --rule B   # runs every request of the frozen plan

``run`` refuses to start unless the frozen commit is an ancestor of HEAD with
no change to the tracked code since, the files match their frozen hashes and
the stack matches. Each request runs in its own process
through the public entry point ``generate_lta_cation_structure_from_specs``
with the modelled f_V, K_max = 20 and a wall-clock cap of 90 min (a timeout
counts as a failure). All requests run in the frozen order; the total budget is
12 h, after which the remaining requests are recorded as unrun (not failed).

Decision rule (D2 §9): rule B is retained iff at most 4 of the 90 fail.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "qualification" / "v716"
CODE_FILES = (
    "lta_cage_salt.py",
    "lta_cage_salt_build.py",
    "lta_cage_salt_data.py",
    "lta_cation_population_grouped.py",
    "lta_cation_density_grouped.py",
    "validation/v716_s12_census.py",
)
CAP_S = 90 * 60
BUDGET_S = 12 * 3600
WORKERS = 4
RINGS = {"Na": dict(n4="Na:2", n6="Na:16", n8="Na:6"), "K": dict(n4="K:2", n6="K:16", n8="K:6")}
SEEDS = (0, 1, 2, 3, 4)
ALPHA_ROWS = [  # per alpha-cage, both alpha-cages the same (D2 §9)
    ("n2", "Cl:2,Li:1,K:1"),
    ("n4", "Cl:4,Li:2,K:2"),
    ("n6", "Cl:6,Li:4,K:2"),
    ("n8", "Cl:8,Li:5,K:3"),
    ("n10", "Cl:10,Li:6,K:4"),
    ("n12", "Cl:12,Li:7,K:5"),
    ("Li5K4Cs1", "Cl:10,Li:5,K:4,Cs:1"),
    ("Li5K3Sr1", "Cl:10,Li:5,K:3,Sr:1"),
    ("Li6K4Ba1", "Cl:12,Li:6,K:4,Ba:1"),
]
BETA = "Cl:1,Li:1"


def requests():
    out = []
    for ring in ("Na", "K"):
        for name, alpha in ALPHA_ROWS:
            for seed in SEEDS:
                out.append({"id": f"{ring}/{name}/s{seed}", "ring": ring, "alpha": alpha, "beta": BETA, "seed": seed})
    return out


def _git(*args):
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()


def stack():
    import lta_cage_salt as S

    st = S.stack_versions()
    st["threads"] = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
    return st


def _plan_path(rule):
    return OUT / ("census-plan.json" if rule == "B" else f"census-plan-{rule}.json")


def freeze(rule):
    sys.path.insert(0, str(REPO))
    import lta_cage_salt as S

    if S.COORDINATION_RULE != rule:
        raise SystemExit(f"COORDINATION_RULE is {S.COORDINATION_RULE!r}, not {rule!r}")
    dirty = _git("status", "--porcelain", "--", *CODE_FILES)
    if dirty:
        raise SystemExit(f"tracked code is modified; commit before freezing:\n{dirty}")
    plan = {
        "kind": "v7.16 S-12 qualification census plan",
        "rule": rule,
        "commit": _git("rev-parse", "HEAD"),
        "code_sha256": {f: hashlib.sha256((REPO / f).read_bytes()).hexdigest() for f in CODE_FILES},
        "frozen_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "entry_point": "lta_cation_population_grouped.generate_lta_cation_structure_from_specs",
        "volume_factor": "modelled f_V (None)",
        "k_max": 20,
        "per_request_cap_s": CAP_S,
        "budget_s": BUDGET_S,
        "workers": WORKERS,
        "thread_env": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
        "decision_rule": "retain the rule iff at most 4 of the 90 requests fail (D2 §9)",
        "d2_parameters": S.parameters_snapshot(20),
        "stack": stack(),
        "requests": requests(),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    _plan_path(rule).write_text(json.dumps(plan, indent=1, sort_keys=True) + "\n")
    print(f"froze {len(plan['requests'])} requests at {plan['commit']}")


WORKER = r"""
import json, sys, time
sys.path.insert(0, {repo!r})
import lta_cation_population_grouped as P
import lta_cage_salt as S
req = json.loads(sys.argv[1])
t = time.time()
out = {{"id": req["id"]}}
try:
    r = P.generate_lta_cation_structure_from_specs(
        {cif!r}, seed=req["seed"], salt_max_starts=20,
        cage_salt=S.CageSaltRequest.from_specs(alpha=req["alpha"], beta=req["beta"]), **req["fw"])
    rec = r.cage_salt_record
    out.update(passed=True, start_used=rec.verification["start_used"], n_starts=len(rec.starts),
               relocations=len(rec.relocations),
               min_slack=rec.verification["detail"]["V2"]["min_slack"],
               floors={{c: rec.floors[c]["floor_by_species"] for c in ("alpha-0", "alpha-1")}},
               n_melt=rec.floors["alpha-0"]["n_melt"])
except S.CageSaltConstructionError as exc:
    rec = exc.record
    out.update(passed=False, reason=str(exc), n_starts=len(rec.starts),
               last_failing=[r["failing"] for s in rec.starts for r in s["rounds"]][-1:] if rec.starts else None)
out["seconds"] = round(time.time() - t, 2)
print("RESULT " + json.dumps(out, sort_keys=True))
"""


def _run_one(req):
    payload = dict(req)
    payload["fw"] = RINGS[req["ring"]]
    code = WORKER.format(repo=str(REPO), cif=str(REPO / "tests" / "real_lta_primitive.cif"))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    t = time.time()
    try:
        proc = subprocess.run([sys.executable, "-c", code, json.dumps(payload)], cwd=REPO, env=env,
                              capture_output=True, text=True, timeout=CAP_S)
    except subprocess.TimeoutExpired:
        return {"id": req["id"], "passed": False, "reason": f"wall-clock cap {CAP_S} s", "seconds": CAP_S}
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    return {"id": req["id"], "passed": False, "reason": "crash: " + proc.stderr[-2000:],
            "seconds": round(time.time() - t, 2)}


def wilson(k, n, z=1.96):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 4), round(c + h, 4)]


def run(rule):
    sys.path.insert(0, str(REPO))
    plan = json.loads(_plan_path(rule).read_text())
    ancestor = subprocess.run(["git", "merge-base", "--is-ancestor", plan["commit"], "HEAD"], cwd=REPO)
    changed = subprocess.run(["git", "diff", "--quiet", plan["commit"], "HEAD", "--", *CODE_FILES], cwd=REPO)
    if ancestor.returncode != 0 or changed.returncode != 0:
        raise SystemExit(f"code changed since the frozen commit {plan['commit']}")
    for f, digest in plan["code_sha256"].items():
        if hashlib.sha256((REPO / f).read_bytes()).hexdigest() != digest:
            raise SystemExit(f"{f} differs from the frozen plan")
    now = stack()
    for key in ("python", "numpy", "scipy", "pymatgen", "platform"):
        if now.get(key) != plan["stack"].get(key):
            raise SystemExit(f"stack {key} differs from the frozen plan")
    log = OUT / f"census-record-{rule}.jsonl"
    log.write_text("")
    start = time.time()
    results = {}
    order = plan["requests"]

    def task(req):
        if time.time() - start > BUDGET_S:
            return {"id": req["id"], "unrun": True}
        res = _run_one(req)
        with open(log, "a") as fh:
            fh.write(json.dumps(res, sort_keys=True) + "\n")
        print(f"{res['id']}: {'PASS' if res.get('passed') else 'FAIL'} ({res.get('seconds')} s)", flush=True)
        return res

    with ThreadPoolExecutor(max_workers=plan["workers"]) as pool:
        for req, res in zip(order, pool.map(task, order)):
            results[req["id"]] = res
    ran = [r for r in results.values() if not r.get("unrun")]
    passed = sum(1 for r in ran if r.get("passed"))
    failed = len(ran) - passed
    unrun = len(order) - len(ran)
    if unrun == 0:
        decision = "retain" if failed <= 4 else "reject"
    elif failed >= 5:
        decision = "reject"
    elif passed >= 86:
        decision = "retain"
    else:
        decision = "undetermined"
    summary = {
        "rule": rule, "plan": str(_plan_path(rule).relative_to(REPO)), "commit": plan["commit"],
        "ran": len(ran), "passed": passed, "failed": failed, "unrun": unrun,
        "pass_rate": round(passed / len(ran), 4) if ran else None, "wilson_95": wilson(passed, len(ran)),
        "decision": decision, "wall_s": round(time.time() - start, 1),
        "failures": sorted(r["id"] for r in ran if not r.get("passed")),
        "by_group": {},
    }
    for ring in RINGS:
        for name, _ in ALPHA_ROWS:
            ids = [f"{ring}/{name}/s{s}" for s in SEEDS]
            rs = [results[i] for i in ids if not results[i].get("unrun")]
            summary["by_group"][f"{ring}/{name}"] = {
                "passed": sum(1 for r in rs if r.get("passed")), "ran": len(rs),
                "max_seconds": max((r.get("seconds", 0) for r in rs), default=None),
                "starts_used": [r.get("start_used") for r in rs],
            }
    (OUT / f"census-record-{rule}.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("ran", "passed", "failed", "unrun", "wilson_95", "decision")}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("freeze", "run"))
    ap.add_argument("--rule", choices=("B", "A"), default="B")
    args = ap.parse_args()
    (freeze if args.step == "freeze" else run)(args.rule)


if __name__ == "__main__":
    main()
