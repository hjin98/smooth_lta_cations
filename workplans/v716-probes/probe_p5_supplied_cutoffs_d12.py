"""Probe 5 (not production, not acceptance evidence): cage-salt feasibility with the
stakeholder-supplied pair table (doc/evidence/molten_salt_pair_distances.csv, rev 2.1)
and the D-12 same-family ring-cation relocation policy.

Hard minima: r_cut column of the supplied table.
  cation-Cl            : (M, Cl) row
  cation-cation        : (M, M') "Cl-bridged" row (no free oxide in scope)
  Cl-Cl                : "around M" rows, weighted by the cage-phase cation mole fraction
  cage ion - framework : no supplied value (pair file section 6.4). Two agent proxies:
      F1 = 0.85 x Shannon sum (P4 continuity; Si 0.26, Al 0.39, O 1.40)
      F2 = free-oxide rows (M-O r_cut; Cl-O "around M" weighted) for O, F1 for Si/Al
Ring-ring pairs are not enforced (existing template geometry); they are reported.

D-12: a ring cation that conflicts with the cage salt may move only to an empty site of
its own ring family; O1a (other positional variant of the same physical site) is a switch.

Coordination (reported only, no CN shaping): first-shell count inside r_min1 of the
opposite-charge pair, counting ring cations as neighbours of cage Cl (O4c proposed).

Run: python probe_p5_supplied_cutoffs_d12.py [grid]   (grid = main | sens | dopant | cn | cnd | cnm | all)
"""
import csv, sys, json, itertools, time, warnings
from pathlib import Path
from multiprocessing import Pool
import numpy as np
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
import lta_cation_population_grouped as P  # noqa: E402

PAIR_CSV = REPO / "doc/evidence/molten_salt_pair_distances.csv"
CATIONS = ("Li", "Na", "K", "Rb", "Cs", "Sr", "Ba")
CHARGE = dict(Li=1, Na=1, K=1, Rb=1, Cs=1, Sr=2, Ba=2)
SHANNON = dict(Li=0.76, Na=1.02, K=1.38, Rb=1.52, Cs=1.67, Sr=1.18, Ba=1.35,
               Cl=1.81, O=1.40, Si=0.26, Al=0.39)
EUTECTIC_XLI = 0.592  # [R1] x_KCl/x_LiCl = 0.69


def load_table():
    t = {}
    with open(PAIR_CSV) as fh:
        for row in csv.DictReader(fh):
            a, b, env = row["ion_i"], row["ion_j"], row["environment"]
            rec = dict(cut=float(row["r_cut_A"]), peak=float(row["r_peak_A"]),
                       min1=float(row["r_min1_A"]))
            t[(a, b, env)] = rec
            t[(b, a, env)] = rec
    return t


TAB = load_table()


def r_cut_pair(a, b, xcage):
    """Hard minimum for two salt-phase/ring ions; KeyError = missing pair (fail closed)."""
    if a == "Cl" and b == "Cl":
        return sum(x * TAB[("Cl", "Cl", f"around {m}")]["cut"] for m, x in xcage.items())
    if "Cl" in (a, b):
        return TAB[(a, b, "")]["cut"]
    return TAB[(a, b, "Cl-bridged")]["cut"]


def r_cut_framework(a, fw, xcage, rule):
    if rule == "F2" and fw == "O":
        if a == "Cl":
            return sum(x * TAB[("Cl", "O", f"around {m}")]["cut"] for m, x in xcage.items())
        return TAB[(a, "O", "")]["cut"]
    return 0.85 * (SHANNON[a] + SHANNON[fw])


def r_min1(a, b):
    return TAB[(a, b, "")]["min1"]


def mic_factory(M):
    Minv = np.linalg.inv(M)

    def mic(d):
        fr = d @ Minv
        fr -= np.round(fr)
        return fr @ M
    return mic


def build_case(case):
    ring_el = case["ring"]
    r = P.generate_lta_cation_structure_from_specs(
        str(REPO / "tests/real_lta_primitive.cif"), n4=f"{ring_el}:2", n6=f"{ring_el}:16",
        n8=f"{ring_el}:6", seed=case["seed"], volume_factor=case.get("vf", 1.0))
    st = r.structure
    L = st.lattice
    fw = [(x.specie.symbol, x.coords) for x in st if x.specie.symbol in ("Si", "Al", "O")]
    ring = [ion for ion in r.placed_ions]
    return r, L, fw, ring


def run_case(case):
    t0 = time.time()
    rng = np.random.default_rng(1000 + case["seed"] * 17 + case["n"])
    r, L, fw, ring = build_case(case)
    M = L.matrix
    mic = mic_factory(M)
    fw_el = [e for e, _ in fw]
    fw_xyz = np.array([c for _, c in fw])
    a_cent = [L.get_cartesian_coords(c) for c in ((.25, .25, .25), (.75, .75, .75))]
    b_cent = [L.get_cartesian_coords(c) for c in ((0, 0, 0), (.5, .5, .5))]
    all_cent = a_cent + b_cent

    # ---- inventory: per alpha cage n Cl + cations; beta per mode
    el, home, mobile_ring_charge = [], [], 0
    comp = case.get("alpha_comp")  # dict cation -> count per alpha
    n = case["n"]
    if comp is None:
        nli = round(EUTECTIC_XLI * n)
        comp = {"Li": nli, "K": n - nli}
    assert sum(CHARGE[m] * k for m, k in comp.items()) == n, comp
    for ci in range(2):
        for m, k in comp.items():
            el += [m] * k
            home += [ci] * k
        el += ["Cl"] * n
        home += [ci] * n
    ring_excess = []
    if case["beta"] in ("B1a", "B1b"):
        for bi in range(2):
            el.append("Cl")
            home.append(2 + bi)
            if case["beta"] == "B1a":
                el.append("Li")
                home.append(2 + bi)
            else:
                ring_excess.append("Li")
    home = np.array(home)
    cage_cations = [e for e in el if e != "Cl"]
    xcage = {}
    for m in set(cage_cations):
        xcage[m] = cage_cations.count(m) / len(cage_cations)
    # ---- ring state: occupied (site_id, variant) per ring ion; add ring-excess on empty 4R
    sites = {s.site_id: s for ring_sites in r.physical_sites_by_ring.values() for s in ring_sites}
    ring_state = [dict(el=ion.element, ring=ion.ring_class, site=ion.site_id,
                       var=ion.variant_id, frac=np.array(ion.frac)) for ion in ring]
    occupied = {s["site"] for s in ring_state}
    empty4 = [sid for sid, s in sites.items() if s.ring_class == "4" and sid not in occupied]
    for e in ring_excess:
        sid = empty4.pop(int(rng.integers(len(empty4))))
        ring_state.append(dict(el=e, ring="4", site=sid, var=0, frac=np.array(sites[sid].variants[0].frac)))
        occupied.add(sid)
    try:
        Dmm = np.array([[r_cut_pair(a, b, xcage) for b in el] for a in el])
        Dmf = np.array([[r_cut_framework(a, f, xcage, case["F"]) for f in fw_el] for a in el])
        ring_els = sorted({s["el"] for s in ring_state})
        Dmr_by = {re: np.array([r_cut_pair(a, re, xcage) for a in el]) for re in ring_els}
    except KeyError as exc:
        return dict(case=case, status="missing_pair", missing=str(exc), secs=time.time() - t0)

    def ring_xyz():
        return np.array([L.get_cartesian_coords(s["frac"]) for s in ring_state])

    def Dmr():
        return np.stack([Dmr_by[s["el"]] for s in ring_state], axis=1)

    rconf = case.get("rconf", 5.7)

    def penalty(x, R, DR):
        X = x.reshape(-1, 3)
        g = np.zeros_like(X)
        e = 0.0
        d = mic(X[:, None] - X[None])
        nrm = np.linalg.norm(d, axis=-1) + 1e-12
        v = np.maximum(Dmm - nrm, 0)
        np.fill_diagonal(v, 0)
        e += 0.5 * (v ** 2).sum()
        g -= 2 * ((v / nrm)[..., None] * d).sum(1)
        for Y, DD in ((fw_xyz, Dmf), (R, DR)):
            d2 = mic(X[:, None] - Y[None])
            n2 = np.linalg.norm(d2, axis=-1) + 1e-12
            v2 = np.maximum(DD - n2, 0)
            e += (v2 ** 2).sum()
            g -= 2 * ((v2 / n2)[..., None] * d2).sum(1)
        C = np.array([all_cent[h] for h in home])
        dc = mic(X - C)
        nc = np.linalg.norm(dc, axis=1) + 1e-12
        if case["C"] == "C1":
            lim = np.where(home < 2, rconf, 1.0)  # beta guests within 1.0 A of beta centre
            vc = np.maximum(nc - lim, 0)
            e += (vc ** 2).sum()
            g += 2 * (vc / nc)[:, None] * dc
        else:  # C2: nearest-centre partition (own centre closer than any other)
            for oc in all_cent:
                do = mic(X - oc)
                no = np.linalg.norm(do, axis=1) + 1e-12
                own = np.array([not np.allclose(oc, all_cent[h]) for h in home])
                vv = np.where(own, np.maximum(nc - no, 0), 0)
                e += (vv ** 2).sum()
                g += 2 * vv[:, None] * (dc / nc[:, None] - do / no[:, None])
        return e, g.ravel()

    def exact_slack(X, R, DR):
        d = np.linalg.norm(mic(X[:, None] - X[None]), axis=-1)
        np.fill_diagonal(d, np.inf)
        s_mm = (d - Dmm).min()
        s_mf = (np.linalg.norm(mic(X[:, None] - fw_xyz[None]), axis=-1) - Dmf).min()
        s_mr = (np.linalg.norm(mic(X[:, None] - R[None]), axis=-1) - DR).min()
        return min(s_mm, s_mf, s_mr), dict(mm=float(s_mm), mf=float(s_mf), mr=float(s_mr))

    def relocate(X):
        """D-12: move conflicting ring cations to the least-conflicting empty same-family slot."""
        moved = 0
        R = ring_xyz()
        DR = Dmr()
        viol = (DR - np.linalg.norm(mic(X[:, None] - R[None]), axis=-1)).max(axis=0)
        for k in np.argsort(-viol):
            if viol[k] <= 1e-3:
                break
            s = ring_state[k]
            cands = []
            for sid, site in sites.items():
                if site.ring_class != s["ring"]:
                    continue
                if sid != s["site"] and sid in occupied:
                    continue
                if sid == s["site"] and not case["O1a"]:
                    continue
                for v in site.variants:
                    if sid == s["site"] and v.variant_id == s["var"]:
                        continue
                    cands.append((sid, v))
            if not cands:
                continue
            col = Dmr_by[s["el"]]
            best, bestv = None, viol[k]
            for sid, v in cands:
                y = L.get_cartesian_coords(v.frac)
                w = (col - np.linalg.norm(mic(X - y), axis=-1)).max()
                if w < bestv - 1e-6:
                    best, bestv = (sid, v), w
            if best is not None:
                occupied.discard(s["site"])
                occupied.add(best[0])
                s.update(site=best[0], var=best[1].variant_id, frac=np.array(best[1].frac))
                moved += 1
        return moved

    starts = case.get("starts", 8)
    ok, best_slack, best_X, total_moves = 0, -np.inf, None, 0
    base_state = [dict(s) for s in ring_state]
    base_occ = set(occupied)
    for t in range(starts):
        ring_state[:] = [dict(s) for s in base_state]
        occupied.clear()
        occupied.update(base_occ)
        C = np.array([all_cent[h] for h in home])
        u = rng.normal(size=C.shape)
        u /= np.linalg.norm(u, axis=1)[:, None]
        rad = np.where(home < 2, rconf, 1.0) * rng.random(len(C)) ** (1 / 3)
        x = (C + u * rad[:, None]).ravel()
        rounds = 4 if case["reloc"] else 1
        for rd in range(rounds):
            R, DR = ring_xyz(), Dmr()
            res = minimize(penalty, x, args=(R, DR), jac=True, method="L-BFGS-B",
                           options=dict(maxiter=4000))
            x = res.x
            X = x.reshape(-1, 3)
            sl, parts = exact_slack(X, R, DR)
            if sl >= -1e-3 or not case["reloc"] or rd == rounds - 1:
                break
            mv = relocate(X)
            total_moves += mv
            if mv == 0:
                break
        if sl >= -1e-3:
            ok += 1
        if sl > best_slack:
            best_slack, best_X, best_parts = sl, X.copy(), parts
            best_ring = [dict(s) for s in ring_state]
    ring_state[:] = best_ring
    R = ring_xyz()
    Xc = best_X
    # ---- optional smooth CN shaping (P6): cation-centred M-Cl targets only (D-13)
    tgt = case.get("cn_target")
    if tgt and best_slack >= -1e-3:
        DR = Dmr()
        w = 0.12
        cat = np.array([i for i, e in enumerate(el) if e != "Cl" and home[i] < 2])
        cls = np.array([i for i, e in enumerate(el) if e == "Cl"])
        T = np.array([tgt[el[i]] for i in cat])
        RM = np.array([[r_min1(el[i], "Cl") for _ in cls] for i in cat])

        def shaped(x, W):
            e, g = penalty(x, R, DR)
            e, g = W * e, W * g
            X = x.reshape(-1, 3)
            d = mic(X[cat][:, None] - X[cls][None])
            nr = np.linalg.norm(d, axis=-1) + 1e-12
            z = np.clip((nr - RM) / w, -50, 50)
            s = 1 / (1 + np.exp(z))
            c = s.sum(1)
            res = c - T
            e += (res ** 2).sum()
            ds = -s * (1 - s) / w  # d s / d r
            coef = (2 * res[:, None] * ds / nr)[..., None] * d  # d/dX_cat
            G = g.reshape(-1, 3)
            np.add.at(G, cat, coef.sum(1))
            np.add.at(G, cls, -coef.sum(0))
            return e, G.ravel()
        x = Xc.ravel()
        for W in (10.0, 100.0, 1000.0):
            x = minimize(shaped, x, args=(W,), jac=True, method="L-BFGS-B", options=dict(maxiter=6000)).x
        x = minimize(penalty, x, args=(R, DR), jac=True, method="L-BFGS-B", options=dict(maxiter=4000)).x
        Xs = x.reshape(-1, 3)
        sl, parts = exact_slack(Xs, R, DR)
        best_slack, best_parts = sl, parts
        Xc = Xs
    # ---- coordination of the reported configuration
    cn = {}
    cl_idx = [i for i, e in enumerate(el) if e == "Cl" and home[i] < 2]
    cat_idx = [i for i, e in enumerate(el) if e != "Cl" and home[i] < 2]
    for i in cat_idx:
        m = el[i]
        d = np.linalg.norm(mic(Xc[[j for j, e in enumerate(el) if e == "Cl"]] - Xc[i]), axis=-1)
        cn.setdefault(f"{m}-Cl", []).append(int((d < r_min1(m, "Cl")).sum()))
    for i in cl_idx:
        k = 0
        for j, e in enumerate(el):
            if e != "Cl":
                k += np.linalg.norm(mic(Xc[j] - Xc[i])) < r_min1(e, "Cl")
        kr = sum(np.linalg.norm(mic(R[q] - Xc[i])) < r_min1(s["el"], "Cl") for q, s in enumerate(ring_state))
        cn.setdefault("Cl-M(cage)", []).append(int(k))
        cn.setdefault("Cl-M(cage+ring)", []).append(int(k + kr))
    cn_summary = {k: dict(mean=round(float(np.mean(v)), 2), min=int(min(v)), max=int(max(v)),
                          hist=np.bincount(v).tolist()) for k, v in cn.items()}
    # ring-ring like-pair check vs Cl-bridged r_cut (not enforced)
    dr = np.linalg.norm(mic(R[:, None] - R[None]), axis=-1)
    np.fill_diagonal(dr, np.inf)
    rr = min((dr[a, b] - r_cut_pair(ring_state[a]["el"], ring_state[b]["el"], xcage), a, b)
             for a in range(len(R)) for b in range(a + 1, len(R)))
    return dict(case=case, status="ok", feasible=ok, starts=starts, best_slack=round(float(best_slack), 3),
                slack_parts={k: round(v, 3) for k, v in best_parts.items()}, ring_moves=total_moves,
                cn=cn_summary, ring_ring_slack=round(float(rr[0]), 3),
                volume=round(float(L.volume), 1), secs=round(time.time() - t0, 1))


def grid(name):
    out = []
    if name in ("main", "all"):
        for ring, seed, n, reloc in itertools.product(("Na", "K", "Li"), (0, 1, 2), (8, 10, 11, 12), (False, True)):
            out.append(dict(grid="main", ring=ring, seed=seed, n=n, reloc=reloc, O1a=True,
                            F="F2", C="C1", beta="B1a"))
    if name in ("sens", "all"):
        for F, C, beta, n, vf in itertools.product(("F1", "F2"), ("C1", "C2"), ("B0", "B1a", "B1b"), (10, 12), (1.0, 0.97)):
            out.append(dict(grid="sens", ring="Na", seed=0, n=n, reloc=True, O1a=True, F=F, C=C, beta=beta, vf=vf))
        for n, o1a in itertools.product((10, 12), (False,)):
            out.append(dict(grid="sens", ring="Na", seed=0, n=n, reloc=True, O1a=o1a, F="F2", C="C1", beta="B1a"))
    if name in ("dopant", "all"):
        for seed, comp in itertools.product((0, 1, 2), ({"Li": 5, "K": 4, "Cs": 1}, {"Li": 5, "K": 3, "Sr": 1},
                                                         {"Li": 6, "K": 4, "Ba": 1})):
            out.append(dict(grid="dopant", ring="Na", seed=seed, n=sum(CHARGE[m] * k for m, k in comp.items()),
                            alpha_comp=comp, reloc=True, O1a=True, F="F2", C="C1", beta="B1a"))
        out.append(dict(grid="dopant", ring="Ca_check", seed=0, n=10, reloc=True, O1a=True, F="F2", C="C1",
                        beta="B1a", note="Ca ring (5A-like) -> expected missing pair"))
    if name in ("cn", "all"):
        targets = {"pure-1a": {"Li": 4.3, "K": 5.2, "Cs": 5.5},  # coordination file section 1a
                   "eutectic-1c": {"Li": 4.2, "K": 7.0, "Cs": 5.5}}  # section 1c (definition-dependent)
        for ring, seed, n, tname in itertools.product(("Na", "K"), (0, 1, 2), (8, 10, 12), targets):
            out.append(dict(grid="cn", ring=ring, seed=seed, n=n, reloc=True, O1a=True, F="F2", C="C1",
                            beta="B1a", cn_target=targets[tname], tname=tname, starts=4))
    if name in ("cnd", "all"):
        # P7: Option A targets scaled linearly by alpha Cl density relative to the 12/alpha cap
        base = {"Li": 4.3, "K": 5.2}
        for ring, seed, n in itertools.product(("Na", "K"), (0, 1, 2), (8, 10, 12)):
            tgt = {m: round(t * n / 12, 3) for m, t in base.items()}
            out.append(dict(grid="cnd", ring=ring, seed=seed, n=n, reloc=True, O1a=True, F="F2", C="C2",
                            beta="B1a", cn_target=tgt, tname="A-density-scaled", starts=4))
    if name in ("cnm", "all"):
        # P8: Option A targets scaled by alpha Cl density relative to the 800 K eutectic melt
        # (n_melt = 15.0 Cl per alpha free volume), floored at 1 (anchoring); low loadings included
        base = {"Li": 4.3, "K": 5.2}
        for ring, seed, n in itertools.product(("Na", "K"), (0, 1, 2), (2, 4, 6, 8, 12)):
            tgt = {m: round(max(1.0, t * n / 15.0), 3) for m, t in base.items()}
            out.append(dict(grid="cnm", ring=ring, seed=seed, n=n, reloc=True, O1a=True, F="F2", C="C2",
                            beta="B1a", cn_target=tgt, tname="A-melt800-scaled", starts=4))
    return out


def run_safe(case):
    try:
        if case["ring"] == "Ca_check":
            r_cut_pair("Ca", "Cl", {"Li": 1.0})
        return run_case(case)
    except KeyError as exc:
        return dict(case=case, status="missing_pair", missing=str(exc))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    cases = grid(which)
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(f"p5_{which}.jsonl")
    with Pool(4) as pool, open(out, "w") as fh:
        for res in pool.imap_unordered(run_safe, cases):
            fh.write(json.dumps(res) + "\n")
            fh.flush()
            c = res["case"]
            print(c["grid"], c["ring"], c["seed"], c["n"], c.get("F"), c.get("C"), c.get("beta"), c.get("vf", 1.0),
                  "reloc" if c.get("reloc") else "fixed", res["status"], res.get("feasible"), res.get("best_slack"),
                  res.get("ring_moves"), flush=True)
