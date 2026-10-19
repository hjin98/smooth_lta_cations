"""v7.16 cage-salt construction numerics (D3 r3 "lta_cage_salt_build.py (numerics)").

Implements D2 r2 §6 (doc/specs/v7_16_cage_salt_numerical_method.md): random sequential
insertion, L-BFGS-B penalty continuation with smooth coordination counts, a final
dead-zone polish, charge-constrained D-12 ring relocation and first-success restarts.
Acceptance is never decided here: every candidate is passed to the caller's verifier,
which runs the exact D2 §7 checks.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize

import lta_cage_salt as S


@dataclass
class SaltIon:
    ion_id: str
    element: str
    cage_id: str


@dataclass
class BuildProblem:
    matrix: np.ndarray                      # final cell (rows = lattice vectors, Å)
    salt: List[SaltIon]
    framework_elements: List[str]
    framework_fracs: np.ndarray
    ring_ions: List[Any]                    # PlacedIon-like (element, frac, ring_class, site_id, variant_id, role)
    sites_by_ring: Mapping[str, Sequence[Any]]
    pairs: S.PairTable
    floors: Mapping[str, S.CageFloors]
    centres: Mapping[str, Sequence[float]]
    request: S.CageSaltRequest
    charges: Mapping[str, int]
    rule: str = S.COORDINATION_RULE


@dataclass
class BuildOutcome:
    success: bool
    fracs: Optional[np.ndarray]
    ring_ions: List[Any]
    start_index: Optional[int]
    starts: List[Dict[str, Any]] = field(default_factory=list)
    relocations: List[Dict[str, Any]] = field(default_factory=list)
    verification: Optional[S.Verification] = None
    stream_tags: List[Optional[str]] = field(default_factory=list)


def _h(u: np.ndarray) -> np.ndarray:
    return np.maximum(u, 0.0)


class _Model:
    """Precomputed target matrices and the analytic objective of D2 §6.1."""

    def __init__(self, prob: BuildProblem, ring_ions: Sequence[Any]):
        self.p = prob
        M = np.asarray(prob.matrix, float)
        self.M = M
        self.Minv = np.linalg.inv(M)
        el = [s.element for s in prob.salt]
        self.el = el
        n = len(el)
        pt = prob.pairs
        db = S.DELTA_B
        self.is_cl = np.array([e == "Cl" for e in el])
        self.R_ss = np.array([[pt.r_cut(a, b) + db for b in el] for a in el]) if n else np.zeros((0, 0))
        self.fw_cart = np.asarray(prob.framework_fracs, float) @ M
        fw_el = prob.framework_elements
        guard = S.GUARD_DISTANCE + db
        # D2 §6.1 keeps the O hard minimum and the 1.5 Å guard as separate terms.
        self.R_so = np.array([[pt.r_cut(a, f) + db if f == "O" else 0.0 for f in fw_el] for a in el]).reshape(
            n, len(fw_el))
        self.R_sg = np.full((n, len(fw_el)), guard)
        # D1 r8 S-13 / D2 r3: every salt ion >= T_min from framework Si/Al
        self.R_st = np.array([[S.T_MIN + db if f in S.T_ELEMENTS else 0.0 for f in fw_el] for _ in el]).reshape(
            n, len(fw_el))
        self.set_ring(ring_ions)
        # opposite-charge salt pairs and their counting radii
        rm = np.zeros((n, n))
        opp = np.zeros((n, n), bool)
        for i in range(n):
            for j in range(n):
                if self.is_cl[i] != self.is_cl[j]:
                    cat = el[j] if self.is_cl[i] else el[i]
                    rm[i, j] = pt.r_min1(cat)
                    opp[i, j] = True
        self.rm_ss, self.opp_ss = rm, opp
        names = list(S.CAGE_IDS)
        self.cen = np.array([prob.centres[c] for c in names], float) @ M
        self.home = np.array([names.index(s.cage_id) for s in prob.salt], int)
        self.lower = np.zeros(n)
        self.upper = np.full(n, np.inf)
        self.groups: List[Tuple[np.ndarray, float]] = []
        alpha_cation = np.array([(not self.is_cl[i]) and prob.salt[i].cage_id in prob.floors for i in range(n)])
        for i, s in enumerate(prob.salt):
            if self.is_cl[i]:
                self.lower[i] = 1.0
            elif alpha_cation[i]:
                fl = prob.floors[s.cage_id]
                self.lower[i] = fl.floor_by_species[s.element] if prob.rule == "B" else 1.0
                self.upper[i] = S.U_M[s.element]
            else:
                self.lower[i] = 1.0
        if prob.rule == "A":
            self._rule_a_groups(alpha_cation)

    def _rule_a_groups(self, alpha_cation: np.ndarray) -> None:
        p = self.p
        species = sorted({p.salt[i].element for i in range(len(p.salt)) if alpha_cation[i]})
        for m in species:
            by_cage = {c: [i for i in range(len(p.salt)) if alpha_cation[i] and p.salt[i].element == m
                           and p.salt[i].cage_id == c] for c in S.ALPHA_IDS}
            pooled: List[int] = []
            for c, idx in by_cage.items():
                if len(idx) >= 3:
                    self.groups.append((np.array(idx), S.A_RULE_A * p.floors[c].c_by_species[m] * len(idx)))
                elif idx:
                    pooled.extend(idx)
            if len(pooled) >= 3:
                req = sum(S.A_RULE_A * p.floors[p.salt[i].cage_id].c_by_species[m] for i in pooled)
                self.groups.append((np.array(pooled), req))

    def set_ring(self, ring_ions: Sequence[Any]) -> None:
        pt = self.p.pairs
        self.ring = list(ring_ions)
        self.ring_cart = np.array([np.asarray(r.frac, float) for r in ring_ions]).reshape(-1, 3) @ self.M
        self.R_sr = np.array([[pt.r_cut(a, r.element) + S.DELTA_B for r in ring_ions] for a in self.el]).reshape(
            len(self.el), len(ring_ions))
        self.rm_sr = np.array([[pt.r_min1(r.element) if a == "Cl" else 0.0 for r in ring_ions] for a in self.el]
                              ).reshape(len(self.el), len(ring_ions))
        self.opp_sr = np.array([[a == "Cl" for _ in ring_ions] for a in self.el], bool).reshape(
            len(self.el), len(ring_ions))

    def mic(self, v: np.ndarray) -> np.ndarray:
        f = v @ self.Minv
        f = f - np.round(f)
        return f @ self.M

    def energy(self, x: np.ndarray, W: float, Wdz: float) -> Tuple[float, np.ndarray]:
        X = x.reshape(-1, 3)
        n = len(X)
        G = np.zeros_like(X)
        E = 0.0
        # salt-salt hard minima
        V = self.mic(X[:, None, :] - X[None, :, :])
        D = np.linalg.norm(V, axis=-1)
        np.fill_diagonal(D, 1.0)
        viol = _h(self.R_ss - D)
        np.fill_diagonal(viol, 0.0)
        E += 0.5 * W * float((viol ** 2).sum())
        G += W * (-2.0 * (viol / D)[..., None] * V).sum(axis=1)
        # salt vs ring and framework
        for Y, R in ((self.ring_cart, self.R_sr), (self.fw_cart, self.R_so), (self.fw_cart, self.R_sg),
                     (self.fw_cart, self.R_st)):
            if len(Y) == 0:
                continue
            V2 = self.mic(X[:, None, :] - Y[None, :, :])
            D2 = np.linalg.norm(V2, axis=-1) + 1e-300
            v2 = _h(R - D2)
            E += W * float((v2 ** 2).sum())
            G += W * (-2.0 * (v2 / D2)[..., None] * V2).sum(axis=1)
        # membership margin
        Vc = self.mic(X[:, None, :] - self.cen[None, :, :])
        Dc = np.linalg.norm(Vc, axis=-1) + 1e-300
        Uc = Vc / Dc[..., None]
        idx = np.arange(n)
        dh = Dc[idx, self.home]
        Do = Dc.copy()
        Do[idx, self.home] = np.inf
        oth = np.argmin(Do, axis=1)
        margin = Do[idx, oth] - dh
        vm = _h(S.EPS_M - margin)
        E += W * float((vm ** 2).sum())
        G += W * (-2.0 * vm)[:, None] * (Uc[idx, oth] - Uc[idx, self.home])
        # dead zone around r_min1 (opposite-charge pairs), only in the polish stage
        if Wdz > 0:
            g = np.abs(D - self.rm_ss)
            vz = np.where(self.opp_ss, _h(S.DELTA_CN - g), 0.0)
            sgn = np.sign(D - self.rm_ss)
            E += Wdz * float((vz ** 2).sum())
            G += Wdz * ((4.0 * vz * (-sgn) / D)[..., None] * V).sum(axis=1)
            if len(self.ring_cart):
                V3 = self.mic(X[:, None, :] - self.ring_cart[None, :, :])
                D3 = np.linalg.norm(V3, axis=-1) + 1e-300
                g3 = np.abs(D3 - self.rm_sr)
                vz3 = np.where(self.opp_sr, _h(S.DELTA_CN - g3), 0.0)
                E += Wdz * float((vz3 ** 2).sum())
                G += Wdz * ((2.0 * vz3 * (-np.sign(D3 - self.rm_sr)) / D3)[..., None] * V3).sum(axis=1)
        # smooth counts
        z = np.clip((self.rm_ss - D) / S.W_S, -50, 50)
        sig = np.where(self.opp_ss, 1.0 / (1.0 + np.exp(-z)), 0.0)
        np.fill_diagonal(sig, 0.0)
        s = sig.sum(axis=1)
        dsig = -sig * (1.0 - sig) / S.W_S  # d sigma / d D
        if len(self.ring_cart):
            V4 = self.mic(X[:, None, :] - self.ring_cart[None, :, :])
            D4 = np.linalg.norm(V4, axis=-1) + 1e-300
            z4 = np.clip((self.rm_sr - D4) / S.W_S, -50, 50)
            sig4 = np.where(self.opp_sr, 1.0 / (1.0 + np.exp(-z4)), 0.0)
            s = s + sig4.sum(axis=1)
            dsig4 = -sig4 * (1.0 - sig4) / S.W_S
        lo = _h(self.lower + S.MU - s)
        up = np.where(np.isfinite(self.upper), _h(s - (self.upper - S.MU)), 0.0)
        E += float((lo ** 2).sum() + (up ** 2).sum())
        dEds = -2.0 * lo + 2.0 * up
        for members, req in self.groups:
            gv = max(0.0, req + S.MU - float(s[members].sum()))
            E += gv * gv
            dEds[members] += -2.0 * gv
        dEdD = (dEds[:, None] + dEds[None, :]) * dsig
        np.fill_diagonal(dEdD, 0.0)
        G += ((dEdD / D)[..., None] * V).sum(axis=1)
        if len(self.ring_cart):
            G += ((dEds[:, None] * dsig4 / D4)[..., None] * V4).sum(axis=1)
        return E, G.ravel()


def _insert(model: _Model, rng: np.random.Generator) -> np.ndarray:
    n = len(model.el)
    X = np.zeros((n, 3))
    placed: List[int] = []
    order = rng.permutation(n)
    for i in order:
        f = rng.random((S.M_INS, 3))
        C = f @ model.M
        tot = np.zeros(S.M_INS)
        Vc = model.mic(C[:, None, :] - model.cen[None, :, :])
        Dc = np.linalg.norm(Vc, axis=-1)
        dh = Dc[:, model.home[i]]
        Do = Dc.copy()
        Do[:, model.home[i]] = np.inf
        tot += _h(S.EPS_M - (Do.min(axis=1) - dh)) ** 2
        for Y, R in ((model.fw_cart, model.R_so[i]), (model.fw_cart, model.R_sg[i]), (model.fw_cart, model.R_st[i]),
                     (model.ring_cart, model.R_sr[i])):
            if len(Y):
                D = np.linalg.norm(model.mic(C[:, None, :] - Y[None, :, :]), axis=-1)
                tot += (_h(R[None, :] - D) ** 2).sum(axis=1)
        if placed:
            Y = X[placed]
            D = np.linalg.norm(model.mic(C[:, None, :] - Y[None, :, :]), axis=-1)
            tot += (_h(model.R_ss[i, placed][None, :] - D) ** 2).sum(axis=1)
        ok = np.flatnonzero(tot == 0.0)
        X[i] = C[ok[0]] if len(ok) else C[int(np.argmin(tot))]
        placed.append(i)
    return X


def _relax(model: _Model, X: np.ndarray, log: Dict[str, Any]) -> np.ndarray:
    x = X.ravel()
    stages = [(w, 0.0) for w in S.W_STAGES] + [(S.W_STAGES[-1], S.W_DZ_POLISH)]
    log["energies"] = []
    for W, Wdz in stages:
        res = minimize(model.energy, x, args=(W, Wdz), jac=True, method="L-BFGS-B",
                       options=dict(S.LBFGS_OPTIONS))
        x = res.x
        log["energies"].append({"W": W, "W_dz": Wdz, "E": float(res.fun), "nit": int(res.nit)})
    return x.reshape(-1, 3)


def _relocate(model: _Model, X: np.ndarray, prob: BuildProblem, ring: List[Any]) -> Tuple[List[Any], List[Dict]]:
    """D2 §6.4: resolving-only, charge-constrained moves within the ring family."""
    pt = prob.pairs
    tau = S.TAU_W
    viols = []
    for k, r in enumerate(ring):
        y = np.asarray(r.frac, float) @ model.M
        D = np.linalg.norm(model.mic(X - y[None, :]), axis=1)
        need = np.array([pt.r_cut(e, r.element) + tau for e in model.el])
        v = float((need - D).max()) if len(D) else -np.inf
        if v > 0:
            viols.append((-v, r.site_id, k))
    viols.sort()
    moves: List[Dict] = []
    ring = list(ring)
    for _, _, k in viols:
        r = ring[k]
        occupied = {q.site_id for q in ring}
        cands = []
        for site in prob.sites_by_ring.get(r.ring_class, []):
            if site.site_id != r.site_id and site.site_id in occupied:
                continue
            for vi, var in enumerate(site.variants):
                vid = int(getattr(var, "variant_id", vi))
                if site.site_id == r.site_id and vid == r.variant_id:
                    continue
                cands.append((site, var, vid))
        best = None
        rho = np.array([pt.r_cut(e, r.element) + S.DELTA_B for e in model.el])
        for site, var, vid in cands:
            y = np.asarray(var.frac, float) @ model.M
            D = np.linalg.norm(model.mic(X - y[None, :]), axis=1)
            slack = float((D - rho).min()) if len(D) else np.inf
            if slack < 0:
                continue
            trial = list(ring)
            trial[k] = replace(r, frac=np.asarray(var.frac, float), site_id=site.site_id, variant_id=vid,
                               site_label=getattr(var, "site_label", r.site_label))
            q_exc, q_fw, zmax = S.cage_charges_quarters(prob.request, trial, prob.charges)
            if any(S.charge_violations(q_exc, q_fw, zmax).values()):
                continue
            key = (-slack, site.site_id, vid)
            if best is None or key < best[0]:
                best = (key, trial[k])
        if best is not None:
            new = best[1]
            moves.append({"ion": f"R{r.site_id}", "to_ion_id": f"R{new.site_id}", "element": r.element, "role": r.role,
                          "from_site": r.site_id, "from_variant": r.variant_id,
                          "to_site": new.site_id, "to_variant": new.variant_id,
                          "min_slack_after": -best[0][0]})
            ring[k] = new
    return ring, moves


def build(prob: BuildProblem, verifier: Callable[[np.ndarray, List[Any]], S.Verification],
          seed: Optional[int], k_max: int = S.K_MAX_DEFAULT) -> BuildOutcome:
    """First-success restarts over per-start substreams R_s,k (D2 §6.5)."""
    outcome = BuildOutcome(False, None, list(prob.ring_ions), None)
    for k in range(1, k_max + 1):
        rng, tag = S.start_stream(seed, k)
        outcome.stream_tags.append(tag)
        ring = list(prob.ring_ions)
        model = _Model(prob, ring)
        log: Dict[str, Any] = {"start": k, "rounds": []}
        X = _insert(model, rng)
        start_moves: List[Dict] = []
        ver = None
        for rnd in range(S.R_RELOC + 1):
            rl: Dict[str, Any] = {}
            X = _relax(model, X, rl)
            fr = S.canonical_frac(X @ model.Minv)
            ver = verifier(fr, ring)
            rl["failing"] = list(ver.failing)
            rl["min_slack"] = ver.detail.get("V2", {}).get("min_slack", {})
            rl["min_membership_margin"] = ver.detail.get("V1", {}).get("min_membership_margin")
            rl["min_distance_to_r_min1"] = ver.detail.get("V3", {}).get("min_distance_to_r_min1")
            per_ion = ver.detail.get("V3", {}).get("per_ion", {})
            rl["cn_deficit"] = sum(max(0, v["floor"] - v["cn"]) for v in per_ion.values())
            rl["cn_excess"] = sum(max(0, v["cn"] - v["cap"]) for v in per_ion.values() if v["cap"] is not None)
            log["rounds"].append(rl)
            if ver.passed:
                break
            if "V2_hard_minima" not in ver.failing or rnd == S.R_RELOC:
                break
            ring, moves = _relocate(model, X, prob, ring)
            rl["relocations"] = moves
            if not moves:
                break
            start_moves.extend(moves)
            model.set_ring(ring)
        log["passed"] = bool(ver and ver.passed)
        outcome.starts.append(log)
        if ver is not None and ver.passed:
            outcome.success = True
            outcome.fracs = S.canonical_frac(X @ model.Minv)
            outcome.ring_ions = ring
            outcome.start_index = k
            outcome.relocations = start_moves
            outcome.verification = ver
            return outcome
        outcome.verification = ver
    return outcome
