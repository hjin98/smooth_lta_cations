"""v7.16 cage-salt occlusion: domain rules (D3 r3 "lta_cage_salt.py (domain)").

Scientific authority: doc/specs/v7_16_molten_salt_occlusion_extension.md (D1 r8).
Numerical authority: doc/specs/v7_16_cage_salt_numerical_method.md (D2 r3).
Architecture: doc/archs/lta_cation_generator_architecture.md, "v7.16 architecture delta".

This module owns request validation, the D1/D2 constants, the window table semantics,
exact charge accounting and rebalancing, the effective-volume estimator, coordination
floors, exact verification (V1-V6), stream derivation and the run record. It must not
import the population or template modules; geometry arrives as plain data.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import numbers
import platform
from dataclasses import dataclass, field, replace
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

# --------------------------------------------------------------------------------------
# D1-bound constants: name -> (value, authority locator, exact literal in that source).
# tests/test_v716_cage_salt.py asserts each literal occurs in the named source.
# --------------------------------------------------------------------------------------

D1_SPEC = "doc/specs/v7_16_molten_salt_occlusion_extension.md"
E_OCC = "doc/evidence/lta_salt_occlusion_evidence.md"

D1_BOUND: Dict[str, Tuple[Any, str, str]] = {
    "T_BULK_Li": (4.3, D1_SPEC + " §7", "4.3 ± 0.5 (§1a, inside r_min1)"),
    "T_BULK_Na": (5.2, D1_SPEC + " §7", "| Na | 5.2 ± 0.5 (§1a, inside r_min1)"),
    "T_BULK_K": (5.2, D1_SPEC + " §7", "| K | 5.2 ± 0.5 (§1a, inside r_min1)"),
    "T_BULK_Rb": (5.4, D1_SPEC + " §7", "5.4 ± 0.5 (§1a; r_min1 here is an E ratio estimate)"),
    "T_BULK_Cs": (5.5, D1_SPEC + " §7", "5.5 ± 0.5 (§1a; r_min1 here is an E ratio estimate)"),
    "T_BULK_Sr": (6.0, D1_SPEC + " §7", "6 ± 1 (E working value, dilute; not defined inside r_min1)"),
    "T_BULK_Ba": (7.0, D1_SPEC + " §7", "7 ± 1 (E working value; 6.4 vs 7.7 conflict open)"),
    "U_M": ({"Li": 6, "Na": 6, "K": 7, "Rb": 7, "Cs": 6, "Sr": 8, "Ba": 9}, D1_SPEC + " §7",
            "| U_M | 6 | 6 | 7 | 7 | 6 | 8 | 9 |"),
    "JANZ_LiCl": ((1.8842, 4.33e-4), D1_SPEC + " §7", "ρ_LiCl = 1.8842 − 4.33×10⁻⁴ T"),
    "JANZ_KCl": ((2.1359, 5.83e-4), D1_SPEC + " §7", "ρ_KCl = 2.1359 − 5.83×10⁻⁴ T"),
    "T_REF_K": (800.0, D1_SPEC + " §7", "At 800 K."),
    "X_LI_FALLBACK": (0.592, D1_SPEC + " §7", "eutectic x_Li = 0.592 is used"),
    "F_REF": (0.97, D1_SPEC + " §9", "f_ref = 0.97 ± 0.03"),
    "N_REF": (10.5, D1_SPEC + " §9", "n_ref = 10.5 Cl per pseudo-cell"),
    "A_RULE_B": (0.4, D1_SPEC + " §7", "**a = 0.4 (stakeholder decision S-11).**"),
    "A_RULE_A": (0.7, D1_SPEC + " §7", "go back to A\n  with a = 0.7."),
    "D_WALL_HALF_RANGE": (0.3, D1_SPEC + " §7", "**Uncertainty ±0.3 Å.**"),
    "DOPANT_FLAG_FRACTION": (0.2, D1_SPEC + " §7", "If dopants exceed 20 % of a cage's"),
    "ALPHA_CL_CAP": (12, D1_SPEC + " §4", "0 <= N_Cl(alpha_k) <= 12"),
    "T_MIN": (3.0, D1_SPEC + " §6", "**T_min = 3.0 Å**"),
    "M_Li": (6.941, E_OCC + " §5", "Li 6.941"),
    "M_K": (39.098, E_OCC + " §5", "K 39.098"),
    "M_Cl": (35.453, E_OCC + " §5", "Cl 35.453"),
}

T_BULK = {m: D1_BOUND[f"T_BULK_{m}"][0] for m in ("Li", "Na", "K", "Rb", "Cs", "Sr", "Ba")}
U_M: Dict[str, int] = dict(D1_BOUND["U_M"][0])
SALT_SPECIES = ("Li", "Na", "K", "Rb", "Cs", "Sr", "Ba")  # D1 §5 set S
F_REF = D1_BOUND["F_REF"][0]
N_REF = D1_BOUND["N_REF"][0]
A_RULE_B = D1_BOUND["A_RULE_B"][0]
A_RULE_A = D1_BOUND["A_RULE_A"][0]
X_LI_FALLBACK = D1_BOUND["X_LI_FALLBACK"][0]
T_REF_K = D1_BOUND["T_REF_K"][0]
ALPHA_CL_CAP = D1_BOUND["ALPHA_CL_CAP"][0]
D_WALL_HALF_RANGE = D1_BOUND["D_WALL_HALF_RANGE"][0]
DOPANT_FLAG_FRACTION = D1_BOUND["DOPANT_FLAG_FRACTION"][0]
T_MIN = D1_BOUND["T_MIN"][0]  # D1 r8 S-13: every cage_salt ion >= 3.0 A from framework Si/Al
T_ELEMENTS = ("Si", "Al")

# Release default coordination rule (D1 §7 S-11; S-12 fallback decided at qualification).
COORDINATION_RULE = "B"

# --------------------------------------------------------------------------------------
# D2-bound constants: the D2 r3 §10 register; each line names its D2 section.
# --------------------------------------------------------------------------------------

DELTA_B = 0.01  # D2 §6.1 build margin
TAU_W = 1.0e-5  # D2 §7 writer band
EPS_M = 0.05  # D2 §6.1 membership margin
DELTA_CN = 0.05  # D2 §6.1 count dead zone
W_S = 0.12  # D2 §6.1 smooth-count width
MU = 0.3  # D2 §6.1 count margin
W_STAGES = (10.0, 100.0, 1000.0)  # D2 §6.3 continuation
W_DZ_POLISH = 10.0  # D2 §6.3 polish dead-zone weight
LBFGS_OPTIONS = {"maxiter": 6000, "maxfun": 15000, "maxcor": 10, "ftol": 2.2e-9, "gtol": 1e-5, "maxls": 20}  # D2 §6.3
M_INS = 2000  # D2 §6.2 insertion candidates
K_MAX_DEFAULT = 20  # D2 §6.5 restarts
R_RELOC = 4  # D2 §6.4 relocation rounds
V_EFF_GRID = 192  # D2 §4.2 grid
TIE_TOL = 1.0e-9  # D2 §4.2 tie sharing
CEIL_GUARD = 1.0e-9  # D2 §2 ceiling guard
CANONICAL_DECIMALS = 10  # D2 §7 canonical coordinates
GUARD_DISTANCE = 1.5  # D1 §6 catastrophic guard (D2 §6.1)
SINGLE_IMAGE_MARGIN = 0.1  # D2 §2 single-image precondition
Z_MAX_NO_CONTRIBUTOR = 1  # D2 §4.3
EIGHTR_EQUIDISTANCE_TOL = 1.0e-6
VARIANT_KEY_DECIMALS = 6
VARIANT_MATCH_TOL = 1.0e-8  # fractional; canonical coordinates carry 10 decimals
REPORT_TOP_ALTERNATIVES = 5  # D2 §5 report

STREAM_TAG_REBALANCE = "v7.16-charge-rebalance:{seed}"  # D2 §3 streams
STREAM_TAG_START = "v7.16-cage-salt:{seed}:start{k}"  # D2 §3 per-start substreams

CAGE_IDS = ("alpha-0", "alpha-1", "beta-0", "beta-1")
ALPHA_IDS = ("alpha-0", "alpha-1")
BETA_IDS = ("beta-0", "beta-1")
RINGS = ("4", "6", "8")


class CageSaltInputError(ValueError):
    """An inadmissible cage-salt request (rejected up front; nothing is written)."""


class CageSaltConstructionError(RuntimeError):
    """Fail-closed construction failure; carries the run record."""

    def __init__(self, message: str, record: "CageSaltRecord"):
        super().__init__(message)
        self.record = record


def _data():
    import lta_cage_salt_data  # generated; lazy so regeneration can bootstrap

    return lta_cage_salt_data


# --------------------------------------------------------------------------------------
# Request model (D1 §4-§5)
# --------------------------------------------------------------------------------------

def parse_salt_spec(spec: Optional[str]) -> Dict[str, int]:
    """Parse ``"Cl:10,Li:6,K:4"`` into exact integer counts."""
    out: Dict[str, int] = {}
    if spec is None:
        return out
    text = str(spec).strip()
    if not text:
        return out
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise CageSaltInputError(f"Salt spec item {item!r} must look like Element:count")
        el, cnt = item.split(":", 1)
        el = el.strip()
        el = el[:1].upper() + el[1:].lower()
        cnt = cnt.strip()
        if not cnt.isdigit():
            raise CageSaltInputError(f"Salt count for {el} must be a nonnegative integer, got {cnt!r}")
        if el in out:
            raise CageSaltInputError(f"Element {el} repeated in salt spec {spec!r}")
        out[el] = int(cnt)
    return out


@dataclass(frozen=True)
class CageSaltRequest:
    """Per-cage salt inventory plus ring-site counterions (D1 §4)."""

    cages: Mapping[str, Mapping[str, int]]
    counterions: Mapping[str, Mapping[str, int]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Unknown keys would otherwise be ignored by cage()/is_active and a request
        # could silently resolve to a dry run.
        for c in self.cages:
            if c not in CAGE_IDS:
                raise CageSaltInputError(f"Unknown cage {c!r}; cages are {CAGE_IDS}")
        for r in self.counterions:
            if r not in RINGS:
                raise CageSaltInputError(f"Unknown ring family {r!r} for salt counterions; families are {RINGS}")
        for where, table in [*((c, self.cages[c]) for c in self.cages),
                             *((f"{r}R counterions", self.counterions[r]) for r in self.counterions)]:
            for el, n in table.items():
                if isinstance(n, bool) or not isinstance(n, numbers.Integral) or n < 0:
                    raise CageSaltInputError(f"{where}: count for {el} must be a nonnegative integer, got {n!r}")

    @classmethod
    def from_specs(
        cls,
        *,
        alpha: Optional[str] = None,
        beta: Optional[str] = None,
        alpha0: Optional[str] = None,
        alpha1: Optional[str] = None,
        beta0: Optional[str] = None,
        beta1: Optional[str] = None,
        counterion_n4: Optional[str] = None,
        counterion_n6: Optional[str] = None,
        counterion_n8: Optional[str] = None,
    ) -> "CageSaltRequest":
        a, b = parse_salt_spec(alpha), parse_salt_spec(beta)
        cages = {
            "alpha-0": parse_salt_spec(alpha0) if alpha0 is not None else dict(a),
            "alpha-1": parse_salt_spec(alpha1) if alpha1 is not None else dict(a),
            "beta-0": parse_salt_spec(beta0) if beta0 is not None else dict(b),
            "beta-1": parse_salt_spec(beta1) if beta1 is not None else dict(b),
        }
        counter = {"4": parse_salt_spec(counterion_n4), "6": parse_salt_spec(counterion_n6),
                   "8": parse_salt_spec(counterion_n8)}
        return cls(cages={k: {e: n for e, n in v.items()} for k, v in cages.items()},
                   counterions={k: {e: n for e, n in v.items()} for k, v in counter.items()})

    def cage(self, cage_id: str) -> Dict[str, int]:
        return {e: int(n) for e, n in self.cages.get(cage_id, {}).items() if int(n)}

    def n_cl(self, cage_id: Optional[str] = None) -> int:
        ids = (cage_id,) if cage_id else CAGE_IDS
        return sum(self.cage(c).get("Cl", 0) for c in ids)

    def salt_cations(self, cage_id: str) -> Dict[str, int]:
        return {e: n for e, n in self.cage(cage_id).items() if e != "Cl"}

    def counterion_counts(self, ring: Optional[str] = None) -> Dict[str, int]:
        rings = (ring,) if ring else RINGS
        out: Dict[str, int] = {}
        for r in rings:
            for e, n in self.counterions.get(r, {}).items():
                if int(n):
                    out[e] = out.get(e, 0) + int(n)
        return out

    @property
    def salt_total(self) -> int:
        return sum(sum(self.cage(c).values()) for c in CAGE_IDS)

    @property
    def counterion_total(self) -> int:
        return sum(self.counterion_counts().values())

    @property
    def is_active(self) -> bool:
        """D1 §3: active iff any cage-salt input is nonzero."""
        return self.salt_total > 0

    def to_dict(self) -> Dict[str, Any]:
        return {"cages": {c: self.cage(c) for c in CAGE_IDS},
                "counterions": {r: {e: int(n) for e, n in self.counterions.get(r, {}).items() if int(n)}
                                for r in RINGS}}


def validate_request(request: CageSaltRequest, charges: Mapping[str, int],
                     framework_species: Iterable[str]) -> None:
    """D1 §4-§5 domain checks that need no geometry (raise CageSaltInputError)."""
    for c in request.cages:
        if c not in CAGE_IDS:
            raise CageSaltInputError(f"Unknown cage {c!r}; cages are {CAGE_IDS}")
    for r in request.counterions:
        if r not in RINGS:
            raise CageSaltInputError(f"Unknown ring family {r!r} for salt counterions; families are {RINGS}")
    if not request.is_active:
        if request.counterion_total:
            raise CageSaltInputError("Salt counterions were given but every salt count is zero")
        return
    for c in CAGE_IDS:
        for el, n in request.cage(c).items():
            if el != "Cl" and el not in SALT_SPECIES:
                raise CageSaltInputError(
                    f"{el} in {c} is outside the cage-salt species set S={SALT_SPECIES} "
                    "(no bound pair data; D1 §5)")
        cl = request.cage(c).get("Cl", 0)
        cations = request.salt_cations(c)
        if c.startswith("alpha") and cl > ALPHA_CL_CAP:
            raise CageSaltInputError(f"{c}: {cl} Cl exceeds the D-2 cap of {ALPHA_CL_CAP} per alpha-cage")
        if c.startswith("beta"):
            if cl > 1:
                raise CageSaltInputError(f"{c}: beta-cages hold 0 or 1 Cl (D-2), got {cl}")
            if sum(cations.values()) > cl:
                raise CageSaltInputError(f"{c}: beta salt cations ({sum(cations.values())}) exceed N_Cl ({cl})")
        if cations and cl == 0:
            raise CageSaltInputError(f"{c}: salt cations without Cl can never be anchored (D1 §4)")
    for ring in RINGS:
        for el in request.counterions.get(ring, {}):
            if el not in SALT_SPECIES:
                raise CageSaltInputError(f"Counterion {el} on {ring}R is outside the species set S")
    for el in framework_species:
        if el not in SALT_SPECIES:
            raise CageSaltInputError(
                f"Framework cation {el} is outside the species set S; with the cage-salt model active "
                "every cation must have bound pair data (D1 §5)")
    q_exc = sum(charges[e] * n for c in CAGE_IDS for e, n in request.salt_cations(c).items())
    q_exc += sum(charges[e] * n for e, n in request.counterion_counts().items())
    if q_exc != request.n_cl():
        raise CageSaltInputError(
            f"Excess charge invariant violated: cations +{q_exc} != {request.n_cl()} Cl- (D1 §4)")


# --------------------------------------------------------------------------------------
# Pair tables (D1 §6, §7)
# --------------------------------------------------------------------------------------

def _pair(i: str, j: str, env: str) -> Tuple[float, float, float]:
    rows = _data().PAIR_ROWS
    if (i, j, env) in rows:
        return rows[(i, j, env)]
    if (j, i, env) in rows:
        return rows[(j, i, env)]
    raise CageSaltInputError(f"No bound pair row for ({i}, {j}, {env!r}); fail closed (D-3)")


def cl_weights(salt_cations: Sequence[str], all_cations: Sequence[str]) -> Dict[str, float]:
    """D1 §6 x_M: mole fractions over all cage-salt cations of the cell, else all cations."""
    pool = list(salt_cations) or list(all_cations)
    if not pool:
        return {"Li": X_LI_FALLBACK, "K": 1.0 - X_LI_FALLBACK}
    out: Dict[str, float] = {}
    for el in pool:
        out[el] = out.get(el, 0.0) + 1.0 / len(pool)
    return out


class PairTable:
    """r_cut / r_min1 lookups for one structure (fixed x_M weights)."""

    def __init__(self, weights: Mapping[str, float]):
        self.weights = dict(weights)
        self._cache: Dict[Tuple[str, str], float] = {}

    def r_cut(self, a: str, b: str) -> float:
        key = (a, b) if a <= b else (b, a)
        if key in self._cache:
            return self._cache[key]
        if a == "Cl" and b == "Cl":
            v = sum(x * _pair("Cl", "Cl", f"around {m}")[0] for m, x in self.weights.items())
        elif "Cl" in (a, b) and "O" in (a, b):
            v = sum(x * _pair("Cl", "O", f"around {m}")[0] for m, x in self.weights.items())
        elif "O" in (a, b):
            m = b if a == "O" else a
            v = _pair(m, "O", "")[0]
        elif "Cl" in (a, b):
            m = b if a == "Cl" else a
            v = _pair(m, "Cl", "")[0]
        else:
            v = _pair(a, b, "Cl-bridged")[0]
        self._cache[key] = float(v)
        return float(v)

    @staticmethod
    def r_min1(cation: str) -> float:
        return float(_pair(cation, "Cl", "")[2])


def d_wall_species() -> Dict[str, float]:
    """D1 §7 d_M = r_peak(Cl–O around M) − r_peak(Cl–Cl around M)/2."""
    return {m: _pair("Cl", "O", f"around {m}")[1] - _pair("Cl", "Cl", f"around {m}")[1] / 2.0
            for m in SALT_SPECIES}


def max_r_min1() -> float:
    return max(PairTable.r_min1(m) for m in SALT_SPECIES)


# --------------------------------------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------------------------------------

def min_image(dfrac: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    d = np.asarray(dfrac, float)
    d = d - np.round(d)
    return d @ matrix


def canonical_frac(frac: np.ndarray) -> np.ndarray:
    """D2 §7: round to 10 decimals, wrap into [0, 1), normalise -0.0."""
    f = np.round(np.asarray(frac, float), CANONICAL_DECIMALS)
    f = f - np.floor(f)
    f = np.round(f, CANONICAL_DECIMALS)
    f[f >= 1.0] = 0.0
    f = f + 0.0
    return f


def variant_set_key(variant_fracs: Iterable[Sequence[float]]) -> Tuple[Tuple[float, float, float], ...]:
    """D3 r3: wrap to [0,1), round to 6 decimals, map 1.0/-0.0 to 0.0, sort."""
    out = []
    for f in variant_fracs:
        g = np.asarray(f, float) % 1.0
        g = np.round(g, VARIANT_KEY_DECIMALS)
        g[g >= 1.0] = 0.0
        g = g + 0.0
        out.append(tuple(float(x) for x in g))
    return tuple(sorted(out))


def perpendicular_widths(matrix: np.ndarray) -> np.ndarray:
    a, b, c = np.asarray(matrix, float)
    vol = abs(float(np.dot(a, np.cross(b, c))))
    return np.array([vol / np.linalg.norm(np.cross(b, c)), vol / np.linalg.norm(np.cross(c, a)),
                     vol / np.linalg.norm(np.cross(a, b))])


def single_image_ok(matrix: np.ndarray) -> Tuple[bool, float, float]:
    half = float(perpendicular_widths(matrix).min() / 2.0)
    need = max_r_min1() + SINGLE_IMAGE_MARGIN
    return half > need, half, need


# --------------------------------------------------------------------------------------
# Window table (D2 §4.3, D3 r3)
# --------------------------------------------------------------------------------------

def window_table() -> Mapping[str, Tuple[str, Tuple[str, ...], Tuple]]:
    return _data().WINDOW_TABLE


def site_weights(site_id: str) -> Dict[str, Tuple[int, int]]:
    """Window weights w(j, c) as exact fractions (numerator, denominator)."""
    _, cages, _ = window_table()[site_id]
    if len(cages) == 1:
        return {cages[0]: (1, 1)}
    return {c: (1, 2) for c in cages}


def site_class(site_id: str) -> str:
    ring, cages, _ = window_table()[site_id]
    return f"{ring}R@" + "+".join(cages)


CLASS_ORDER = (
    "4R@alpha-0", "4R@alpha-1",
    "6R@beta-0+alpha-0", "6R@beta-0+alpha-1", "6R@beta-1+alpha-0", "6R@beta-1+alpha-1",
    "8R@alpha-0+alpha-1",
)


def class_weights_quarters(cls: str) -> Dict[str, int]:
    """Weights in units of 1/4 (w=1 -> 4, w=1/2 -> 2)."""
    cages = cls.split("@", 1)[1].split("+")
    return {c: (4 if len(cages) == 1 else 2) for c in cages}


def check_site_join(physical_sites_by_ring: Mapping[str, Sequence[Any]], matrix: np.ndarray,
                    alpha_centres: Mapping[str, Sequence[float]]) -> None:
    """Fail closed unless the loaded sites equal the window table (D3 r3 template gate)."""
    table = window_table()
    seen = set()
    for ring in RINGS:
        for site in physical_sites_by_ring.get(ring, []):
            entry = table.get(site.site_id)
            if entry is None:
                raise CageSaltInputError(
                    f"Site {site.site_id} is not in the canonical window table; the cage-salt model "
                    "requires the canonical template and default grouping options")
            key = variant_set_key([v.frac for v in site.variants])
            if entry[0] != ring or tuple(entry[2]) != key:
                raise CageSaltInputError(
                    f"Site {site.site_id} variants do not match the canonical window table")
            seen.add(site.site_id)
            if ring == "8":
                a0 = np.asarray(alpha_centres["alpha-0"], float)
                a1 = np.asarray(alpha_centres["alpha-1"], float)
                for v in site.variants:
                    d0 = np.linalg.norm(min_image(np.asarray(v.frac) - a0, matrix))
                    d1 = np.linalg.norm(min_image(np.asarray(v.frac) - a1, matrix))
                    if abs(d0 - d1) > EIGHTR_EQUIDISTANCE_TOL:
                        raise CageSaltInputError(
                            f"8R site {site.site_id} is not equidistant from alpha-0/alpha-1 "
                            f"({d0:.6f} vs {d1:.6f} A)")
    if seen != set(table):
        missing = sorted(set(table) - seen)[:5]
        raise CageSaltInputError(f"Loaded sites do not cover the window table (missing e.g. {missing})")


# --------------------------------------------------------------------------------------
# Charge (D1 §4, D2 §4.3) in exact units of 1/4 e
# --------------------------------------------------------------------------------------

def cage_charges_quarters(
    request: CageSaltRequest,
    ring_ions: Sequence[Any],
    charges: Mapping[str, int],
    *,
    include_framework: bool = True,
) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]:
    """Return (q_exc, q_fw, z_max) per cage; charges in units of e/4."""
    q_exc = {c: 4 * (sum(charges[e] * n for e, n in request.salt_cations(c).items())
                     - request.cage(c).get("Cl", 0)) for c in CAGE_IDS}
    fw_raw = {c: 0 for c in CAGE_IDS}
    zmax = {c: max([charges[e] for e in request.salt_cations(c)] or [0]) for c in CAGE_IDS}
    for ion in ring_ions:
        z = charges[ion.element]
        for c, (num, den) in site_weights(ion.site_id).items():
            contrib = z * num * (4 // den)
            if ion.role == "occlusion_counterion":
                q_exc[c] += contrib
                zmax[c] = max(zmax[c], z)
            elif include_framework:
                fw_raw[c] += contrib
                zmax[c] = max(zmax[c], z)
    q_fw = {c: 0 for c in CAGE_IDS}
    if include_framework:
        for group in (ALPHA_IDS, BETA_IDS):
            tot = sum(fw_raw[c] for c in group)
            for c in group:
                # q_fw = raw - mean; with two cages: (raw_c - raw_other)/2, exact in 1/8 -> keep 1/4 by
                # construction: raw values are multiples of 2 (w = 1/2 or 1 -> 2 or 4 quarters).
                q_fw[c] = fw_raw[c] - tot // 2 if tot % 2 == 0 else None  # type: ignore[assignment]
                if q_fw[c] is None:
                    raise AssertionError("framework charge share is not a multiple of e/4")
    zmax = {c: (z if z > 0 else Z_MAX_NO_CONTRIBUTOR) for c, z in zmax.items()}
    return q_exc, q_fw, zmax


def charge_violations(q_exc: Mapping[str, int], q_fw: Mapping[str, int], zmax: Mapping[str, int]) -> Dict[str, int]:
    """Return per-cage excess of |q| over q_max (quarters); 0 means satisfied."""
    return {c: max(0, abs(q_exc[c] + q_fw[c]) - 4 * zmax[c]) for c in CAGE_IDS}


# --------------------------------------------------------------------------------------
# Exact rebalancing (D2 §5)
# --------------------------------------------------------------------------------------

def _compositions(n: int, k: int) -> Iterable[Tuple[int, ...]]:
    if k == 1:
        yield (n,)
        return
    for first in range(n + 1):
        for rest in _compositions(n - first, k - 1):
            yield (first,) + rest


@dataclass
class RebalanceResult:
    feasible: bool
    vector: Dict[str, Dict[str, int]]
    objective: Tuple[int, int]
    alternatives: List[Dict[str, Any]]
    moves: List[Dict[str, Any]]
    best_violation_quarters: int
    blocking_sites: List[str]
    placed: List[Any]


def _classes_by_family() -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {r: [] for r in RINGS}
    for cls in CLASS_ORDER:
        out[cls[0]].append(cls)
    return out


def _search(request: CageSaltRequest, counter_groups: Mapping[Tuple[str, int], int],
            capacity: Mapping[str, int], base_q: Mapping[str, int], base_z: Mapping[str, int],
            n0: Mapping[Tuple[str, int], Mapping[str, int]], charges: Mapping[str, int]):
    fam_classes = _classes_by_family()
    groups = sorted(counter_groups)  # by (family, charge)
    per_family: Dict[str, List[Tuple[Dict[Tuple[str, int], Tuple[int, ...]], Dict[str, int]]]] = {}
    for fam in RINGS:
        classes = fam_classes[fam]
        fam_groups = [g for g in groups if g[0] == fam]
        options = []
        split_lists = [list(_compositions(counter_groups[g], len(classes))) for g in fam_groups]
        for combo in itertools.product(*split_lists) if fam_groups else [()]:
            used = [sum(c[i] for c in combo) for i in range(len(classes))]
            if any(used[i] > capacity.get(classes[i], 0) for i in range(len(classes))):
                continue
            q = {c: 0 for c in CAGE_IDS}
            zc = {c: 0 for c in CAGE_IDS}
            for g, split in zip(fam_groups, combo):
                for cls, cnt in zip(classes, split):
                    if cnt:
                        for cage, wq in class_weights_quarters(cls).items():
                            q[cage] += g[1] * wq * cnt
                            zc[cage] = max(zc[cage], g[1])
            options.append(({g: s for g, s in zip(fam_groups, combo)}, q, zc))
        per_family[fam] = options
    best = None
    ranked: List[Tuple] = []
    best_violation = None
    worst_cages: List[str] = []
    for choice in itertools.product(*(per_family[f] for f in RINGS)):
        q = dict(base_q)
        z = dict(base_z)
        splits: Dict[Tuple[str, int], Tuple[int, ...]] = {}
        for sp, qf, zf in choice:
            splits.update(sp)
            for c in CAGE_IDS:
                q[c] += qf[c]
                z[c] = max(z[c], zf[c])
        z = {c: (v if v > 0 else Z_MAX_NO_CONTRIBUTOR) for c, v in z.items()}
        # D2 §5 reports min over vectors of max_c(|q_c| - q_max(c)); feasible iff it is <= 0
        viol = max(abs(q[c]) - 4 * z[c] for c in CAGE_IDS)
        if best_violation is None or viol < best_violation:
            best_violation = viol
            worst_cages = [c for c in CAGE_IDS if abs(q[c]) - 4 * z[c] > 0]
        if viol > 0:
            continue
        moved = 0
        for g, split in splits.items():
            classes = fam_classes[g[0]]
            for cls, cnt in zip(classes, split):
                moved += max(0, n0.get(g, {}).get(cls, 0) - cnt)
        sumsq = sum(v * v for v in q.values())
        vec_key = tuple(cnt for g in groups for cnt in splits[g])
        key = (moved, sumsq, vec_key)
        ranked.append((key, splits, q))
        if best is None or key < best[0]:
            best = (key, splits, q)
    ranked.sort(key=lambda t: t[0])
    return best, ranked, max(0, best_violation if best_violation is not None else 0), worst_cages


def upfront_charge_check(request: CageSaltRequest, charges: Mapping[str, int]) -> None:
    """D1 §4 deliberately stricter input rule: balanceable on excess charge alone (D2 §5)."""
    if not request.is_active:
        return
    table = window_table()
    capacity: Dict[str, int] = {}
    for sid in table:
        cls = site_class(sid)
        capacity[cls] = capacity.get(cls, 0) + 1
    groups: Dict[Tuple[str, int], int] = {}
    for ring in RINGS:
        for el, n in request.counterions.get(ring, {}).items():
            if n:
                g = (ring, charges[el])
                groups[g] = groups.get(g, 0) + int(n)
    base_q = {c: 4 * (sum(charges[e] * n for e, n in request.salt_cations(c).items())
                      - request.cage(c).get("Cl", 0)) for c in CAGE_IDS}
    base_z = {c: max([charges[e] for e in request.salt_cations(c)] or [0]) for c in CAGE_IDS}
    best, _, viol, cages = _search(request, groups, capacity, base_q, base_z, {}, charges)
    if best is None:
        if cages:
            detail = f"best excess {viol / 4:.2f} e in {', '.join(cages)}"
        else:
            detail = "no class vector fits the ring-site capacity"
        raise CageSaltInputError(
            "No placement of the requested salt counterions can keep every cage charge within "
            f"q_max ({detail}); rejected up front (D1 §4, D2 §5)")


def rebalance(
    request: CageSaltRequest,
    placed: Sequence[Any],
    physical_sites_by_ring: Mapping[str, Sequence[Any]],
    charges: Mapping[str, int],
    rng: np.random.Generator,
) -> RebalanceResult:
    """Step 4: exhaustive class-count search, then deterministic-up-to-R_q realization."""
    sites = {s.site_id: s for ring in RINGS for s in physical_sites_by_ring.get(ring, [])}
    occupied_fw = {ion.site_id for ion in placed if ion.role == "framework"}
    counter = [ion for ion in placed if ion.role == "occlusion_counterion"]
    capacity: Dict[str, int] = {}
    for sid in sites:
        if sid not in occupied_fw:
            cls = site_class(sid)
            capacity[cls] = capacity.get(cls, 0) + 1
    groups: Dict[Tuple[str, int], int] = {}
    n0: Dict[Tuple[str, int], Dict[str, int]] = {}
    for ion in counter:
        g = (ion.ring_class, charges[ion.element])
        groups[g] = groups.get(g, 0) + 1
        cls = site_class(ion.site_id)
        n0.setdefault(g, {})
        n0[g][cls] = n0[g].get(cls, 0) + 1
    q_exc_salt = {c: 4 * (sum(charges[e] * n for e, n in request.salt_cations(c).items())
                          - request.cage(c).get("Cl", 0)) for c in CAGE_IDS}
    fw_only = [ion for ion in placed if ion.role == "framework"]
    _, q_fw, z_fw = cage_charges_quarters(CageSaltRequest(cages={}), fw_only, charges)
    base_q = {c: q_exc_salt[c] + q_fw[c] for c in CAGE_IDS}
    base_z = {c: max([charges[e] for e in request.salt_cations(c)] + [0]) for c in CAGE_IDS}
    for ion in fw_only:
        for c in site_weights(ion.site_id):
            base_z[c] = max(base_z[c], charges[ion.element])
    best, ranked, viol, _ = _search(request, groups, capacity, base_q, base_z, n0, charges)
    if best is None:
        # framework-held sites of the counterion families: freeing any of them widens the search
        families = {g[0] for g in groups}
        blocking = sorted(ion.site_id for ion in fw_only if ion.ring_class in families)
        return RebalanceResult(False, {}, (0, 0), [], [], viol, blocking, list(placed))
    key, splits, q = best
    fam_classes = _classes_by_family()
    target: Dict[Tuple[str, int], Dict[str, int]] = {
        g: {cls: cnt for cls, cnt in zip(fam_classes[g[0]], split)} for g, split in splits.items()}
    new_placed = list(placed)
    occupied = {ion.site_id for ion in new_placed}
    moves: List[Dict[str, Any]] = []
    for fam in RINGS:
        fam_groups = sorted(g for g in target if g[0] == fam)
        leaving: Dict[Tuple[str, int], List[int]] = {}
        before: Dict[Tuple[str, int], Dict[str, int]] = {}
        for g in fam_groups:
            members = [i for i, ion in enumerate(new_placed)
                       if ion.role == "occlusion_counterion" and ion.ring_class == fam and charges[ion.element] == g[1]]
            by_cls: Dict[str, List[int]] = {}
            for i in members:
                by_cls.setdefault(site_class(new_placed[i].site_id), []).append(i)
            before[g] = {cls: len(v) for cls, v in by_cls.items()}
            leaving[g] = []
            for cls in fam_classes[fam]:
                cur = sorted(by_cls.get(cls, []), key=lambda i: (new_placed[i].site_id, new_placed[i].element))
                surplus = len(cur) - target[g].get(cls, 0)
                if surplus > 0:
                    leaving[g].extend(cur[:surplus])
        # Phase 1: every surplus counterion of the family vacates its site, so a deficit
        # class is never blocked by another charge group that leaves later (D4 realization
        # of D2 §5; groups of one family share the same physical sites).
        for g in fam_groups:
            for i in leaving[g]:
                occupied.discard(new_placed[i].site_id)
        # Phase 2: fill deficits in canonical (group, class) order with R_q draws.
        for g in fam_groups:
            for cls in fam_classes[fam]:
                deficit = target[g].get(cls, 0) - before[g].get(cls, 0)
                for _ in range(max(0, deficit)):
                    i = leaving[g].pop(0)
                    empty = sorted(sid for sid in sites if site_class(sid) == cls and sid not in occupied)
                    if not empty:
                        raise RuntimeError(f"rebalancing realization found no empty {cls} site (D2 §5)")
                    sid = empty[int(rng.integers(len(empty)))]
                    site = sites[sid]
                    vidx = int(rng.integers(len(site.variants)))
                    var = site.variants[vidx]
                    old = new_placed[i]
                    occupied.add(sid)
                    new_placed[i] = replace(old, frac=np.asarray(var.frac, float), site_id=sid,
                                            site_label=getattr(var, "site_label", old.site_label),
                                            variant_id=int(getattr(var, "variant_id", vidx)))
                    moves.append({"ion": f"R{old.site_id}", "to_ion_id": f"R{sid}", "element": old.element,
                                  "from_site": old.site_id, "to_site": sid,
                                  "to_variant": int(getattr(var, "variant_id", vidx)),
                                  "from_class": site_class(old.site_id), "to_class": cls})
    alts = [{"objective": [k[0], k[1]], "q_quarters": dict(qq),
             "vector": {f"{g[0]}R/z{g[1]}": dict(zip(fam_classes[g[0]], sp)) for g, sp in sorted(spl.items())}}
            for k, spl, qq in ranked[1:1 + REPORT_TOP_ALTERNATIVES]]
    vector = {f"{g[0]}R/z{g[1]}": dict(target[g]) for g in sorted(target)}
    return RebalanceResult(True, vector, (key[0], key[1]), alts, moves, 0, [], new_placed)


# --------------------------------------------------------------------------------------
# Effective volume V_eff (D2 §4.2) with a d-independent cached field
# --------------------------------------------------------------------------------------

_FIELD_CACHE: Dict[str, Dict[str, Tuple[np.ndarray, np.ndarray]]] = {}


def _field(matrix: np.ndarray, o_fracs: np.ndarray, centres: Mapping[str, Sequence[float]]):
    from scipy.spatial import cKDTree

    key = hashlib.sha256(np.ascontiguousarray(np.round(o_fracs, 12)).tobytes()
                         + np.ascontiguousarray(matrix).tobytes()).hexdigest()
    if key in _FIELD_CACHE:
        return _FIELD_CACHE[key], key
    G = V_EFF_GRID
    g = (np.arange(G) + 0.5) / G
    F = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    X = F @ matrix
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)
    o_img = ((np.asarray(o_fracs, float)[None, :, :] + shifts[:, None, :]).reshape(-1, 3)) @ matrix
    d_o, _ = cKDTree(o_img).query(X)
    names = list(CAGE_IDS)
    c_fr = np.array([centres[c] for c in names], float)
    c_img = ((c_fr[None, :, :] + shifts[:, None, :]).reshape(-1, 3)) @ matrix
    labels = np.tile(np.arange(len(names)), len(shifts))
    dc, ic = cKDTree(c_img).query(X, k=4)
    lab = labels[ic]
    near = np.abs(dc - dc[:, :1]) <= TIE_TOL
    mult = near.sum(axis=1)
    out: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for idx, name in enumerate(names):
        w = ((lab == idx) & near).sum(axis=1) / mult
        sel = w > 0
        order = np.argsort(d_o[sel], kind="stable")
        dvals = d_o[sel][order]
        cum = np.cumsum(w[sel][order][::-1])[::-1]  # cum[i] = sum of weights with index >= i
        out[name] = (dvals, cum)
    _FIELD_CACHE[key] = out
    return out, key


def v_eff(ref_matrix: np.ndarray, o_fracs: np.ndarray, centres: Mapping[str, Sequence[float]],
          cage: str, d: float, scale: float = 1.0) -> float:
    """V_eff(scaled cell; d) = s^3 * V_eff(ref cell; d/s) (D2 §4.2 scale identity)."""
    fld, _ = _field(np.asarray(ref_matrix, float), np.asarray(o_fracs, float), centres)
    dvals, cum = fld[cage]
    idx = int(np.searchsorted(dvals, d / scale, side="right"))
    total = float(cum[idx]) if idx < len(cum) else 0.0
    vol = abs(float(np.linalg.det(ref_matrix)))
    return scale ** 3 * vol * total / V_EFF_GRID ** 3


def v_pair_800k(x_li: float) -> float:
    """Volume per ion pair (Å^3) of a Li/K chloride melt at 800 K, ideal mixing (D1 §7)."""
    (a_li, b_li), (a_k, b_k) = D1_BOUND["JANZ_LiCl"][0], D1_BOUND["JANZ_KCl"][0]
    m_li = D1_BOUND["M_Li"][0] + D1_BOUND["M_Cl"][0]
    m_k = D1_BOUND["M_K"][0] + D1_BOUND["M_Cl"][0]
    rho_li = a_li - b_li * T_REF_K
    rho_k = a_k - b_k * T_REF_K
    vm = x_li * m_li / rho_li + (1.0 - x_li) * m_k / rho_k  # cm^3 / mol
    return vm / 6.02214076e23 * 1e24


@dataclass
class CageFloors:
    cage: str
    n_cl: int
    x_li: float
    x_li_fallback: bool
    dopant_fraction: float
    dopant_flag: bool
    d_wall: float
    v_eff: float
    v_eff_range: Tuple[float, float]
    n_melt: float
    n_melt_range: Tuple[float, float]
    phi: float
    c_by_species: Dict[str, float]
    floor_by_species: Dict[str, int]
    a: float
    rule: str

    def to_dict(self) -> Dict[str, Any]:
        return {k: (list(v) if isinstance(v, tuple) else v) for k, v in self.__dict__.items()}


def cage_floors(request: CageSaltRequest, ref_matrix: np.ndarray, o_fracs: np.ndarray,
                centres: Mapping[str, Sequence[float]], scale: float,
                rule: str = COORDINATION_RULE) -> Dict[str, CageFloors]:
    a = A_RULE_B if rule == "B" else A_RULE_A
    dm = d_wall_species()
    out: Dict[str, CageFloors] = {}
    for c in ALPHA_IDS:
        cat = request.salt_cations(c)
        total = sum(cat.values())
        if total:
            d_wall = sum(dm[m] * n for m, n in cat.items()) / total
        else:
            d_wall = X_LI_FALLBACK * dm["Li"] + (1 - X_LI_FALLBACK) * dm["K"]
        lik = cat.get("Li", 0) + cat.get("K", 0)
        if lik:
            x_li, fb = cat.get("Li", 0) / lik, False
        else:
            x_li, fb = X_LI_FALLBACK, True
        dop = (total - lik) / total if total else 0.0
        vp = v_pair_800k(x_li)
        v0 = v_eff(ref_matrix, o_fracs, centres, c, d_wall, scale)
        vlo = v_eff(ref_matrix, o_fracs, centres, c, d_wall + D_WALL_HALF_RANGE, scale)
        vhi = v_eff(ref_matrix, o_fracs, centres, c, d_wall - D_WALL_HALF_RANGE, scale)
        n_melt = v0 / vp
        n_cl = request.cage(c).get("Cl", 0)
        phi = min(1.0, n_cl / n_melt) if n_melt > 0 else 1.0
        cvals = {m: T_BULK[m] * phi for m in cat}
        floors = {m: max(1, int(math.ceil(a * cvals[m] - CEIL_GUARD))) for m in cat}
        out[c] = CageFloors(c, n_cl, x_li, fb, dop, dop > DOPANT_FLAG_FRACTION, d_wall, v0, (vlo, vhi),
                            n_melt, (vlo / vp, vhi / vp), phi, cvals, floors, a, rule)
    return out


# --------------------------------------------------------------------------------------
# Streams (D2 §3)
# --------------------------------------------------------------------------------------

def _stream(tag: str) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(int.from_bytes(hashlib.sha256(tag.encode()).digest()[:8], "big")))


def rebalance_stream(seed: Optional[int]) -> Tuple[np.random.Generator, Optional[str]]:
    if seed is None:
        return np.random.Generator(np.random.PCG64()), None
    tag = STREAM_TAG_REBALANCE.format(seed=seed)
    return _stream(tag), tag


def start_stream(seed: Optional[int], k: int) -> Tuple[np.random.Generator, Optional[str]]:
    if seed is None:
        return np.random.Generator(np.random.PCG64()), None
    tag = STREAM_TAG_START.format(seed=seed, k=k)
    return _stream(tag), tag


# --------------------------------------------------------------------------------------
# Verification V1-V6 (D2 §7)
# --------------------------------------------------------------------------------------

@dataclass
class Atom:
    """Plain per-atom data handed to verification."""

    element: str
    role: str  # "framework_atom" | "framework" | "occlusion_counterion" | "cage_salt"
    ion_id: str
    cage_id: Optional[str] = None
    site_id: Optional[str] = None


@dataclass
class Verification:
    passed: bool
    checks: Dict[str, bool]
    detail: Dict[str, Any]
    failing: List[str]


def _floor_for(atom: Atom, floors: Mapping[str, CageFloors], rule: str) -> Tuple[int, Optional[int]]:
    if atom.element == "Cl":
        return 1, None
    if atom.cage_id in floors and rule == "B":
        return floors[atom.cage_id].floor_by_species[atom.element], U_M[atom.element]
    if atom.cage_id in floors:
        return 1, U_M[atom.element]
    return 1, None


def verify(
    matrix: np.ndarray,
    fracs: np.ndarray,
    atoms: Sequence[Atom],
    request: CageSaltRequest,
    *,
    pairs: PairTable,
    floors: Mapping[str, CageFloors],
    centres: Mapping[str, Sequence[float]],
    charges: Mapping[str, int],
    expected_ring_counts: Mapping[Tuple[str, str, str], int],
    site_variants: Mapping[str, Sequence[Sequence[float]]],
    rule: str = COORDINATION_RULE,
    tau: float = TAU_W,
    ring_ions: Optional[Sequence[Any]] = None,
) -> Verification:
    """Exact verification on the given coordinates (D2 §7; brute force, float64)."""
    M = np.asarray(matrix, float)
    F = np.asarray(fracs, float)
    checks: Dict[str, bool] = {}
    detail: Dict[str, Any] = {}
    salt = [i for i, a in enumerate(atoms) if a.role == "cage_salt"]
    ring = [i for i, a in enumerate(atoms) if a.role in ("framework", "occlusion_counterion")]
    fwat = [i for i, a in enumerate(atoms) if a.role == "framework_atom"]
    names = list(CAGE_IDS)
    cen = np.array([centres[c] for c in names], float)

    # V6 single image
    ok6, half, need = single_image_ok(M)
    checks["V6_single_image"] = ok6
    detail["V6"] = {"half_width": half, "required": need}

    # V1 inventory and membership
    counts: Dict[str, Dict[str, int]] = {c: {} for c in names}
    margins = []
    v1 = True
    for i in salt:
        d = np.linalg.norm(min_image(F[i][None, :] - cen, M), axis=1)
        order = np.argsort(d)
        home = names.index(atoms[i].cage_id)
        margin = float(min(d[j] for j in range(len(names)) if j != home) - d[home])
        margins.append(margin)
        if order[0] != home or margin < tau:
            v1 = False
        counts[atoms[i].cage_id][atoms[i].element] = counts[atoms[i].cage_id].get(atoms[i].element, 0) + 1
    for c in names:
        if counts[c] != request.cage(c):
            v1 = False
    ring_counts: Dict[Tuple[str, str, str], int] = {}
    sites_used: Dict[str, int] = {}
    for i in ring:
        a = atoms[i]
        fam = window_table()[a.site_id][0]
        key = (a.role, fam, a.element)
        ring_counts[key] = ring_counts.get(key, 0) + 1
        sites_used[a.site_id] = sites_used.get(a.site_id, 0) + 1
        # periodic, absolute-only comparison (no relative tolerance)
        if not any(float(np.abs((np.asarray(v, float) - F[i] + 0.5) % 1.0 - 0.5).max()) <= VARIANT_MATCH_TOL
                   for v in site_variants[a.site_id]):
            v1 = False
    if any(n > 1 for n in sites_used.values()):
        v1 = False
    if {k: v for k, v in ring_counts.items() if v} != {k: v for k, v in expected_ring_counts.items() if v}:
        v1 = False
    checks["V1_inventory"] = v1
    detail["V1"] = {"min_membership_margin": min(margins) if margins else None,
                    "salt_counts": counts, "ring_counts": {"/".join(k): v for k, v in ring_counts.items()}}

    # V2 hard minima
    v2 = True
    slack: Dict[str, float] = {}

    def upd(cls: str, s: float):
        slack[cls] = min(slack.get(cls, np.inf), s)

    for n_i, i in enumerate(salt):
        ei = atoms[i].element
        for j in salt[n_i + 1:]:
            d = float(np.linalg.norm(min_image(F[j] - F[i], M)))
            upd("salt-salt", d - pairs.r_cut(ei, atoms[j].element))
        if ring:
            dr = np.linalg.norm(min_image(F[ring] - F[i], M), axis=1)
            for d, j in zip(dr, ring):
                upd("salt-ring", float(d) - pairs.r_cut(ei, atoms[j].element))
        if fwat:
            df = np.linalg.norm(min_image(F[fwat] - F[i], M), axis=1)
            for d, j in zip(df, fwat):
                upd("guard", float(d) - GUARD_DISTANCE)
                if atoms[j].element == "O":
                    upd("salt-O", float(d) - pairs.r_cut(ei, "O"))
                elif atoms[j].element in T_ELEMENTS:
                    upd("salt-T", float(d) - T_MIN)
    for cls, s in slack.items():
        if s < tau:
            v2 = False
    checks["V2_hard_minima"] = v2
    detail["V2"] = {"min_slack": slack}

    # V3 anchoring, coordination, counting dead band
    v3 = True
    cn: Dict[str, Dict[str, Any]] = {}
    band_min = np.inf
    cl_idx = [i for i in salt if atoms[i].element == "Cl"]
    cat_idx = [i for i in range(len(atoms)) if atoms[i].role != "framework_atom" and atoms[i].element != "Cl"]
    per_group: Dict[Tuple[str, str], List[int]] = {}
    for i in salt:
        a = atoms[i]
        others = cat_idx if a.element == "Cl" else cl_idx
        n_nb = 0
        for j in others:
            rmin = pairs.r_min1(atoms[j].element if a.element == "Cl" else a.element)
            d = float(np.linalg.norm(min_image(F[j] - F[i], M)))
            band_min = min(band_min, abs(d - rmin))
            if d < rmin:
                n_nb += 1
        lo, hi = _floor_for(a, floors, rule)
        ok = n_nb >= lo and (hi is None or n_nb <= hi)
        if not ok:
            v3 = False
        cn[a.ion_id] = {"cn": n_nb, "floor": lo, "cap": hi, "ok": ok}
        if a.element != "Cl" and a.cage_id in floors:
            per_group.setdefault((a.cage_id, a.element), []).append(n_nb)
    if rule == "A":
        groups_ok, detail_a = rule_a_check(per_group, floors)
        v3 = v3 and groups_ok
        detail["V3_rule_A"] = detail_a
    if band_min < tau:
        v3 = False
    checks["V3_coordination"] = v3
    detail["V3"] = {"per_ion": cn, "min_distance_to_r_min1": float(band_min) if np.isfinite(band_min) else None,
                    "rule": rule}

    # V4 cage charge, V5 invariants
    if ring_ions is not None:
        q_exc, q_fw, zmax = cage_charges_quarters(request, ring_ions, charges)
        viol = charge_violations(q_exc, q_fw, zmax)
        checks["V4_cage_charge"] = all(v == 0 for v in viol.values())
        detail["V4"] = {"q_exc": {c: q_exc[c] / 4 for c in names}, "q_fw": {c: q_fw[c] / 4 for c in names},
                        "q_max": {c: zmax[c] for c in names}}
        q_fw_total = sum(charges[r.element] for r in ring_ions if r.role == "framework")
        q_exc_total = sum(charges[r.element] for r in ring_ions if r.role == "occlusion_counterion")
        q_exc_total += sum(charges[atoms[i].element] for i in salt if atoms[i].element != "Cl")
        n_cl = len(cl_idx)
        checks["V5_invariants"] = (q_fw_total == 24) and (q_exc_total == n_cl)
        detail["V5"] = {"Q_fw": q_fw_total, "Q_exc": q_exc_total, "N_Cl": n_cl}
    failing = [k for k, v in checks.items() if not v]
    return Verification(not failing, checks, detail, failing)


def rule_a_check(per_group: Mapping[Tuple[str, str], List[int]], floors: Mapping[str, CageFloors]):
    """Rule A (S-12 fallback): per-cage mean with pooling of small populations (D2 §6.6)."""
    detail: Dict[str, Any] = {}
    ok = True
    species = sorted({m for (_, m) in per_group})
    for m in species:
        alone, pooled = [], []
        for c in ALPHA_IDS:
            vals = per_group.get((c, m), [])
            if not vals:
                continue
            (alone if len(vals) >= 3 else pooled).append(c)
        for c in alone:
            vals = per_group[(c, m)]
            need = A_RULE_A * floors[c].c_by_species[m] * len(vals)
            good = sum(vals) >= need - CEIL_GUARD
            ok = ok and good
            detail[f"{c}/{m}"] = {"sum_cn": sum(vals), "required": need, "ok": good}
        if pooled:
            vals = [v for c in pooled for v in per_group[(c, m)]]
            if len(vals) >= 3:
                need = sum(A_RULE_A * floors[c].c_by_species[m] * len(per_group[(c, m)]) for c in pooled)
                good = sum(vals) >= need - CEIL_GUARD
                ok = ok and good
                detail[f"pooled/{m}"] = {"cages": pooled, "sum_cn": sum(vals), "required": need, "ok": good}
            else:
                detail[f"pooled/{m}"] = {"cages": pooled, "note": "fewer than 3 centres: anchoring and U_M only"}
    return ok, detail


# --------------------------------------------------------------------------------------
# Run record (D1 §14, D2 §8.1)
# --------------------------------------------------------------------------------------

def stack_versions() -> Dict[str, Any]:
    out: Dict[str, Any] = {"python": platform.python_version(), "platform": platform.platform(),
                           "numpy": np.__version__}
    try:
        import scipy

        out["scipy"] = scipy.__version__
    except ImportError:  # pragma: no cover
        out["scipy"] = None
    try:
        import pymatgen.core

        out["pymatgen"] = pymatgen.core.__version__
    except Exception:  # pragma: no cover
        out["pymatgen"] = None
    try:
        cfg = np.show_config(mode="dicts")  # numpy >= 1.25
        out["blas"] = cfg.get("Build Dependencies", {}).get("blas", {})
    except Exception:  # pragma: no cover
        out["blas"] = None
    import os

    out["threads"] = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
    return out


@dataclass
class CageSaltRecord:
    """Canonical realized record for one cage-salt run (schema_version 1)."""

    status: str
    request: Dict[str, Any]
    seed: Optional[int]
    rule: str
    volume: Dict[str, Any] = field(default_factory=dict)
    rebalance: Dict[str, Any] = field(default_factory=dict)
    floors: Dict[str, Any] = field(default_factory=dict)
    starts: List[Dict[str, Any]] = field(default_factory=list)
    relocations: List[Dict[str, Any]] = field(default_factory=list)
    verification: Dict[str, Any] = field(default_factory=dict)
    histograms: Dict[str, Any] = field(default_factory=dict)
    density_after: Dict[str, Any] = field(default_factory=dict)
    streams: Dict[str, Any] = field(default_factory=dict)
    parameters: Dict[str, Any] = field(default_factory=dict)
    stack: Dict[str, Any] = field(default_factory=dict)
    evidence_sha256: Dict[str, str] = field(default_factory=dict)
    message: str = ""
    traceback: Optional[str] = None
    schema_version: int = 1

    def to_json(self, atom_rows: Optional[List[Dict[str, Any]]] = None,
                structure_sha256: Optional[str] = None, structure_written: bool = False) -> str:
        payload = dict(self.__dict__)
        payload["structure_written"] = structure_written
        payload["structure_sha256"] = structure_sha256
        if not structure_written:
            payload["warning"] = ("No structure was written by this run; any file at the output path "
                                  "was not produced by it.")
        payload["atoms"] = atom_rows or []
        return json.dumps(payload, indent=1, sort_keys=True, default=_json_default) + "\n"


def _json_default(obj):
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, tuple):
        return list(obj)
    if obj == np.inf:
        return None
    raise TypeError(f"not JSON serializable: {type(obj)}")


def parameters_snapshot(k_max: int) -> Dict[str, Any]:
    return {"delta_b": DELTA_B, "tau_w": TAU_W, "eps_m": EPS_M, "delta_cn": DELTA_CN, "w_s": W_S, "mu": MU,
            "w_stages": list(W_STAGES), "w_dz_polish": W_DZ_POLISH, "lbfgs": dict(LBFGS_OPTIONS),
            "m_ins": M_INS, "k_max": k_max, "r_reloc": R_RELOC, "v_eff_grid": V_EFF_GRID,
            "a_rule_b": A_RULE_B, "a_rule_a": A_RULE_A, "f_ref": F_REF, "n_ref": N_REF, "t_min": T_MIN}


def histogram(values: Sequence[float], lo: float, hi: float, width: float) -> Dict[str, Any]:
    edges = np.arange(lo, hi + width / 2, width)
    counts, _ = np.histogram(np.asarray(values, float), bins=edges)
    return {"edges": [round(float(e), 6) for e in edges], "counts": [int(c) for c in counts]}


def f_v(n_cl_total: int) -> float:
    """D1 §9 default volume factor; n_bar = N_Cl / 2 Cl per pseudo-cell."""
    n_bar = n_cl_total / 2.0
    return 1.0 - (1.0 - F_REF) * min(n_bar / N_REF, 1.0)
