"""v7.16 cage-salt model: D4 stages S1-S4 and D2 §9 oracles 1-7.

Authorities: D1 doc/specs/v7_16_molten_salt_occlusion_extension.md (r8), D2
doc/specs/v7_16_cage_salt_numerical_method.md (r3), D3 "v7.16 architecture delta"
(r3). The independent verifier at the end of this file is written from D1/D2 and the
evidence CSV directly; it does not call ``lta_cage_salt``.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import re
import subprocess
import sys
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import lta_cage_salt as S
import lta_cage_salt_build as B
import lta_cation_population_grouped as P

ROOT = Path(__file__).resolve().parents[1]
REAL = ROOT / "tests" / "real_lta_primitive.cif"
FW = dict(n4="Na:2", n6="Na:16", n8="Na:6")
ALPHA = "Cl:12,Li:7,K:5"
BETA = "Cl:1,Li:1"
CENTRES = {c.cage_id: c.frac for c in P.CANONICAL_CAGE_SITES}


def request(alpha=ALPHA, beta=BETA, **kw):
    return S.CageSaltRequest.from_specs(alpha=alpha, beta=beta, **kw)


def generate(req=None, seed=3, **kw):
    params = dict(FW)
    params.update(kw)
    return P.generate_lta_cation_structure_from_specs(REAL, cage_salt=req or request(), seed=seed, **params)


@pytest.fixture(scope="module")
def salt_result():
    return generate()


@pytest.fixture(scope="module")
def counterion_result():
    # 4R counterions must be rebalanced: each alpha-cage salt carries -2 e.
    req = request(alpha="Cl:10,Li:4,K:4", counterion_n4="Li:4")
    return generate(req, seed=5)


# ----------------------------------------------------------------------------- S1


def test_generated_data_module_is_byte_identical_to_regeneration():
    sys.path.insert(0, str(ROOT / "validation"))
    try:
        import regenerate_cage_salt_data as R
    finally:
        sys.path.pop(0)
    assert R.render() == (ROOT / "lta_cage_salt_data.py").read_text()


def test_evidence_sha256_match_the_bound_csvs():
    for name, digest in S._data().EVIDENCE_CSV_SHA256.items():
        assert hashlib.sha256((ROOT / "doc" / "evidence" / name).read_bytes()).hexdigest() == digest


def _section(path: Path, locator: str) -> str:
    text = path.read_text()
    num = locator.split("§", 1)[1].strip()
    m = re.search(rf"^## {re.escape(num)}\. .*?(?=^## |\Z)", text, re.M | re.S)
    assert m, f"section {locator} not found in {path}"
    return m.group(0)


def _numbers(value):
    if isinstance(value, dict):
        return [x for v in value.values() for x in _numbers(v)]
    if isinstance(value, (tuple, list)):
        return [x for v in value for x in _numbers(v)]
    return [value]


def test_d1_bound_literals_occur_in_their_sections_and_parse_to_the_value():
    for name, (value, locator, literal) in S.D1_BOUND.items():
        source, sec = locator.split(" §")
        section = _section(ROOT / source, "§" + sec)
        assert literal in section, f"{name}: literal not found in {locator}"
        flat = literal.replace("×10⁻⁴", "e-4").replace("−", "-")
        found = [float(x) for x in re.findall(r"\d+(?:\.\d+)?(?:e-?\d+)?", flat)]
        found += [float(x) / 100.0 for x in re.findall(r"(\d+(?:\.\d+)?) %", flat)]
        for v in _numbers(value):
            assert any(abs(f - float(v)) <= 1e-12 * max(1.0, abs(float(v))) for f in found), (name, v, literal)


# ----------------------------------------------------------------------------- S2 request


@pytest.mark.parametrize(
    "kw,needle",
    [
        (dict(alpha="Cl:13,Li:13"), "cap of 12"),
        (dict(beta="Cl:2,Li:2"), "0 or 1 Cl"),
        (dict(beta="Cl:1,Li:2"), "exceed N_Cl"),
        (dict(alpha="Cl:2,Mg:1"), "species set"),
        (dict(alpha="Li:1", beta=None), "never be anchored"),
        (dict(alpha="Cl:3,Li:1"), "Excess charge invariant"),
        (dict(alpha="Cl:2,Li:2", counterion_n4="Mg:1"), "outside the species set"),
    ],
)
def test_request_domain_validation(kw, needle):
    args = dict(alpha=ALPHA, beta=BETA)
    args.update(kw)
    with pytest.raises(S.CageSaltInputError, match=re.escape(needle)):
        S.validate_request(request(**args), P.CATION_CHARGES, ["Na"])


def test_framework_cations_outside_s_are_rejected():
    with pytest.raises(S.CageSaltInputError, match="Framework cation Ca"):
        S.validate_request(request(), P.CATION_CHARGES, ["Ca"])


def test_salt_spec_parsing():
    assert S.parse_salt_spec("Cl:10, li:6,K:4") == {"Cl": 10, "Li": 6, "K": 4}
    for bad in ("Cl", "Cl:-1", "Cl:1,Cl:2", "Cl:1.5"):
        with pytest.raises(S.CageSaltInputError):
            S.parse_salt_spec(bad)
    r = S.CageSaltRequest.from_specs(alpha="Cl:2,Li:2", alpha1="Cl:4,K:4")
    assert r.cage("alpha-0") == {"Cl": 2, "Li": 2} and r.cage("alpha-1") == {"Cl": 4, "K": 4}
    assert r.n_cl() == 6 and r.is_active


# ----------------------------------------------------------------------------- S2 window table


def _sites(**kw):
    _, framework, sites, _ = P.load_lta_template(REAL, **kw)
    return framework, sites


def test_window_table_matches_loaded_sites_one_to_one():
    _, sites = _sites()
    table = S.window_table()
    assert len(table) == 46
    classes = Counter(S.site_class(sid) for sid in table)
    assert classes == {"4R@alpha-0": 12, "4R@alpha-1": 12, "6R@beta-0+alpha-0": 4, "6R@beta-0+alpha-1": 4,
                       "6R@beta-1+alpha-0": 4, "6R@beta-1+alpha-1": 4, "8R@alpha-0+alpha-1": 6}
    lattice = P.read_cif_minimal(REAL).lattice
    S.check_site_join(sites, np.asarray(lattice.matrix), CENTRES)


def test_window_table_join_fails_closed_with_nondefault_grouping():
    _, sites = _sites(variant_group_cutoffs={"4": 0.0, "6": 0.0, "8": 0.0})
    assert sum(len(v) for v in sites.values()) == 64
    lattice = P.read_cif_minimal(REAL).lattice
    with pytest.raises(S.CageSaltInputError):
        S.check_site_join(sites, np.asarray(lattice.matrix), CENTRES)


def test_8r_equidistance_check_fails_closed_on_a_strained_metric():
    _, sites = _sites()
    M = np.asarray(P.read_cif_minimal(REAL).lattice.matrix).copy()
    M[0] *= 1.02
    with pytest.raises(S.CageSaltInputError, match="equidistant"):
        S.check_site_join(sites, M, CENTRES)


def test_shifted_template_is_rejected_up_front(tmp_path: Path):
    shifted = tmp_path / "shifted.cif"
    shifted.write_text(REAL.read_text().replace(
        "SI1   Si  0.27815000  0.09183000  0.90817000", "SI1   Si  0.28815000  0.09183000  0.90817000"))
    with pytest.raises(ValueError, match="canonical Fm-3c primitive LTA template identity"):
        P.generate_lta_cation_structure_from_specs(shifted, cage_salt=request(), seed=1, **FW)


# ----------------------------------------------------------------------------- S2 charge


def test_cage_charge_sums_to_zero_for_random_placements():
    _, sites = _sites()
    req = request(alpha="Cl:10,Li:4,K:4", counterion_n4="Li:4")
    for seed in range(5):
        res = P.generate_lta_cation_structure_from_specs(
            REAL, cage_salt=req, seed=seed, _defer_cage_salt=True, **FW)
        q_exc, q_fw, _ = S.cage_charges_quarters(req, res.placed_ions, P.CATION_CHARGES)
        assert sum(q_exc.values()) == 0 and sum(q_fw.values()) == 0


def _brute_force_best(req, placed, sites_by_ring, charges):
    """Brute force over individual site assignments of the counterions (D2 §9 oracle 4)."""
    fw = [i for i in placed if i.role == "framework"]
    counter = [i for i in placed if i.role == "occlusion_counterion"]
    occupied = {i.site_id for i in fw}
    groups = sorted({(i.ring_class, charges[i.element]) for i in counter})
    n0 = Counter((i.ring_class, charges[i.element], S.site_class(i.site_id)) for i in counter)
    per_group = []
    for g in groups:
        n = sum(1 for i in counter if (i.ring_class, charges[i.element]) == g)
        empty = [s.site_id for s in sites_by_ring[g[0]] if s.site_id not in occupied]
        per_group.append([(g, combo) for combo in itertools.combinations(empty, n)])
    best = None
    for choice in itertools.product(*per_group):
        used = [sid for _, combo in choice for sid in combo]
        if len(set(used)) != len(used):
            continue
        ions = list(fw)
        for (fam, z), combo in choice:
            el = next(i.element for i in counter if (i.ring_class, charges[i.element]) == (fam, z))
            ions += [replace(counter[0], element=el, site_id=sid, ring_class=fam) for sid in combo]
        q_exc, q_fw, zmax = S.cage_charges_quarters(req, ions, charges)
        if any(S.charge_violations(q_exc, q_fw, zmax).values()):
            continue
        cls = Counter((fam, z, S.site_class(sid)) for (fam, z), combo in choice for sid in combo)
        moved = sum(max(0, n0[k] - cls.get(k, 0)) for k in n0)
        sumsq = sum((q_exc[c] + q_fw[c]) ** 2 for c in S.CAGE_IDS)
        vec = tuple(cls.get((g[0], g[1], c), 0) for g in groups for c in S.CLASS_ORDER if c[0] == g[0])
        key = (moved, sumsq, vec)
        best = key if best is None or key < best else best
    return best


def _vector_tuple(rb):
    out = []
    for g in sorted(rb.vector, key=lambda k: (k.split("R/z")[0], int(k.split("R/z")[1]))):
        fam = g.split("R/z")[0]
        out.extend(rb.vector[g].get(c, 0) for c in S.CLASS_ORDER if c[0] == fam)
    return tuple(out)


REBALANCE_CASES = [
    # one 4R group of one charge
    (FW, dict(alpha="Cl:10,Li:4,K:4", counterion_n4="Li:4"), range(3)),
    (FW, dict(alpha="Cl:11,Li:5,K:5", counterion_n4="Li:2"), range(3)),
    # two charge groups sharing the 6R family (review finding B1: needs a cross-group swap)
    (dict(n4="Na:4", n6="Na:14", n8="Na:6"),
     dict(alpha=None, beta=None, alpha0="Cl:3,Li:1", alpha1="Cl:2,Li:1", counterion_n6="Li:1,Sr:1"), range(40)),
    # two families, mixed charge
    (dict(n4="Na:4", n6="Na:14", n8="Na:6"),
     dict(alpha="Cl:4,Li:1,K:1", beta=None, counterion_n4="Li:2", counterion_n6="Sr:1"), range(10)),
]


@pytest.mark.parametrize("fw,spec,seeds", REBALANCE_CASES)
def test_rebalancing_matches_brute_force(fw, spec, seeds):
    req = request(**spec)
    S.validate_request(req, P.CATION_CHARGES, ["Na"])
    for seed in seeds:
        res = P.generate_lta_cation_structure_from_specs(
            REAL, cage_salt=req, seed=seed, _defer_cage_salt=True, _skip_cage_salt_validation=True, **fw)
        rb = S.rebalance(req, res.placed_ions, res.physical_sites_by_ring, P.CATION_CHARGES,
                         S.rebalance_stream(seed)[0])
        brute = _brute_force_best(req, res.placed_ions, res.physical_sites_by_ring, P.CATION_CHARGES)
        assert rb.feasible == (brute is not None), seed
        if brute is None:
            assert rb.best_violation_quarters > 0 and rb.blocking_sites
            continue
        assert (*rb.objective, _vector_tuple(rb)) == brute, seed
        # the realization reaches the chosen vector with one ion per site
        assert len({i.site_id for i in rb.placed}) == len(rb.placed)
        realized = Counter((i.ring_class, P.CATION_CHARGES[i.element], S.site_class(i.site_id))
                           for i in rb.placed if i.role == "occlusion_counterion")
        for g, row in rb.vector.items():
            fam, z = g.split("R/z")
            for cls, n in row.items():
                assert realized.get((fam, int(z), cls), 0) == n
        q_exc, q_fw, zmax = S.cage_charges_quarters(req, rb.placed, P.CATION_CHARGES)
        assert not any(S.charge_violations(q_exc, q_fw, zmax).values())
        assert sum(q_exc[c] + q_fw[c] for c in S.CAGE_IDS) == 0


def test_upfront_charge_check_rejects_unbalanceable_requests():
    # all salt charge sits in alpha-0 (-10 e); capacity suffices (6 8R + 4 4R sites) but 8R
    # counterions split their charge between both alpha-cages, so no vector is feasible
    req = S.CageSaltRequest.from_specs(alpha0="Cl:12,Li:2", alpha1="Cl:0", counterion_n8="Li:6",
                                       counterion_n4="Li:4")
    S.validate_request(req, P.CATION_CHARGES, ["Na"])
    with pytest.raises(S.CageSaltInputError, match=r"rejected up front") as exc:
        S.upfront_charge_check(req, P.CATION_CHARGES)
    assert "best excess 0.00" not in str(exc.value)


def test_upfront_rejection_names_the_cages():
    req = S.CageSaltRequest.from_specs(alpha0="Cl:12,Li:2", alpha1="Cl:0", counterion_n8="Li:6",
                                       counterion_n4="Li:4")
    with pytest.raises(S.CageSaltInputError, match=r"in alpha-0"):
        S.upfront_charge_check(req, P.CATION_CHARGES)


@pytest.mark.parametrize("cages,needle", [
    ({"alpha-0": {"Cl": 2.7, "Li": 1, "K": 1.9}}, "nonnegative integer"),
    ({"alpha-0": {"Cl": 2, "Li": 3, "K": -1}}, "nonnegative integer"),
    ({"alpha-0": {"Cl": True, "Li": 1}}, "nonnegative integer"),
])
def test_api_requests_must_carry_nonnegative_integer_counts(cages, needle):
    with pytest.raises(S.CageSaltInputError, match=needle):
        S.validate_request(S.CageSaltRequest(cages=cages), P.CATION_CHARGES, ["Na"])
    with pytest.raises(S.CageSaltInputError, match=needle):
        generate(S.CageSaltRequest(cages=cages))


def test_unknown_cage_or_ring_keys_are_rejected_at_construction():
    with pytest.raises(S.CageSaltInputError, match="Unknown cage"):
        S.CageSaltRequest(cages={"gamma-0": {"Cl": 1, "Li": 1}})
    with pytest.raises(S.CageSaltInputError, match="Unknown ring family"):
        S.CageSaltRequest(cages={}, counterions={"5": {"Li": 1}})


def test_v1_compares_ring_counts_with_the_input_allocations(salt_result, ref_cell):
    """D1 §4: an equal-charge species swap on a ring site fails V1 (and only V1)."""
    M, F, atoms, req, kw = _verify_inputs(salt_result, ref_cell[0])
    exp = {("framework", "4", "Na"): 2, ("framework", "6", "Na"): 16, ("framework", "8", "Na"): 6}
    assert S.verify(M, F, atoms, req, **{**kw, "expected_ring_counts": exp}).passed
    swapped = dict(exp)
    swapped[("framework", "8", "Na")] = 5
    swapped[("framework", "8", "K")] = 1
    assert S.verify(M, F, atoms, req, **{**kw, "expected_ring_counts": swapped}).failing == ["V1_inventory"]


def test_cl_only_cages_round_trip_through_the_independent_verifier(tmp_path):
    """x_M falls back to all cations when the cell has no cage-salt cation (D1 §6)."""
    r = P.generate_lta_cation_structure_from_specs(
        REAL, cage_salt=request(alpha="Cl:4", beta=None, counterion_n4="Li:8"), seed=6, **FW)
    path = tmp_path / "clonly.cif"
    r.write(path)
    assert independent_verify(path) == {"V1": True, "V2": True, "V3": True, "V4": True, "V5": True, "V6": True}
    # a Cl-Cl clash the fallback must catch: move one Cl 2.0 A from another (r_cut ~2.97 A)
    side = json.loads(Path(str(path) + ".salt.json").read_text())
    cl = [k for k, row in enumerate(side["atoms"]) if row["element"] == "Cl"]
    M = np.asarray(r.structure.lattice.matrix)
    a, b = np.asarray(side["atoms"][cl[0]]["frac"]), np.asarray(side["atoms"][cl[1]]["frac"])
    side["atoms"][cl[1]]["frac"] = list(a + np.array([2.0, 0, 0]) @ np.linalg.inv(M))
    from pymatgen.core import Structure
    st = Structure(r.structure.lattice, [row["element"] for row in side["atoms"]],
                   [row["frac"] for row in side["atoms"]])
    bad = tmp_path / "clash.cif"
    P.write_structure_by_suffix(st, bad)
    side["structure_sha256"] = hashlib.sha256(bad.read_bytes()).hexdigest()
    Path(str(bad) + ".salt.json").write_text(json.dumps(side))
    assert independent_verify(bad)["V2"] is False


# ----------------------------------------------------------------------------- S2 volume (oracle 1)


@pytest.fixture(scope="module")
def ref_cell():
    res = P.generate_lta_cation_structure_from_specs(REAL, cage_salt=request(), seed=3, _defer_cage_salt=True, **FW)
    M = np.asarray(res.structure.lattice.matrix, float)
    o = np.array([a.frac for a in res.framework_atoms if a.element == "O"], float)
    return M, o


def test_v_eff_is_v_over_4_at_zero_and_symmetric(ref_cell):
    M, o = ref_cell
    vol = abs(np.linalg.det(M))
    for c in S.CAGE_IDS:
        assert S.v_eff(M, o, CENTRES, c, 0.0) == pytest.approx(vol / 4, rel=1e-12)
    for d in (1.0, 1.62, 2.2):
        assert S.v_eff(M, o, CENTRES, "alpha-0", d) == S.v_eff(M, o, CENTRES, "alpha-1", d)
        assert S.v_eff(M, o, CENTRES, "beta-0", d) == S.v_eff(M, o, CENTRES, "beta-1", d)


def test_v_eff_scale_identity(ref_cell):
    M, o = ref_cell
    s = 0.97 ** (1 / 3)
    S._FIELD_CACHE.clear()
    direct = S.v_eff(M * s, o, CENTRES, "alpha-0", 1.62)
    via = S.v_eff(M, o, CENTRES, "alpha-0", 1.62, scale=s)
    assert via == pytest.approx(direct, rel=1e-9)
    assert S.v_eff(M, o, CENTRES, "alpha-0", 1.62 * s, scale=s) == pytest.approx(
        s ** 3 * S.v_eff(M, o, CENTRES, "alpha-0", 1.62), rel=1e-12)


def test_v_eff_cold_and_warm_cache_identical(ref_cell):
    M, o = ref_cell
    S._FIELD_CACHE.clear()
    cold = [S.v_eff(M, o, CENTRES, c, 1.6) for c in S.CAGE_IDS]
    warm = [S.v_eff(M, o, CENTRES, c, 1.6) for c in S.CAGE_IDS]
    assert cold == warm


def test_v_eff_agrees_with_monte_carlo_at_d_wall(ref_cell):
    from scipy.spatial import cKDTree

    M, o = ref_cell
    d_wall = S.cage_floors(request(), M, o, CENTRES, 1.0)["alpha-0"].d_wall
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3)), float)
    o_img = (o[None] + shifts[:, None]).reshape(-1, 3) @ M
    names = list(S.CAGE_IDS)
    c_img = (np.array([CENTRES[c] for c in names])[None] + shifts[:, None]).reshape(-1, 3) @ M
    labels = np.tile(np.arange(4), len(shifts))
    to, tc = cKDTree(o_img), cKDTree(c_img)
    rng = np.random.default_rng(716)
    hits, total = 0, 0
    for _ in range(10):
        X = rng.random((1_000_000, 3)) @ M
        d_o, _ = to.query(X)
        _, ic = tc.query(X)
        hits += int(np.count_nonzero((labels[ic] == 0) & (d_o > d_wall)))
        total += len(X)
    mc = abs(np.linalg.det(M)) * hits / total
    assert S.v_eff(M, o, CENTRES, "alpha-0", d_wall) == pytest.approx(mc, rel=5e-3)


def test_floors_follow_rule_b(ref_cell):
    M, o = ref_cell
    fl = S.cage_floors(request(), M, o, CENTRES, 0.97 ** (1 / 3))["alpha-0"]
    v_pair = S.v_pair_800k(7 / 12)
    assert fl.n_melt == pytest.approx(fl.v_eff / v_pair)
    phi = min(1.0, 12 / fl.n_melt)
    for m in ("Li", "K"):
        assert fl.c_by_species[m] == pytest.approx(S.T_BULK[m] * phi)
        assert fl.floor_by_species[m] == max(1, int(np.ceil(0.4 * S.T_BULK[m] * phi - 1e-9)))
    assert S.v_pair_800k(0.592) == pytest.approx(57.35, abs=0.01)
    assert S.f_v(26) == pytest.approx(0.97) and S.f_v(10) == pytest.approx(1 - 0.03 * 5 / 10.5)


def test_streams_are_tagged_and_seed_none_uses_os_entropy():
    a, tag = S.start_stream(7, 2)
    b, _ = S.start_stream(7, 2)
    assert tag == "v7.16-cage-salt:7:start2" and a.random() == b.random()
    assert S.rebalance_stream(7)[1] == "v7.16-charge-rebalance:7"
    n1, t1 = S.start_stream(None, 1)
    n2, _ = S.start_stream(None, 1)
    assert t1 is None and n1.random() != n2.random()


# ----------------------------------------------------------------------------- S3 build


def _problem(rule):
    req = request()
    res = P.generate_lta_cation_structure_from_specs(REAL, cage_salt=req, seed=3, _defer_cage_salt=True, **FW)
    M = np.asarray(res.structure.lattice.matrix, float)
    o = np.array([a.frac for a in res.framework_atoms if a.element == "O"], float)
    floors = S.cage_floors(req, M, o, CENTRES, 1.0, rule)
    salt = [B.SaltIon(f"S-{c}-{e}-{k}", e, c) for c in S.CAGE_IDS for e, n in req.cage(c).items()
            for k in range(1, n + 1)]
    pairs = S.PairTable(S.cl_weights([x.element for x in salt if x.element != "Cl"], []))
    return B.BuildProblem(M, salt, [a.element for a in res.framework_atoms],
                          np.array([a.frac for a in res.framework_atoms]), res.placed_ions,
                          res.physical_sites_by_ring, pairs, floors, CENTRES, req, P.CATION_CHARGES, rule)


@pytest.mark.parametrize("rule", ["B", "A"])
@pytest.mark.parametrize("W,Wdz", [(10.0, 0.0), (1000.0, 10.0)])
def test_objective_gradient_matches_finite_differences(rule, W, Wdz):
    prob = _problem(rule)
    model = B._Model(prob, prob.ring_ions)
    if rule == "A":
        assert model.groups
    rng = np.random.default_rng(1)
    x = (B._insert(model, rng) + rng.normal(scale=0.3, size=(len(prob.salt), 3))).ravel()
    E, G = model.energy(x, W, Wdz)
    h = 1e-6
    idx = rng.choice(len(x), size=40, replace=False)
    num = np.array([(model.energy(x + h * np.eye(len(x))[i], W, Wdz)[0]
                     - model.energy(x - h * np.eye(len(x))[i], W, Wdz)[0]) / (2 * h) for i in idx])
    assert E > 0
    assert np.allclose(G[idx], num, rtol=1e-5, atol=1e-5 * max(1.0, np.abs(G).max()))


# ----------------------------------------------------------------------------- S2/S3 verification (oracle 3)


def _verify_inputs(result, ref_matrix):
    _, meta = P.assemble_structure(result.structure.lattice, result.framework_atoms, result.placed_ions,
                                   (), result.cage_salt_ions)
    atoms = P._atoms_for_verify(meta)
    F = S.canonical_frac(np.asarray(result.structure.frac_coords))
    M = np.asarray(result.structure.lattice.matrix)
    req = S.CageSaltRequest(cages=result.cage_salt_record.request["cages"],
                            counterions=result.cage_salt_record.request["counterions"])
    salt_cat = [i.element for i in result.cage_salt_ions if i.element != "Cl"]
    pairs = S.PairTable(S.cl_weights(salt_cat, [i.element for i in result.placed_ions]))
    o = np.array([a.frac for a in result.framework_atoms if a.element == "O"], float)
    s_rel = result.cage_salt_record.volume["scale_from_V_ref"]
    floors = S.cage_floors(req, ref_matrix, o, CENTRES, s_rel)
    expected = Counter((i.role, i.ring_class, i.element) for i in result.placed_ions)
    site_variants = {s.site_id: [v.frac for v in s.variants]
                     for r in P.RING_CLASSES for s in result.physical_sites_by_ring[r]}
    kw = dict(pairs=pairs, floors=floors, centres=CENTRES, charges=dict(P.CATION_CHARGES),
              expected_ring_counts=dict(expected), site_variants=site_variants, ring_ions=list(result.placed_ions))
    return M, F, atoms, req, kw


def test_manufactured_feasible_configuration_passes(salt_result, ref_cell):
    M, F, atoms, req, kw = _verify_inputs(salt_result, ref_cell[0])
    assert S.verify(M, F, atoms, req, **kw).passed


def test_each_single_violation_is_rejected_by_exactly_that_check(salt_result, ref_cell, monkeypatch):
    M, F, atoms, req, kw = _verify_inputs(salt_result, ref_cell[0])
    base = S.verify(M, F, atoms, req, **kw)
    # V1: a ring cation off its template variant by 1e-6 fractional
    ring_i = next(i for i, a in enumerate(atoms) if a.role == "framework")
    F1 = F.copy()
    F1[ring_i] += 1e-6
    assert S.verify(M, F1, atoms, req, **kw).failing == ["V1_inventory"]
    # V2: a margin just above the achieved minimum slack
    tight = min(base.detail["V2"]["min_slack"].values()) + 1e-4
    assert tight < base.detail["V3"]["min_distance_to_r_min1"] and tight < base.detail["V1"]["min_membership_margin"]
    assert S.verify(M, F, atoms, req, **{**kw, "tau": tight}).failing == ["V2_hard_minima"]
    # V3: an unreachable rule-B floor
    floors = {c: replace(f, floor_by_species={m: 6 for m in f.floor_by_species}) for c, f in kw["floors"].items()}
    assert S.verify(M, F, atoms, req, **{**kw, "floors": floors}).failing == ["V3_coordination"]
    # V5: one framework 8R cation relabelled as a counterion (weights 1/2 + 1/2 keep V4)
    k = next(n for n, r in enumerate(kw["ring_ions"]) if r.ring_class == "8")
    ring2 = list(kw["ring_ions"])
    ring2[k] = replace(ring2[k], role="occlusion_counterion")
    exp = Counter((r.role, r.ring_class, r.element) for r in ring2)
    atoms2 = list(atoms)
    j = next(i for i, a in enumerate(atoms) if a.ion_id == f"R{ring2[k].site_id}")
    atoms2[j] = replace(atoms2[j], role="occlusion_counterion")
    out = S.verify(M, F, atoms2, req, **{**kw, "ring_ions": ring2, "expected_ring_counts": dict(exp)})
    assert out.failing == ["V5_invariants"]
    # V6: the single-image requirement raised above the half width
    monkeypatch.setattr(S, "max_r_min1", lambda: 100.0)
    assert S.verify(M, F, atoms, req, **kw).failing == ["V6_single_image"]


def test_salt_ions_keep_t_min_from_framework_si_al(salt_result, ref_cell, monkeypatch):
    """D1 r8 S-13: every cage_salt ion is >= 3.0 A from every Si/Al, and V2 enforces it alone."""
    M, F, atoms, req, kw = _verify_inputs(salt_result, ref_cell[0])
    salt = [i for i, a in enumerate(atoms) if a.role == "cage_salt"]
    tet = [i for i, a in enumerate(atoms) if a.element in ("Si", "Al")]
    d_min = min(float(np.linalg.norm(S.min_image(F[j] - F[i], M))) for i in salt for j in tet)
    assert S.T_MIN == 3.0 and d_min >= 3.0 + S.TAU_W
    base = S.verify(M, F, atoms, req, **kw)
    assert base.detail["V2"]["min_slack"]["salt-T"] == pytest.approx(d_min - 3.0)
    monkeypatch.setattr(S, "T_MIN", d_min + 1e-4)
    assert S.verify(M, F, atoms, req, **kw).failing == ["V2_hard_minima"]


def test_v4_cage_charge_violation_is_rejected_alone():
    # alpha-0 carries +3 e and alpha-1 -3 e of salt (framework q_fw is at most 1 e);
    # Q_exc = N_Cl still holds, so V5 passes
    req = S.CageSaltRequest.from_specs(alpha0="Cl:1,Li:4", alpha1="Cl:4,Li:1", beta="Cl:0")
    res = P.generate_lta_cation_structure_from_specs(REAL, cage_salt=req, seed=2, _defer_cage_salt=True,
                                                     _skip_cage_salt_validation=True, **FW)
    prob = _problem("B")
    o = np.array([a.frac for a in res.framework_atoms if a.element == "O"], float)
    M = np.asarray(res.structure.lattice.matrix, float)
    floors = S.cage_floors(req, M, o, CENTRES, 1.0)
    salt = [B.SaltIon(f"S-{c}-{e}-{k}", e, c) for c in S.CAGE_IDS for e, n in req.cage(c).items()
            for k in range(1, n + 1)]
    prob = replace(prob, matrix=M, salt=salt, ring_ions=res.placed_ions, floors=floors, request=req,
                   pairs=S.PairTable(S.cl_weights(["Li"] * 5, [])))
    fw_atoms = [S.Atom(a.element, "framework_atom", f"F{k}") for k, a in enumerate(res.framework_atoms)]
    expected = Counter((i.role, i.ring_class, i.element) for i in res.placed_ions)
    site_variants = {s.site_id: [v.frac for v in s.variants]
                     for r in P.RING_CLASSES for s in res.physical_sites_by_ring[r]}

    def run(fr, ring):
        atoms = fw_atoms + [S.Atom(r.element, r.role, f"R{r.site_id}", site_id=r.site_id) for r in ring] + [
            S.Atom(x.element, "cage_salt", x.ion_id, cage_id=x.cage_id) for x in salt]
        Fa = np.vstack([prob.framework_fracs, np.array([r.frac for r in ring]), fr])
        return S.verify(M, Fa, atoms, req, pairs=prob.pairs, floors=floors, centres=CENTRES,
                        charges=P.CATION_CHARGES, expected_ring_counts=dict(expected),
                        site_variants=site_variants, ring_ions=ring)

    def verifier(fr, ring):
        v = run(fr, ring)
        rest = [f for f in v.failing if f != "V4_cage_charge"]
        return S.Verification(not rest, v.checks, v.detail, rest)

    out = B.build(prob, verifier, seed=4, k_max=5)
    assert out.success
    assert run(out.fracs, out.ring_ions).failing == ["V4_cage_charge"]


def test_rule_a_pooling():
    fl = {c: S.CageFloors(c, 4, 0.5, False, 0.0, False, 1.6, 800.0, (700.0, 900.0), 14.0, (12.0, 16.0), 1.0,
                          {"Li": 4.3, "Cs": 5.5}, {"Li": 2, "Cs": 3}, 0.7, "A") for c in S.ALPHA_IDS}
    ok, detail = S.rule_a_check({("alpha-0", "Li"): [4, 3, 3], ("alpha-1", "Li"): [2], ("alpha-0", "Cs"): [1],
                                 ("alpha-1", "Cs"): [9]}, fl)
    assert detail["alpha-0/Li"]["required"] == pytest.approx(0.7 * 4.3 * 3)
    assert "fewer than 3" in detail["pooled/Li"]["note"] and "fewer than 3" in detail["pooled/Cs"]["note"]
    assert ok
    ok, detail = S.rule_a_check({("alpha-0", "Li"): [1, 1, 1], ("alpha-1", "Li"): [2, 2]}, fl)
    assert not ok and detail["alpha-0/Li"]["ok"] is False


# ----------------------------------------------------------------------------- S4 resolver


def test_resolver_precedence_and_classification():
    salt = request()
    with pytest.raises(ValueError, match="sparse-v715"):
        P.resolve_occlusion_mode("cage-salt", alpha_cl=1)
    with pytest.raises(ValueError, match="not valid with --occlusion-model sparse-v715"):
        P.resolve_occlusion_mode("sparse-v715", alpha_cl=1, cage_salt=salt)
    with pytest.raises(ValueError, match="counterions were given"):
        P.resolve_occlusion_mode("cage-salt", cage_salt=S.CageSaltRequest.from_specs(counterion_n4="Li:1"))
    assert isinstance(P.resolve_occlusion_mode("sparse-v715"), P.OcclusionDry)
    assert isinstance(P.resolve_occlusion_mode("cage-salt"), P.OcclusionDry)
    assert isinstance(P.resolve_occlusion_mode("cage-salt", cage_salt=S.CageSaltRequest.from_specs(alpha="Cl:0")),
                      P.OcclusionDry)
    leg = P.resolve_occlusion_mode("sparse-v715", alpha_cl=1, occlusion_n4="Cs:1")
    assert leg == P.OcclusionLegacySparse(1, 0, "Cs:1", None, None, 1.5)
    assert isinstance(P.resolve_occlusion_mode(cage_salt=salt), P.OcclusionCageSalt)
    with pytest.raises(ValueError):
        P.resolve_occlusion_mode("other")


def test_volume_factor_resolution():
    assert P.resolve_volume_factor(None, P.OcclusionDry()) == 1.0
    assert P.resolve_volume_factor(None, P.OcclusionCageSalt(request())) == pytest.approx(0.97)
    assert P.resolve_volume_factor(0.99, P.OcclusionCageSalt(request())) == 0.99
    with pytest.raises(ValueError):
        P.resolve_volume_factor(0.0, P.OcclusionDry())


def test_dry_and_legacy_imports_load_no_cage_salt_module():
    code = (
        "import sys\n"
        "import lta_cation_population_grouped as P, generate_smooth_lta_cations_grouped\n"
        "P.generate_lta_cation_structure_from_specs('tests/real_lta_primitive.cif', n4='Na:2', n6='Na:16', n8='Na:6', seed=0)\n"
        "P.generate_lta_cation_structure_from_specs('tests/real_lta_primitive.cif', n4='Na:2', n6='Na:16', n8='Na:6', "
        "seed=0, occlusion_model='sparse-v715', alpha_cl=1, occlusion_n4='Cs:1')\n"
        "bad = sorted(m for m in sys.modules if m.startswith('lta_cage_salt'))\n"
        "print(bad)\n"
    )
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, check=True)
    # pymatgen itself imports scipy.optimize on this stack, so the "no SciPy" half of the
    # D3 claim cannot be observed here; only the cage-salt modules are asserted absent.
    assert out.stdout.strip() == "[]"


# ----------------------------------------------------------------------------- S4 integration


def test_cage_salt_result_derived_values(salt_result):
    r = salt_result
    assert r.cage_salt_record.status == "passed" and r.cage_salt_pending is None
    assert r.chloride_count == 26 and r.alpha_chloride_count == 24 and r.beta_chloride_count == 2
    assert r.cation_counts["Li"] == 16 and r.cation_counts["K"] == 10 and r.cation_counts["Na"] == 24
    assert r.occlusion_positive_charge == r.chloride_count
    assert r.volume_factor == pytest.approx(0.97)
    assert r.target_volume_ang3 == pytest.approx(r.reference_volume_ang3 * 0.97)
    assert float(r.structure.volume) == pytest.approx(r.target_volume_ang3)
    ids = [i.ion_id for i in r.cage_salt_ions]
    assert len(ids) == len(set(ids)) == 52
    assert "S-alpha-0-Cl-12" in ids and "S-beta-1-Li-1" in ids
    text = P.make_population_report(r)
    assert "Cage-salt occlusion (v7.16" in text and "V5_invariants=PASS" in text


def test_counterion_rebalancing_end_to_end(counterion_result):
    r = counterion_result
    q = r.cage_salt_record.verification["detail"]["V4"]
    for c in S.CAGE_IDS:
        assert abs(q["q_exc"][c] + q["q_fw"][c]) <= q["q_max"][c]
    assert r.occlusion_counterion_counts["Li"] == 4
    assert r.cage_salt_record.rebalance["feasible"]


def test_pending_intermediate_refuses_to_write_and_reports_pending(tmp_path):
    res = P.generate_lta_cation_structure_from_specs(REAL, cage_salt=request(), seed=3, _defer_cage_salt=True, **FW)
    assert res.cage_salt_pending is not None and res.volume_factor == 1.0
    with pytest.raises(ValueError, match="pending"):
        res.write(tmp_path / "x.cif")
    assert "cage-salt pending" in P.make_population_report(res)


def test_significant_figures_floor(salt_result, tmp_path):
    with pytest.raises(ValueError, match="significant_figures >= 8"):
        salt_result.write(tmp_path / "x.cif", significant_figures=6)
    assert not (tmp_path / "x.cif").exists()


def test_oracle5_step3_snapshot_independent_of_salt_and_f_v():
    def snap(req, vf=None):
        res = P.generate_lta_cation_structure_from_specs(
            REAL, cage_salt=req, seed=11, volume_factor=vf, _defer_cage_salt=True, **FW)
        ring = [(i.element, i.role, i.site_id, i.variant_id, tuple(i.frac)) for i in res.placed_ions]
        return ring, np.asarray(res.structure.lattice.matrix).tobytes()

    base = snap(request(alpha="Cl:5,Li:2,K:2", counterion_n4="Li:2"))
    assert snap(request(alpha="Cl:12,Li:6,K:5", counterion_n4="Li:2")) == base
    assert snap(request(alpha="Cl:5,Li:2,K:2", counterion_n4="Li:2"), vf=0.9) == base


def test_oracle6_determinism_same_seed_same_bytes(tmp_path):
    a, b = tmp_path / "a.cif", tmp_path / "b.cif"
    generate(request(alpha="Cl:6,Li:3,K:3"), seed=21).write(a)
    generate(request(alpha="Cl:6,Li:3,K:3"), seed=21).write(b)
    assert a.read_bytes() == b.read_bytes()
    ja = json.loads(Path(str(a) + ".salt.json").read_text())
    jb = json.loads(Path(str(b) + ".salt.json").read_text())
    assert ja["atoms"] == jb["atoms"] and ja["structure_sha256"] == jb["structure_sha256"]


def test_oracle7_single_image_precondition_fails_closed():
    with pytest.raises(S.CageSaltConstructionError, match="Single-image precondition") as exc:
        generate(request(alpha="Cl:4,Li:2,K:2"), volume_factor=0.3)
    assert exc.value.record.status == "failed"


def test_forced_ring_conflict_is_resolved_by_relocation_end_to_end(monkeypatch, tmp_path):
    """D2 §6.4 through finalize: the first relaxation is followed by a forced Cl/ring clash."""
    real_relax = B._relax
    calls = {"n": 0}

    def clashing_relax(model, X, log):
        X = real_relax(model, X, log)
        calls["n"] += 1
        if calls["n"] == 1:
            ring = next(r for r in model.ring if r.ring_class == "4" and r.role == "framework")
            cage = S.site_class(ring.site_id).split("@")[1]
            i = next(k for k, s_ in enumerate(model.p.salt) if s_.element == "Cl" and s_.cage_id == cage)
            X = X.copy()
            X[i] = np.asarray(ring.frac, float) @ model.M + np.array([0.3, 0.0, 0.0])
        return X

    monkeypatch.setattr(B, "_relax", clashing_relax)
    # the move must be logged on start 1; a relocated start may still fail later checks
    # (then a later start wins and the record keeps no move), so try a few seeds
    for seed in range(8, 20):
        calls["n"] = 0
        r = generate(request(alpha="Cl:6,Li:3,K:3"), seed=seed)
        rec = r.cage_salt_record
        assert rec.starts[0]["rounds"][0].get("relocations"), seed
        if rec.relocations:
            break
    assert rec.status == "passed" and rec.relocations and rec.verification["start_used"] == 1
    move = rec.relocations[0]
    assert move["from_site"] != move["to_site"] and move["to_site"].startswith("4R")
    assert move["to_ion_id"] == "R" + move["to_site"]
    # D1 §4: each move reports its charge effect
    assert set(move["q_before_e"]) == set(move["q_after_e"]) == set(S.CAGE_IDS)
    assert sum(move["q_after_e"].values()) == pytest.approx(0.0)
    assert any(i.site_id == move["to_site"] for i in r.placed_ions)
    assert not any(i.site_id == move["from_site"] for i in r.placed_ions)
    path = tmp_path / "reloc.cif"
    r.write(path)
    assert independent_verify(path) == {"V1": True, "V2": True, "V3": True, "V4": True, "V5": True, "V6": True}


def test_relocation_picks_the_largest_min_slack_resolving_candidate():
    """D2 §6.4 choice rule, recomputed independently of B._relocate, over several configurations."""
    prob = _problem("B")
    model = B._Model(prob, prob.ring_ions)
    ring = list(prob.ring_ions)
    k = next(n for n, r in enumerate(ring) if r.ring_class == "4" and r.role == "framework")
    target = ring[k]
    cage = S.site_class(target.site_id).split("@")[1]
    i = next(n for n, x in enumerate(prob.salt) if x.element == "Cl" and x.cage_id == cage)
    moved = 0
    for seed in range(10):
        X = B._insert(model, np.random.default_rng(seed))
        X[i] = np.asarray(target.frac, float) @ model.M + np.array([0.3, 0.0, 0.0])
        _, moves = B._relocate(model, X, prob, ring)
        # independent candidate scan: empty 4R sites and the other variants of the own site
        occupied = {r.site_id for r in ring}
        best = None
        for site in prob.sites_by_ring["4"]:
            if site.site_id != target.site_id and site.site_id in occupied:
                continue
            for vid, var in enumerate(site.variants):
                if site.site_id == target.site_id and vid == target.variant_id:
                    continue
                y = np.asarray(var.frac, float) @ model.M
                slack = min(float(np.linalg.norm(S.min_image((X[j] - y) @ np.linalg.inv(model.M), model.M)))
                            - (prob.pairs.r_cut(x.element, target.element) + S.DELTA_B)
                            for j, x in enumerate(prob.salt))
                if slack < 0:
                    continue
                trial = list(ring)
                trial[k] = replace(target, frac=np.asarray(var.frac, float), site_id=site.site_id, variant_id=vid)
                if any(S.charge_violations(*S.cage_charges_quarters(prob.request, trial, P.CATION_CHARGES)).values()):
                    continue
                key = (-slack, site.site_id, vid)
                best = key if best is None or key < best else best
        mine = [m for m in moves if m["from_site"] == target.site_id]
        if best is None:
            assert not mine, seed  # no resolving candidate: the cation stays
            continue
        moved += 1
        assert (mine[0]["to_site"], mine[0]["to_variant"]) == (best[1], best[2]), seed
        assert mine[0]["min_slack_after"] == pytest.approx(-best[0])
    assert moved >= 1


def test_rebalancing_realization_follows_d2_departure_and_draw_order():
    """Replay R_q against an independent realization of D2 §5 (two-phase per family)."""
    # a two-species group (Li, K both z = 1) makes the departure order observable
    cases = [(FW, dict(alpha="Cl:10,Li:4,K:4", counterion_n4="Li:2,K:2"), range(60)),
             (FW, dict(alpha="Cl:10,Li:4,K:4", counterion_n4="Li:4"), range(12)),
             (dict(n4="Na:4", n6="Na:14", n8="Na:6"),
              dict(alpha=None, beta=None, alpha0="Cl:3,Li:1", alpha1="Cl:2,Li:1", counterion_n6="Li:1,Sr:1"),
              range(12))]
    z = P.CATION_CHARGES
    for fw, spec, seeds in cases:
        req = request(**spec)
        for seed in seeds:
            res = P.generate_lta_cation_structure_from_specs(
                REAL, cage_salt=req, seed=seed, _defer_cage_salt=True, _skip_cage_salt_validation=True, **fw)
            rb = S.rebalance(req, res.placed_ions, res.physical_sites_by_ring, z, S.rebalance_stream(seed)[0])
            if not rb.feasible:
                continue
            sites = {s_.site_id: s_ for r in res.physical_sites_by_ring.values() for s_ in r}
            placed = list(res.placed_ions)
            occupied = {i.site_id for i in placed}
            rng = S.rebalance_stream(seed)[0]
            for fam in P.RING_CLASSES:
                classes = [c for c in S.CLASS_ORDER if c[0] == fam]
                groups = sorted({(fam, int(g.split("R/z")[1])) for g in rb.vector if g.startswith(fam + "R")})
                leave, need = {}, {}
                for g in groups:
                    tgt = rb.vector[f"{fam}R/z{g[1]}"]
                    mine = [n for n, i in enumerate(placed)
                            if i.role == "occlusion_counterion" and i.ring_class == fam and z[i.element] == g[1]]
                    leave[g], need[g] = [], []
                    for c in classes:
                        cur = sorted((n for n in mine if S.site_class(placed[n].site_id) == c),
                                     key=lambda n: (placed[n].site_id, placed[n].element))
                        leave[g] += cur[:max(0, len(cur) - tgt.get(c, 0))]
                        need[g] += [c] * max(0, tgt.get(c, 0) - len(cur))
                for g in groups:
                    for n in leave[g]:
                        occupied.discard(placed[n].site_id)
                for g in groups:
                    for n, c in zip(leave[g], need[g]):
                        empty = sorted(sid for sid in sites if S.site_class(sid) == c and sid not in occupied)
                        sid = empty[int(rng.integers(len(empty)))]
                        vid = int(rng.integers(len(sites[sid].variants)))
                        occupied.add(sid)
                        placed[n] = replace(placed[n], site_id=sid, variant_id=sites[sid].variants[vid].variant_id)
            assert [(i.element, i.site_id, i.variant_id) for i in rb.placed] == \
                [(i.element, i.site_id, i.variant_id) for i in placed], seed


def test_internal_errors_fail_closed_with_the_record(monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("synthetic")

    monkeypatch.setattr(S, "rebalance", boom)
    with pytest.raises(S.CageSaltConstructionError, match="internal construction error") as exc:
        generate(request(alpha="Cl:4,Li:2,K:2"))
    assert exc.value.record.status == "failed"


def test_seed_level_rebalancing_failure_writes_the_failure_record(tmp_path):
    args = ["--n4", "Na:4", "--n6", "Na:14", "--n8", "Na:6", "--alpha0-salt", "Cl:3,Li:1",
            "--alpha1-salt", "Cl:2,Li:1", "--salt-counterion-n6", "Li:1,Sr:1"]
    req = request(alpha=None, beta=None, alpha0="Cl:3,Li:1", alpha1="Cl:2,Li:1", counterion_n6="Li:1,Sr:1")
    seed = None
    for s_ in range(40):
        res = P.generate_lta_cation_structure_from_specs(
            REAL, cage_salt=req, seed=s_, _defer_cage_salt=True, n4="Na:4", n6="Na:14", n8="Na:6")
        if not S.rebalance(req, res.placed_ions, res.physical_sites_by_ring, P.CATION_CHARGES,
                           S.rebalance_stream(s_)[0]).feasible:
            seed = s_
            break
    assert seed is not None
    out = _cli("lta_cation_population_grouped.py", tmp_path / "f.cif", *args, "--seed", str(seed))
    assert out.returncode == 1, out.stdout + out.stderr
    side = json.loads((tmp_path / "f.cif.salt.json").read_text())
    assert side["status"] == "failed" and side["rebalance"]["feasible"] is False
    assert side["rebalance"]["best_violation_e"] > 0 and side["rebalance"]["blocking_sites"]
    assert not (tmp_path / "f.cif").exists()


# ----------------------------------------------------------------------------- independent verifier (oracle 2)


def _csv_rows():
    rows = {}
    with open(ROOT / "doc" / "evidence" / "molten_salt_pair_distances.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            rows[(r["ion_i"], r["ion_j"], r["environment"])] = (float(r["r_cut_A"]), float(r["r_min1_A"]))
            rows[(r["ion_j"], r["ion_i"], r["environment"])] = (float(r["r_cut_A"]), float(r["r_min1_A"]))
    return rows


def _read_structure(path: Path):
    from pymatgen.core import Structure
    from pymatgen.io.cif import CifParser
    from pymatgen.io.lammps.data import LammpsData

    name = path.name.lower()
    if name.endswith(".cif"):
        return CifParser(str(path), frac_tolerance=0).parse_structures(primitive=False)[0]
    if name.endswith(".data"):
        return LammpsData.from_file(str(path), atom_style="atomic").structure
    return Structure.from_file(str(path))


def _written_elements(path: Path):
    """Element order as written in the file itself (the sidecar order contract)."""
    text = path.read_text().splitlines()
    name = path.name.lower()
    if name.endswith(".cif"):
        out = []
        for line in text:
            tok = line.split()
            if len(tok) == 7 and tok[2] == "1" and tok[1].startswith(tok[0]):
                out.append(tok[0])
        return out
    if name.endswith(".data"):
        types = {}
        for line in text:
            m = re.match(r"^\s*(\d+)\s+[\d.]+\s+#\s*(\w+)", line)
            if m:
                types[m.group(1)] = m.group(2)
        start = next(i for i, line in enumerate(text) if line.strip() == "Atoms")
        rows = []
        for line in text[start + 1:]:
            tok = line.split()
            if not tok:
                if rows:
                    break
                continue
            rows.append((int(tok[0]), types[tok[1]]))
        assert [r[0] for r in rows] == list(range(1, len(rows) + 1))
        return [r[1] for r in rows]
    symbols, counts = text[5].split(), [int(x) for x in text[6].split()]
    return [el for el, n in zip(symbols, counts) for _ in range(n)]


def _align_to_sidecar(st, rows):
    """Fractional coordinates of the parsed structure, in sidecar order (match within 1e-6)."""
    F = np.asarray(st.frac_coords, float)
    sym = [s.specie.symbol for s in st]
    out = np.zeros((len(rows), 3))
    used = set()
    for k, r in enumerate(rows):
        d = np.abs((F - np.asarray(r["frac"]) + 0.5) % 1.0 - 0.5).max(axis=1)
        j = int(np.argmin(d))
        assert d[j] < 1e-6 and sym[j] == r["element"] and j not in used
        used.add(j)
        out[k] = F[j]
    return out


def independent_verify(path: Path, tau: float = 0.0):
    """Brute-force V1-V6 from D1/D2 on a re-parsed file plus its sidecar (no lta_cage_salt).

    Limit: the per-ion floors are taken from the record (they rest on the V_eff
    estimator, which oracle 1 checks separately), not recomputed here.
    """
    side = json.loads(Path(str(path) + ".salt.json").read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == side["structure_sha256"]
    st = _read_structure(path)
    rows = side["atoms"]
    assert _written_elements(path) == [r["element"] for r in rows]
    M = np.asarray(st.lattice.matrix, float)
    F = _align_to_sidecar(st, rows)
    csvr = _csv_rows()
    charges = {"Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Sr": 2, "Ba": 2}
    salt = [i for i, r in enumerate(rows) if r["role"] == "cage_salt"]
    ring = [i for i, r in enumerate(rows) if r["role"] in ("framework", "occlusion_counterion")]
    fw = [i for i, r in enumerate(rows) if r["role"] == "framework_atom"]
    salt_cat = [rows[i]["element"] for i in salt if rows[i]["element"] != "Cl"]
    # D1 §6: mole fractions over the cage-salt cations of the cell, else over all cations
    pool = salt_cat or [rows[i]["element"] for i in ring]
    x = {m: pool.count(m) / len(pool) for m in set(pool)}

    def dist(i, j):
        d = F[j] - F[i]
        d -= np.round(d)
        return float(np.linalg.norm(d @ M))

    def r_cut(a, b):
        if a == b == "Cl":
            return sum(w * csvr[("Cl", "Cl", f"around {m}")][0] for m, w in x.items())
        if {a, b} == {"Cl", "O"}:
            return sum(w * csvr[("Cl", "O", f"around {m}")][0] for m, w in x.items())
        if "O" in (a, b):
            return csvr[(a if b == "O" else b, "O", "")][0]
        if "Cl" in (a, b):
            return csvr[(a if b == "Cl" else b, "Cl", "")][0]
        return csvr[(a, b, "Cl-bridged")][0]

    centres = {"alpha-0": (.25, .25, .25), "alpha-1": (.75, .75, .75), "beta-0": (0, 0, 0), "beta-1": (.5, .5, .5)}
    checks = {}
    # V1 membership and inventory
    ok = True
    for i in salt:
        d = {c: float(np.linalg.norm((((F[i] - np.array(f)) + 0.5) % 1.0 - 0.5) @ M)) for c, f in centres.items()}
        home = rows[i]["cage_id"]
        ok &= min(v for c, v in d.items() if c != home) - d[home] >= tau
    inv = Counter((rows[i]["cage_id"], rows[i]["element"]) for i in salt)
    want = Counter({(c, e): n for c, d in side["request"]["cages"].items() for e, n in d.items()})
    checks["V1"] = ok and inv == want and len({rows[i]["site_id"] for i in ring}) == len(ring)
    # V2 hard minima
    ok = True
    for n, i in enumerate(salt):
        ei = rows[i]["element"]
        for j in salt[n + 1:] + ring:
            ok &= dist(i, j) >= r_cut(ei, rows[j]["element"]) + tau
        for j in fw:
            d = dist(i, j)
            ok &= d >= 1.5 + tau
            if rows[j]["element"] == "O":
                ok &= d >= r_cut(ei, "O") + tau
            elif rows[j]["element"] in ("Si", "Al"):
                ok &= d >= 3.0 + tau  # D1 r8 S-13 T_min
    checks["V2"] = ok
    # V3 counts (rule B), anchoring and caps
    floors = side["floors"]
    u_m = {"Li": 6, "Na": 6, "K": 7, "Rb": 7, "Cs": 6, "Sr": 8, "Ba": 9}
    ok = True
    cations = [i for i in salt + ring if rows[i]["element"] != "Cl"]
    cl = [i for i in salt if rows[i]["element"] == "Cl"]
    for i in salt:
        e, c = rows[i]["element"], rows[i]["cage_id"]
        others = cations if e == "Cl" else cl
        n_nb = 0
        for j in others:
            rmin = csvr[(rows[j]["element"] if e == "Cl" else e, "Cl", "")][1]
            d = dist(i, j)
            ok &= abs(d - rmin) >= tau
            n_nb += d < rmin
        if e == "Cl" or c.startswith("beta"):
            ok &= n_nb >= 1
        else:
            ok &= floors[c]["floor_by_species"][e] <= n_nb <= u_m[e]
    checks["V3"] = ok
    # V4 window charges from geometry (nearest cages), V5 invariants
    q = {c: 0.0 for c in centres}
    zmax = {c: 0 for c in centres}
    for c, d in side["request"]["cages"].items():
        q[c] += sum(charges[e] * n for e, n in d.items() if e != "Cl") - d.get("Cl", 0)
        zmax[c] = max([charges[e] for e in d if e != "Cl"] or [0])
    fw_raw = {c: 0.0 for c in centres}
    for i in ring:
        e, fam = rows[i]["element"], rows[i]["site_id"][0]
        d = {c: float(np.linalg.norm((((F[i] - np.array(f)) + 0.5) % 1.0 - 0.5) @ M)) for c, f in centres.items()}
        alpha = min(("alpha-0", "alpha-1"), key=d.get)
        beta = min(("beta-0", "beta-1"), key=d.get)
        w = {"4": {alpha: 1.0}, "6": {beta: 0.5, alpha: 0.5}, "8": {"alpha-0": 0.5, "alpha-1": 0.5}}[fam]
        for c, wc in w.items():
            zmax[c] = max(zmax[c], charges[e])
            if rows[i]["role"] == "occlusion_counterion":
                q[c] += charges[e] * wc
            else:
                fw_raw[c] += charges[e] * wc
    for grp in (("alpha-0", "alpha-1"), ("beta-0", "beta-1")):
        mean = sum(fw_raw[c] for c in grp) / 2
        for c in grp:
            q[c] += fw_raw[c] - mean
    checks["V4"] = all(abs(q[c]) <= (zmax[c] or 1) + 1e-12 for c in centres)
    q_fw = sum(charges[rows[i]["element"]] for i in ring if rows[i]["role"] == "framework")
    q_exc = sum(charges[rows[i]["element"]] for i in ring + salt
                if rows[i]["role"] != "framework" and rows[i]["element"] != "Cl")
    checks["V5"] = q_fw == 24 and q_exc == len(cl)
    # V6 single image: half the smallest perpendicular width exceeds max r_min1 + 0.1 A
    a, b, c = M
    vol = abs(float(np.dot(a, np.cross(b, c))))
    widths = [vol / np.linalg.norm(np.cross(b, c)), vol / np.linalg.norm(np.cross(c, a)),
              vol / np.linalg.norm(np.cross(a, b))]
    r_min1_max = max(v[1] for k, v in csvr.items() if k[1] == "Cl" and k[2] == "" and k[0] in charges)
    checks["V6"] = min(widths) / 2 > r_min1_max + 0.1
    # ring cations sit on a template variant of their site (the template, not the cage-salt tables)
    _, _, sites, _ = P.load_lta_template(REAL)
    variants = {s_.site_id: [np.asarray(v.frac, float) for v in s_.variants] for r in sites.values() for s_ in r}
    checks["V1"] = checks["V1"] and all(
        min(np.abs((v - F[i] + 0.5) % 1.0 - 0.5).max() for v in variants[rows[i]["site_id"]]) < 1e-6 for i in ring)
    return checks


@pytest.mark.parametrize("name", ["out.cif", "out.POSCAR", "out.data"])
def test_independent_verifier_on_reparsed_output(counterion_result, tmp_path, name):
    path = tmp_path / name
    counterion_result.write(path)
    side = json.loads(Path(str(path) + ".salt.json").read_text())
    assert side["schema_version"] == 1 and side["status"] == "passed" and side["structure_written"] is True
    assert [r["index"] for r in side["atoms"]] == list(range(len(side["atoms"])))
    assert independent_verify(path) == {"V1": True, "V2": True, "V3": True, "V4": True, "V5": True, "V6": True}


def test_sidecar_sha256_detects_a_stale_structure(salt_result, tmp_path):
    path = tmp_path / "s.cif"
    salt_result.write(path)
    path.write_text(path.read_text() + "\n")
    with pytest.raises(AssertionError):
        independent_verify(path)


# ----------------------------------------------------------------------------- CLI


def _cli(script, out, *args, cwd=ROOT):
    return subprocess.run([sys.executable, str(ROOT / script), str(REAL), str(out), *args],
                          cwd=cwd, capture_output=True, text=True)


def test_cli_end_to_end_and_failure_paths(tmp_path):
    base = ["--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6"]
    ok = _cli("lta_cation_population_grouped.py", tmp_path / "ok.data", *base,
              "--alpha-salt", "Cl:8,Li:5,K:3", "--beta-salt", "Cl:1,Li:1", "--seed", "4")
    assert ok.returncode == 0, ok.stdout + ok.stderr
    assert independent_verify(tmp_path / "ok.data") == {"V1": True, "V2": True, "V3": True, "V4": True, "V5": True, "V6": True}
    # up-front validation error: exit 2, nothing written
    bad = _cli("lta_cation_population_grouped.py", tmp_path / "bad.cif", *base, "--alpha-salt", "Cl:13,Li:13")
    assert bad.returncode == 2 and "cap of 12" in bad.stderr
    assert not (tmp_path / "bad.cif").exists() and not (tmp_path / "bad.cif.salt.json").exists()
    # legacy flags without the selector
    leg = _cli("lta_cation_population_grouped.py", tmp_path / "leg.cif", *base, "--alpha-cl", "1")
    assert leg.returncode == 2 and "sparse-v715" in leg.stderr
    # construction failure: failure sidecar only, existing structure untouched
    keep = tmp_path / "fail.cif"
    keep.write_text("old\n")
    fail = _cli("lta_cation_population_grouped.py", keep, *base, "--alpha-salt", "Cl:4,Li:2,K:2",
                "--volume-factor", "0.3", "--seed", "1")
    assert fail.returncode == 1
    assert keep.read_text() == "old\n"
    side = json.loads((tmp_path / "fail.cif.salt.json").read_text())
    assert side["status"] == "failed" and side["structure_written"] is False and "warning" in side
    # rejected modes
    for extra, needle in ((["--unique-batch-size", "2"], "unique-batch-size"),
                          (["--auto-rebalance"], "auto-rebalance"),
                          (["--write-candidate", str(tmp_path / "c.cif")], "write-candidate"),
                          (["--population-mode", "maximum-symmetry"], "maximum-symmetry")):
        r = _cli("lta_cation_population_grouped.py", tmp_path / "r.cif", *base, "--alpha-salt", "Cl:2,Li:2", *extra)
        assert r.returncode == 2 and needle in r.stderr


def test_smoothing_attempt_seed_reproduces_bytes(tmp_path):
    base = ["--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-salt", "Cl:6,Li:3,K:3",
            "--beta-salt", "Cl:1,Li:1", "--max-attempts", "3", "--mode", "regenerate"]
    first = _cli("generate_smooth_lta_cations_grouped.py", tmp_path / "a.cif", *base, "--seed", "9")
    assert first.returncode == 0, first.stdout + first.stderr
    seed = re.search(r"^Random seed: (\d+)$", first.stdout, re.M).group(1)
    second = _cli("generate_smooth_lta_cations_grouped.py", tmp_path / "b.cif", *base, "--attempt-seed", seed)
    assert second.returncode == 0, second.stdout + second.stderr
    assert (tmp_path / "a.cif").read_bytes() == (tmp_path / "b.cif").read_bytes()
    assert independent_verify(tmp_path / "b.cif")["V2"]
