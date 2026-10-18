---
toc: true
toc-depth: 3
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - |
    \usepackage{booktabs}
    \usepackage{longtable}
    \usepackage{microtype}
  - |
    \setlength{\emergencystretch}{3em}
---

# v7.16 cage-salt construction: numerical method (D2)

**Lifecycle: ACCEPTED-CURRENT D2 (r2; human acknowledgement 2026-10-18).**
The stakeholder acknowledged the §12 items and the D1 readings in chat:
"Sounds ok. These are very specific items. I'll need to run the program later
myself in order to figure out the problems." An AI agent recorded this.

- **Validation route.** The stakeholder's own runs are an expected
  validation route. Problems they find reopen this D2 through its owner.
- **What it does not do.** It does not self-qualify an implementation.

History of the draft:
- **Drafted:** 2026-10-17 by an AI agent (Claude Code session
  https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi), under workplan
  `LTA-V7.16-MOLTEN-SALT-OCCLUSION`, SSDP 7.2.1.
- **Parent:** the accepted-current D1
  `doc/specs/v7_16_molten_salt_occlusion_extension.md` (r7, human-ratified
  2026-10-17), cited as **D1 §n**.
- **r1 folds in independent D2 review 1 (§13).** That review raised two
  SERIOUS CHALLENGES and 15 findings.
- **r2 folds in D2 review 2,** which raised no SERIOUS CHALLENGE and eight
  new findings, N1-N8 (§13).
- **Status of choices:** every parameter value is an agent-proposed D2
  default (§10) until this document is accepted.
- **Acceptance gate:** independent review, then human acknowledgement of the
  output-changing choices (§12).

## 1. Scope and what D2 must preserve

These D1 invariants hold throughout:
- **Inventory.** Inventory and charge invariants are exact (D1 §4), and each
  `cage_salt` ion stays in its own cage (S-10).
- **Hard minima.** d_ij ≥ r_cut(i, j) holds on the written coordinates (D1
  §6), and so does the 1.5 Å guard.
- **Coordination.** Anchoring holds, together with
  max(1, ⌈a·C_M(c)⌉) ≤ N_i(Cl) ≤ U_M (D1 §7, S-11, a = 0.4). Under the S-12
  fallback, rule A with a = 0.7 applies instead.
- **Cage charge.** |q_c| ≤ q_max(c) (D1 §4).
- **Ring relocation.** D-12 limits apply (D1 §8): a ring cation moves only to
  a site or variant that *resolves* its conflict; otherwise it stays.
- **Order and streams.** The construction order and stream isolation of D1
  §11 apply.
- **Fail-closed.** No structure is written unless every check passes.
- **No sampling claim.** No sampling law is claimed (D-9).

**D2 readings of D1. The D1 owner must confirm both before D2 is
accepted.**
- **(R-1) Cage membership applies to `cage_salt` ions only.** D1 §4 requires
  every "role, cage and species count" to equal its input. Ring roles have
  per-family inputs, not per-cage inputs, so for them the exact check is the
  count per (role, family, species), plus "each occupied site holds one
  cation at a template variant of that site".
  - **This is a narrowing of D1 §3.** It waives D1 §3's "Ties are a
    construction failure" for ring cations.
  - Why it is needed: read literally, D1 would fail every structure with an
    occupied 8R site. Every 8R site anchor and variant is equidistant from α₀
    and α₁ (anchors 6.139/6.139 Å; variants 6.2668/6.2668 Å).
- **(R-2) Window weights w(j, c) are fixed per physical site** by ring-window
  topology on the canonical, gated template (§4.3). They are never derived
  from distances to individual variants.

## 2. Definitions

- **Cell and metric.**
  - L is the lattice matrix of the cell in use (rows are lattice vectors, Å),
    and r = f L.
  - d(r₁, r₂) = ‖(f − round(f)) L‖ with f = (r₁ − r₂) L⁻¹ is the
    minimum-image distance. It is single-image exact when
    max r_min1 + 0.1 Å < ½ min_k(perpendicular width k).
  - That precondition holds on the canonical cell (7.09 Å against 5.11 Å;
    the 0.1 Å margin is an agent choice). It is checked at run time, and a
    violation fails closed.
  - The cut-off is a volume factor of about 0.37 for Na-A, and it depends on
    the population through V_ref.
  - Nearest-centre and membership results remain exact even where centre
    distances exceed half the width. The farthest any point lies from its
    nearest centre is 6.84 Å, below 7.09 Å, and that ratio is unchanged by
    scaling.
- **Ions and sets.**
  - S: `cage_salt` ions.
  - R: ring cations (`framework` and `occlusion_counterion`).
  - F_O: framework O atoms; F: all framework atoms.
- **Pair tables.**
  - r_cut and r_min1 come from the D1-bound E-pair CSV, with the D1 §6 rows
    and weights x_M.
  - Notation: w(j, c) denotes the D1 window weights. The smooth-count width
    is written w_s.
  - The build target is ρ_cut = r_cut + δ_b.
- **Cage centres and membership.**
  - c ∈ {α₀, α₁, β₀, β₁}, at the D1 §3 canonical fractional positions.
  - home(i) is the cage of i ∈ S.
  - margin(i) = min_{c' ≠ home(i)} d(i, c') − d(i, home(i)).
- **Coordination count.**
  - j is a first neighbour of i iff they have opposite charge and
    d(i, j) < r_min1(i, j).
  - N_i counts first neighbours over all ions, in any role and any cage.
- **Floors and caps.**
  - For α `cage_salt` cations: L_i = max(1, ⌈a·C_M(c) − 10⁻⁹⌉) and
    U_i = U_M.
  - For β `cage_salt` cations and all Cl: L_i = 1 and there is no cap.
- **Writer tolerance band.** τ_w = 10⁻⁵ Å (§7).

## 3. Pipeline (concretizes D1 §11; scaling order fixed after review finding 9)

| Step | Operation | Cell | Stream |
|---|---|---|---|
| 1 | framework ring placement | V_ref cell (volume factor 1) | legacy `random.Random(seed)` |
| 2 | `occlusion_counterion` placement | V_ref cell | v7.15 secondary (`random.Random`, SHA-256 `v7.15-occlusion:{seed}`) |
| 3 | cation density smoothing (smoothing entry point only) | V_ref cell | existing |
| 4 | charge rebalancing (§5) | topology only | R_q |
| 5 | isotropic scaling to V_target = f·V_ref (§4.1) | → final cell | none |
| 6 | V_eff, n_melt, floors (§4.2) | final cell | none |
| 7 | cage-salt build with D-12 relocation (§6) | final cell | R_s,k for start k |
| 8 | exact verification (§7) | final cell | none |
| 9 | write | final cell | none |

- **Smoothing on the V_ref cell.** In cage-salt mode, smoothing always runs
  on the V_ref cell, which is the v7.14 dry reference at volume factor 1. The
  salt-dependent f_V, or an explicit override, is applied only at step 5.
  - Consequence: the step-3 state is independent of every cage-salt input
    other than the population and counterion allocations. That state covers
    the placed ring ions, sites, variants and the V_ref cell.
  - With an explicit `--volume-factor`, steps 1-3 run at factor 1, not at the
    override. This is the D2 reading of D1 §11 "must not perturb".
  - This differs from v7.15, which scales before smoothing. Legacy and dry
    paths keep their existing behaviour.
- **Smoothing entry point.** Steps 4-9 run once, on the population the
  attempt loop selected, using that attempt's integer seed. That seed is
  drawn from `master_rng` and logged, so smoothing-entry-point runs are
  reproducible from the logged seed even when the user gives no seed.
- **Derived streams.**
  - seed_int(tag) = int.from_bytes(SHA-256(tag)[:8], "big").
  - R_q = `numpy.random.Generator(PCG64(seed_int("v7.16-charge-rebalance:{seed}")))`.
  - R_s,k = the same construction with tag `v7.16-cage-salt:{seed}:start{k}`.
    Each restart has its own substream, so start k does not depend on how
    many draws earlier starts consumed.
  - With seed = None, R_q and R_s,k are seeded from OS entropy, not from a
    tag containing "None", and the run is not reproducible, as in v7.15.
  - Reproducing a smoothing-entry-point run from its logged attempt seed needs
    a way to run that single attempt seed. That is a D3 obligation (§11).

## 4. Deterministic quantities

### 4.1 Volume

- n̄ = N_Cl/2 and f_V(n̄) = 1 − 0.03·min(n̄/10.5, 1).
- An explicit `--volume-factor` replaces f_V (D1 §9).
- The cell is scaled isotropically by s = f^{1/3}, and fractional
  coordinates are unchanged.

### 4.2 Effective volume (revised after review finding 5)

**Estimator.** Midpoint quadrature on a *fixed* fractional grid of size
G = 192, the same for every cell scale:

```
f_ijk = ((i + 1/2)/G, (j + 1/2)/G, (k + 1/2)/G)
V_eff(c; d) = (V_cell / G^3) * sum_ijk  w_ijk(c) * [ min_{o in F_O} d(f_ijk, o) > d ]
```

- **Tie sharing.** w_ijk(c) = 1/m if c is one of the m cage centres
  equidistant (within 10⁻⁹ Å) from f_ijk at the minimum distance, and 0
  otherwise. Ties are shared, never assigned by order.
- **Why G is commensurate.** G is a multiple of 8, matching the
  (¼, ¼, ¼) centre lattice. Symmetric cages then get identical volumes. A
  non-commensurate G aliases: G = 174 gave ±1 % at d = 0.
- **Scale consistency.** The grid is fixed in fractional coordinates, so
  V_eff(scaled cell; d) = s³·V_eff(V_ref cell; d/s) holds exactly. n_melt
  therefore has no grid jitter as the loading changes f_V.
- **Accuracy (measured; agent, 2026-10-17, fixture at V_ref).**
  - At d = 0 each cell must be exactly V/4 = 925.823 Å³; the estimator gives
    925.82 Å³ for all four cages.
  - At d = 1.62 Å it gives α₀ = α₁ = 809.05 Å³, against a 2×10⁷-point
    Monte Carlo reference of 808.6 ± 0.34 Å³ (review 1), i.e. +0.06 %.
  - The requirement is |error| ≤ 0.5 %. It is checked by these two oracles at
    qualification, not by grid refinement, because the error is not monotone
    in G.
- **Cost.** About 8 s for G = 192 with a KD-tree over periodic images. It is
  computed once per structure, or cached per V_ref cell and scaled through
  the identity.
- **Derived quantities.** d_wall(c), x_Li(c), v_pair, n_melt(c), the range,
  φ(c), C_M(c) and L_i follow exactly as D1 §7 defines them. They are fixed
  before step 7.

### 4.3 Window classes and charge (revised after review findings 2-3)

**Window table (fixed per physical site, from the canonical template's
`anchor_frac`; verified in review 2).** This is a property of the gated template, so it is stored
as a table keyed by site_id and checked at load.

| Family | Rule | Weights | Classes on the template |
|---|---|---|---|
| 4R | the nearest **α** centre (the D4R-face anchor is nearer to β, at 6.445 Å, than to its α, at 6.719 Å, so "nearest centre" would be wrong) | w = 1 | 4R@α₀ (12 sites), 4R@α₁ (12); α gap ≥ 2.18 Å |
| 6R | the nearest β and the nearest α | ½ / ½ | 6R@(β_k, α_l), 4 sites each; gaps ≥ 6.4 Å |
| 8R | the two α sharing the window | ½ / ½ | 8R@(α₀, α₁), 6 sites |

**Charge and threshold.**
- q_c^exc and q_c^fw are as in D1 §4, held as exact fractions. q_c^fw
  subtracts a same-type mean, so q_c is a multiple of ¼.
- q_max(c) = z_max(c)·e. z_max(c) is the largest formal charge among the
  cations that contribute to c: its `cage_salt` cations and every ring cation
  with w(j, c) > 0. It is evaluated *on the configuration under test*, so a
  counterion moving into a window can raise it.
- **Excess-only check.** In the up-front check (§5), framework cations are
  excluded, both from q and from z_max.
- **No contributing cation.** If no cation contributes to c, z_max(c) = 1.
  This is a D2 default filling an undefined case; it is the smallest
  supported charge, so it is conservative.

## 5. Charge rebalancing (step 4): exact

**Decision variables.**
- Counterions are aggregated into groups (family, charge z ∈ {1, 2}).
  Species of equal charge are interchangeable for q_c and z_max.
- The decision is the integer vector n[g, class]: how many counterions of
  group g occupy each class of g's family.
- The capacity of a class is its number of sites that are empty after step 3
  or hold an `occlusion_counterion`. It is shared across the charge groups of
  the family: Σ_z n[(f, z), class] ≤ capacity(class).

**Search size (corrected).** The number of candidate vectors is
Π_g C(n_g + k_g − 1, k_g − 1), with k = 2, 4 and 1 classes for 4R, 6R and 8R.
- Over all admissible inputs (Q_exc ≤ 26), the worst case is 163,072
  vectors, or 35,840 once per-class capacities are applied (both recomputed
  by brute force in review 2). Agent enumeration found it at
  (4R: 3 singly and 1 doubly charged; 6R: 11 singly and 5 doubly charged).
- Each vector costs O(1) to evaluate, so exhaustive enumeration is exact and
  costs well under a second.

**Objective.** Among vectors meeting every capacity and every
|q_c| ≤ q_max(c), choose lexicographically by:
1. the number of moved counterions,
   Σ_g Σ_class max(0, n₀[g, class] − n[g, class]), where n₀ is the step-3
   vector;
2. Σ_c q_c²;
3. the lexicographically smallest vector in canonical order. Groups are
   ordered by family (4, 6, 8), then charge; classes are ordered 4R@α₀,
   4R@α₁, 6R@(β₀,α₀), 6R@(β₀,α₁), 6R@(β₁,α₀), 6R@(β₁,α₁), 8R.

**Infeasible.** If no vector is feasible, this is a seed-level failure. The
run reports min over vectors of max_c(|q_c| − q_max(c)) and the framework
sites blocking it.

**Up-front check.** The same search, with every site of each family
available and framework cations excluded, is the deliberately stricter input
rule of D1 §4.

**Realization (deterministic up to R_q).** Within each group, for each class
with a surplus:
- the surplus counterions leave in ascending (site_id, species) order;
- each is assigned to a class with a deficit, in canonical class order;
- its target site is drawn uniformly from that class's empty sites with R_q,
  and its variant uniformly from that site's variants with R_q, as v7.15
  placement does.

The report keeps only the chosen vector, its objective, and the next 5 best
vectors.

## 6. Cage-salt build (step 7)

### 6.1 Variables and objective

- **Variables.** The positions x ∈ ℝ^{3|S|}. Ring and framework atoms are
  fixed except through §6.4.
- **Penalty function.** h(u) = max(0, u)².
- **Smooth count.** s_i = Σ_{j opposite} σ((r_min1(i, j) − d_ij)/w_s).

```
E_W(x) = W * [ sum_{i<j in S} h(rho_cut(i,j) - d_ij)
             + sum_{i in S, k in R} h(rho_cut(i,k) - d_ik)
             + sum_{i in S, o in F_O} h(rho_cut(i,o) - d_io)
             + sum_{i in S, k in F} h(1.5 + delta_b - d_ik)
             + sum_{i in S} h(eps_m - margin(i)) ]
         + W_dz * sum_{i in S, j opposite to i} h(delta_cn - |d_ij - r_min1(i,j)|)
         + sum_{i in S, cation} [ h(L_i + mu - s_i) + h(s_i - (U_i - mu)) ]
         + sum_{i in S, Cl}      h(1 + mu - s_i)
```

W_dz = 0 in the three continuation stages. It is W_dz = 10 only in a final
polish stage (§6.3).

- **Dead-zone term (added after review 1 finding 6; reweighted after review
  2 N4).** It pushes every opposite-charge pair at least δ_cn away from its
  counting boundary r_min1, so smooth and discrete counts agree.
  - It carries its own weight and acts only in the final polish. Weighted by W
    from the start, it would have frozen the counts after the first stage.
  - It depletes a ±0.05 Å band around r_min1 in every opposite-charge pair
    distribution. That changes the output; it is reported (§8.1) and listed
    for acknowledgement (§12).
  - It is not covered by probes P6-P8.
- **The count terms remain a heuristic.** For example, a smooth count of
  2.59 is reachable with one true neighbour and four just outside r_min1.
  Only the exact verification (§7) decides acceptance.
- **Smoothness.** E_W is piecewise C¹: h is C¹. It is not differentiable:
  - where margin(i) switches the minimizing centre;
  - where the minimum image switches;
  - at d = 0;
  - at d = r_min1 when W_dz > 0, where the dead-zone gradient jumps. L-BFGS-B tolerates these measure-zero sets in
  practice, and §7 catches any failure.

### 6.2 Initialization (start k, stream R_s,k)

- **Order.** Ions are inserted in an order drawn from R_s,k.
- **Candidates.** Each position is drawn uniformly in the cell. A candidate
  is accepted if:
  - margin ≥ ε_m;
  - it is at least ρ_cut from every framework O, ring cation and placed salt
    ion;
  - it is at least 1.5 + δ_b from every framework atom.
- **Fallback.** After M_ins = 2000 rejections, the candidate with the
  smallest sum of squared violations is kept.

### 6.3 Relaxation and shaping

- **Minimizer.** L-BFGS-B (SciPy) minimises E_W at W = 10, 100 and 1000 in
  turn, with W_dz = 0. A final polish stage follows at W = 1000 and W_dz = 10.
  Each stage starts from the previous result.
- **Options:** maxiter 6000, maxfun 15000, maxcor 10, ftol 2.2×10⁻⁹,
  gtol 10⁻⁵, maxls 20.
- The minimiser's success flag is ignored; §7 decides.

### 6.4 D-12 relocation (corrected after review finding 4)

If §7 fails on a hard minimum involving a ring cation k:
- **Candidates.** k is offered its D1 §8 candidates: empty sites of its
  family, and the other variants of its own site.
- **Acceptance.** A candidate is accepted only if it *resolves* every
  violation between k and S at the current x (all d ≥ ρ_cut, guard
  satisfied) and keeps every |q_c| ≤ q_max(c).
- **Choice.** Among accepted candidates, the one with the largest minimum
  slack against S is taken; ties go to ascending (site_id, variant_id). If
  none qualifies, k stays.
- **Order.** Conflicting ring cations are processed in descending order of
  their maximum violation; ties go by site_id.
- **Rounds.** After any move, §6.3 is re-run from x. At most R_reloc = 4
  rounds are allowed per start.

### 6.5 Restarts and fail-closed

- Starts k = 1..K_max with K_max = 20 are run in order, and the first that
  passes §7 is returned (first-success). The output of a successful run is
  therefore independent of K_max.
- If no start passes, the run fails closed. The report gives, per start, the
  best slack per constraint class and the CN deficit and excess.
- **Cost (review estimate).** Typically 10-20 s per start, and a worst case
  of about 3.7 min with all relocation rounds. A failing request at
  K_max = 20 can take about 75 min. That is acceptable for an
  initial-structure generator, and it is reported.

### 6.6 Fallback rule A (S-12, only if adopted at qualification)

The cation count terms become a smooth version of the per-cage, per-species
mean condition:

```
sum_{i in M, c} N_i >= sum_{i in M, c} 0.7 * C_M(c)
```

- **Pooling for small populations (D1 §7).**
  - A species is evaluated per cage wherever it has at least 3 centres there.
  - The α-cages where it has fewer than 3 are pooled together, each ion
    carrying its own cage's C. A cage evaluated alone is never pooled, so no
    ion is counted twice.
  - If the pooled set still has fewer than 3 centres, only anchoring and U_M
    apply. In the census this branch covers n = 2 and every Cs/Sr/Ba dopant.
- Anchoring and U_M stay per ion, and everything else is unchanged.

## 7. Exact verification (step 8; revised after review SC-2)

**What is verified.**
- The canonical coordinates are the fractional coordinates rounded to 10
  decimals and wrapped into [0, 1), together with the lattice.
- All writers serialize from these. The existing writers quantize further:
  CIF and POSCAR to 8 decimals, and LAMMPS to 8 decimals in its rotated
  Cartesian frame. Including the lattice parameters and the rotation, the
  worst case is about 2.6×10⁻⁷ Å per atom and 5.3×10⁻⁷ Å per pair, so
  τ_w = 10⁻⁵ Å covers it (measured in review 2).
- Verification is brute-force O(N²) in float64 on the canonical coordinates,
  with the margins below. Every check therefore keeps its verdict after any
  writer's quantization.

| # | Check | Margin |
|---|---|---|
| V1 | `cage_salt` inventory per cage and species, by nearest centre; ring cations as counts per (role, family, species) equal to the input, with each occupied site holding one cation at a template variant (R-1) | margin(i) ≥ τ_w |
| V2 | d_ij ≥ r_cut(i, j) + τ_w for every pair with at least one `cage_salt` ion, against ring cations and framework O; d ≥ 1.5 + τ_w to every atom | +τ_w (stricter than D1) |
| V3 | Anchoring and L_i ≤ N_i ≤ U_i (rule A if adopted), with every opposite-charge pair |d − r_min1| ≥ τ_w | exact counts |
| V4 | Every |q_c| ≤ q_max(c) | exact fractions |
| V5 | Q_fw = 24 and Q_exc = N_Cl | exact |
| V6 | Single-image precondition (§2) | — |

D4 adds a writer round-trip test: re-parse each written file and re-run V1-V5
with τ_w = 0. The verdicts must agree. Two conditions apply (review 2):
- **Non-snapping reader.** The round-trip uses a reader that does not snap
  coordinates, for example pymatgen `CifParser(frac_tolerance=0)`. The
  default reader snaps near-common fractions, moving framework O by about
  1.2×10⁻⁴ Å. Users re-reading the CIF should be warned of this.
- **Role and site sidecar.** Roles and site ids are not representable in CIF,
  POSCAR or LAMMPS, and the writer reorders atoms by element. A per-atom
  sidecar in written order is therefore required (D3, §11).

## 8. Stochastic, error and reproducibility semantics

- **Random elements.**
  - Insertion order and positions, R_s,k.
  - The choice of site and variant within a class during rebalancing, R_q.
  - Nothing else: relocation choices are deterministic.
- **No sampling law.** No sampling law is claimed. Outputs are biased toward
  pair distances near ρ_cut and counts just inside the floor and cap
  margins, because penalty relaxation settles on constraint boundaries. This
  is accepted under D-9 and made visible (§8.1).
- **No probability bound.** D2 gives no bound on the probability of a CN miss
  or a failure outside the S-12 grid. That grid is a census, not a sample.
  - This is a recorded deviation from the D1 §16 hand-off, which lists "how
    the CN miss probability is bounded" as D2-owned.
  - D2 bounds it only by the census, and asks the D1 owner to acknowledge
    this (§12).
- **Reproducibility envelope.** Output is byte-identical for the same
  inputs, seed and recorded stack:
  - Python version;
  - NumPy (Generator stream stability is not guaranteed across versions,
    NEP 19), SciPy (L-BFGS-B implementation) and pymatgen (writers);
  - BLAS library and thread count;
  - platform.

  Across stacks, outputs may differ but still satisfy §7. The S-12 census is
  re-run on any stack change in a release. The legacy model and the dry path
  are unaffected.
- **Numerical error.**

  | Source | Bound | Effect |
  |---|---|---|
  | V_eff | ≤ 0.5 % (measured +0.06 %) | floors L_i can flip only if a·C lies within that band of an integer; reported with its distance |
  | float64 | ≪ τ_w | none |
  | non-convergence | — | restart or failure only, never a wrong output |

### 8.1 Numerical inspectability (O1)

Fields added to the D1 §14 per-run report:
- **Rebalancing:** the chosen class vector, its objective, the next 5 best,
  and the moves.
- **Volume:** G and V_eff at d_wall and ± 0.3 Å, and each floor's distance
  from an integer threshold.
- **Per start:**
  - E_W after each stage;
  - the best slack per constraint class (salt–salt, salt–ring, salt–O, guard,
    membership, count dead zone);
  - CN deficit and excess;
  - relocation moves;
  - the failing check.
- **Final structure:**
  - the V1-V6 table;
  - histograms of d − r_cut per pair class (0-0.5 Å, 0.05 Å bins);
  - a histogram of d − r_min1 for opposite-charge pairs (±0.5 Å, 0.01 Å
    bins), which shows the dead-zone band;
  - N_i − L_i per species;
  - the start used, and K_max;
  - the stack versions.

**Reader and questions.**
- Did it pass, and on which start?
- How tight was the packing?
- Was the result boundary-hugging?
- Which ring cations moved?
- Were any floors near a threshold?

**Product surface.** These fields extend the D1 §14 product surface. Technical
acceptance does not accept them; the human acknowledgement (§12) does.

## 9. Evidence and oracles (produced in D4; defined here first)

1. **V_eff:**
   - d = 0 gives exactly V/4 per cage;
   - at d = d_wall, the Monte Carlo reference (≥ 10⁷ points) agrees within
     0.5 %;
   - α₀ = α₁ and β₀ = β₁ by symmetry;
   - the scale identity V(s·L; d) = s³·V(L; d/s) holds.
2. **An independent brute-force verifier**, written separately from the
   production checks, re-checks V1-V6 on re-parsed CIF, POSCAR and LAMMPS
   output.
3. **Manufactured feasibility.** A configuration known to satisfy every
   check passes. Single-check violations, one per check, are rejected by
   exactly that check.
4. **Rebalancing.** The exhaustive search is compared with brute force over
   individual site assignments on small synthetic problems. Σ_c q_c = 0
   always holds.
5. **Streams.** The step-3 state snapshot (ring ions, sites, variants, V_ref
   cell) depends only on the population and counterion allocations, the seed
   and the smoothing settings. It does not depend on the α/β salt inventories
   or on f_V. Changing R_s,k leaves steps 1-6 unchanged.
6. **Determinism.** Same seed and stack give identical bytes.
7. **Precondition.** A volume factor below the single-image limit fails
   closed.

**S-12 qualification census (pre-specified; all parameters frozen before
it is run; any re-tuning re-runs the whole census).**
- **Entry point.** The population entry point
  `generate_lta_cation_structure_from_specs`, with integer seeds 0-4 and the
  modelled f_V (no override).
- **Rings.** Na:2/16/6 and K:2/16/6, with no `occlusion_counterion`s.
- **α loadings,** per α-cage, both α-cages the same, Li/K split by
  round(0.592·n):

  | n | 2 | 4 | 6 | 8 | 10 | 12 |
  |---|---|---|---|---|---|---|
  | Li/K | 1/1 | 2/2 | 4/2 | 5/3 | 6/4 | 7/5 |

- **Dopant rows,** per α-cage, at seeds 0-4:
  - Li5K4Cs1 with 10 Cl;
  - Li5K3Sr1 with 10 Cl;
  - Li6K4Ba1 with 12 Cl.
- **β content.** Each β holds 1 Cl and 1 Li in every row.
- **Size.** 2 rings × 5 seeds × 9 compositions = 90 requests.
- **Decision rule.**
  - Rule B is retained iff at most 4 of the 90 fail.
  - Otherwise rule A (a = 0.7) becomes the release default and must pass the
    same criterion. If A also fails, that reopens D1.
  - Both pass rates are recorded, with Wilson intervals for information.

**Prior exploratory evidence.** P5-P8 (D1 §12) are probes of an earlier
variant of this method. They motivate it but are not acceptance evidence.

## 10. Parameter register (agent-proposed D2 defaults)

| Parameter | Value | Basis |
|---|---|---|
| δ_b build margin | 0.01 Å | ≫ the penalty residual (about 0.002 Å, review estimate) and τ_w |
| τ_w writer band | 10⁻⁵ Å | above the worst writer quantization (< 10⁻⁶ Å) |
| ε_m membership margin | 0.05 Å | agent choice, no external basis |
| δ_cn count dead zone | 0.05 Å | keeps smooth and discrete counts consistent |
| w_s smooth-count width | 0.12 Å | P6-P8 |
| μ count margin | 0.3 | agent choice, no external basis |
| W continuation | 10, 100, 1000 | P6-P8 |
| W_dz dead-zone weight | 0 in continuation; 10 in the final polish (W = 1000) | review 2 N4; agent choice, not probed |
| Single-image precondition margin | 0.1 Å | agent choice |
| L-BFGS-B options | maxiter 6000, maxfun 15000, maxcor 10, ftol 2.2e-9, gtol 1e-5, maxls 20 | SciPy defaults, except maxiter |
| M_ins and insertion fallback metric | 2000; minimum sum of squared violations | agent choice, no external basis |
| Relocation choice | resolving candidates only; maximum minimum slack; (site_id, variant_id) | D1 §8 |
| K_max, R_reloc | 20, 4 | agent choice; cost in §6.5 |
| V_eff grid | G = 192, fixed fractional, tie sharing within 10⁻⁹ Å | §4.2 |
| Ceiling guard | 10⁻⁹ | §2 |
| Canonicalisation | round to 10 decimals, then wrap into [0, 1) | §7 |
| Streams | SHA-256 tags, first 8 bytes big-endian, PCG64; per-start substreams | §3 |
| z_max with no contributor | 1 | §4.3 |
| Report sizes | 5 best alternatives; histogram bins 0.05 Å over 0-0.5 Å | §8.1 |
| S-12 census | §9 | stakeholder (rule) and agent (grid) |

Changing a value is a D2 change unless D3 exposes it as a recorded user
option. Only the seed and K_max are proposed as options.

## 11. D2 → D3 hand-off

D3 must preserve:
- the §3 order, including smoothing on the V_ref cell and scaling at step 5,
  and the stream derivation with per-start substreams;
- §4.2: a fixed G = 192 grid with tie sharing, and the scale identity;
- §4.3: the window table keyed by site_id and validated at load, exact
  fractional charge, and z_max(c) on the configuration under test;
- §5: aggregation by (family, charge), exhaustive search, the lexicographic
  objective and canonical orders;
- §6: the objective, including the dead zone; continuation; L-BFGS-B
  options; resolving-only relocation; first-success restarts;
- §7: canonical coordinates as the single source for every writer, the τ_w
  margins, brute-force verification, and the single-image precondition;
- §8.1: the report fields;
- SciPy and NumPy as runtime dependencies, and stack versions recorded in
  qualification;
- a per-atom role and site-id sidecar in written order, for the round-trip
  verification;
- a way to run a single logged attempt seed;
- a non-snapping reader in the round-trip tests.

D3 owns modules, data structures, neighbour lists outside verification,
vectorisation, CLI and schema, errors and report format.

## 12. Acceptance request

- **Requested:** human acknowledgement of the output-changing choices (r1
  has had an independent re-review, review 2):
  - first-success restarts with K_max = 20, and their runtime;
  - the boundary-hugging bias;
  - V1, V2 and V3 margins stricter than D1 by τ_w;
  - the dead-zone depletion band of ±0.05 Å around r_min1;
  - smoothing on the V_ref cell in cage-salt mode, which differs from
    v7.15's scale-then-smooth;
  - the S-12 census and decision rule;
  - the added report fields.
- **D1 owner confirmations required.**
  - R-1 narrows D1 §3's tie rule to `cage_salt` ions.
  - R-2: window weights are fixed per physical site.
  - The D1 §16 deviation: the CN miss probability is bounded only by the
    S-12 census.
- **Tensions.** Searched D1 §13 (T1-T22), the workplans, the D1 git log and
  the repository issues (none). Relevant: T8 (anchoring), T17/T21 (V_eff,
  handled deterministically) and T22 (fallback, §6.6). None is closed here.

## 13. Review record

**Independent D2 review 1** (2026-10-17). AI reviewer: a fresh subagent
context with the same model family. It is process-independent and carries
common-mode risk.

**SERIOUS CHALLENGES.**
- **SC-1: ring cations on 8R are equidistant from α₀ and α₁, so V1 failed by
  construction.** Fixed by R-1 and R-2.
- **SC-2: writers cannot emit the canonical coordinates exactly.** Fixed by
  the τ_w margins and a round-trip test.

**Findings and dispositions.**
- **Major:**
  1. The search bound was false → aggregation by (family, charge), worst case
     163,072 vectors.
  2. The feasible set was underdefined → z_max(c) on the configuration under
     test; excess-only up-front check; default 1; quarter-integers.
  3. Realization was ambiguous → canonical departure order, R_q for site and
     variant, canonical class order.
  4. Relocation conflicted with D1 → resolving-only, with a defined choice.
  5. The V_eff claims were false → fixed commensurate G = 192, tie sharing,
     the scale identity, and oracles in place of refinement. Measured: exact
     at d = 0 and +0.06 % at d_wall.
  6. The smooth count did not guarantee the discrete count → dead zone added,
     labelled a heuristic.
  8. Stream semantics → seeding specified, per-start substreams, seed = None
     corrected, envelope expanded.
  9. Scaling order → smoothing on the V_ref cell, scaling at step 5,
     oracle 5 reworded.
  10. S-12 was underspecified → fully specified, with rule A pooling defined.
  11. The parameter register was incomplete → completed.
- **Minor:**
  7. Differentiability wording → fixed.
  12. Single-image precondition → added as V6.
  13. Report size and product surface → top 5, with human acknowledgement.
  14. No probability bound outside the grid → stated.
  15. Insertion guard → added.

**Verified by the reviewer:**
- the window classes on the fixture;
- the up-front check wording;
- floors and caps cannot conflict;
- the penalty residual is below δ_b;
- the 10-decimal rounding bound;
- the S-12 arithmetic;
- first-success is independent of K_max;
- S-10 is preserved for `cage_salt` ions;
- cross-window counting;
- relocation cannot change the floors.

**Independent D2 review 2** (2026-10-17/18, of r1). AI reviewer: a fresh
subagent context, with the same common-mode caveat. No SERIOUS CHALLENGE.

**Closure of review 1.** The re-review closed findings 1-5, 12, 13 and 15.
It judged the rest partly closed, and finding 7 stale; the gaps are
addressed below.

**New findings and dispositions.**
- **N1:** the 4R rule must be "nearest α" → fixed, with anchor distances.
- **N2:** the default CIF reader snaps coordinates → non-snapping reader
  required; the "8 decimals" wording corrected.
- **N3:** roles are not in the written files → sidecar added to the D3
  hand-off.
- **N4:** the W-weighted dead zone froze the counts → separate W_dz in a
  final polish stage; reported and acknowledged.
- **N5:** the oracle 5 comparator did not exist → step-3 snapshot oracle;
  explicit-override reading stated.
- **N6:** rule A pooling incomplete → defined.
- **N7:** V1 wording → fixed.
- **N8:** workplan records, the tracked `.pyc`, the w/w_s notation → fixed.

**Partial closures addressed:**
- R-1 is now stated as a narrowing that needs D1-owner confirmation;
- the tier-1 objective and shared capacity are defined;
- the dead-zone cusp is added to the smoothness note;
- the seed = None semantics are stated;
- single-attempt reproduction is a D3 obligation;
- the register is completed;
- the precondition cut-off is population-dependent, with the reason
  nearest-centre stays exact;
- the D1 §16 deviation is recorded.
