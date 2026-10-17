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

**Lifecycle: PROPOSED (D2 draft, not accepted).**
- **Drafted:** 2026-10-17 by an AI agent (Claude Code session
  https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi), under workplan
  `LTA-V7.16-MOLTEN-SALT-OCCLUSION`, SSDP 7.2.1.
- **Parent:** the accepted-current D1
  `doc/specs/v7_16_molten_salt_occlusion_extension.md` (r7, human-ratified
  2026-10-17), cited as **D1 §n**.
- **Status of choices:** every parameter value below is an agent-proposed
  D2 default (§10) until this document is accepted.
- **Acceptance gate:** independent falsification review. The method can
  change which structures are produced and when a request fails, so human
  acknowledgement is also requested (§12).

## 1. Scope and what D2 must preserve

D2 turns the D1 §16 hand-off into a method. These D1 invariants hold
throughout, and nothing below may weaken them:
- **Inventory.** Inventory and charge invariants are exact (D1 §4). Each
  `cage_salt` ion stays in its own cage (S-10).
- **Hard minima.** The D1 §6 minima d_ij ≥ r_cut(i, j) hold on the written
  coordinates, with no tolerance below r_cut. The 1.5 Å catastrophic guard
  also holds.
- **Coordination.** Anchoring holds, together with the per-ion coordination
  rule max(1, ⌈a·C_M(c)⌉) ≤ N_i(Cl) ≤ U_M (D1 §7, S-11, a = 0.4). Under the
  S-12 fallback, the mean rule (A) with a = 0.7 applies instead.
- **Cage charge.** |q_c| ≤ q_max(c) = z_max(c)·e (D1 §4).
- **Ring relocation.** Ring-family limits apply (D1 §8, D-12).
- **Order and streams.** The construction order and stream isolation of D1
  §11 apply. Steps 1-3 (legacy framework stream, v7.15 counterion stream,
  smoothing) are untouched.
- **Fail-closed.** No structure is written unless every check passes.
- **No sampling claim.** No uniform, equilibrium or Boltzmann sampling is
  claimed (D1 §2, D-9).

**Out of scope.** Architecture (modules, data structures, CLI, schema) is D3.
D2 fixes only the numerical identities those layers must honour (§11).

## 2. Definitions

- **Cell and metric.**
  - L is the 3×3 lattice matrix of the realized, scaled cell (rows are
    lattice vectors, Å).
  - Fractional and Cartesian coordinates are related by r = f L.
  - The minimum-image distance is d(r₁, r₂) = ‖(f − round(f)) L‖, with
    f = (r₁ − r₂) L⁻¹.
  - Because max r_min1 = 5.01 Å is less than half the smallest
    perpendicular cell width (about 7.1 Å), a single image suffices for every
    pair query below (verified in D1 review 1).
- **Ions and sets.**
  - S is the set of `cage_salt` ions.
  - R is the set of ring cations (roles `framework` and
    `occlusion_counterion`).
  - F_O is the set of framework O atoms, and F the set of all framework
    atoms.
  - Each ion i has a formal charge z_i and a species.
- **Pair tables.** From the D1-bound E-pair CSV:
  - r_cut(i, j) and r_min1(i, j) are read exactly as D1 §6 and §7 define;
  - the Cl–Cl and Cl–O "around M" rows use the D1 §6 cell-wide weights x_M.
  - ρ_cut(i, j) = r_cut(i, j) + δ_b is the build target, with build margin
    δ_b > 0 (§10). Verification uses r_cut itself.
- **Cage centres.**
  - c ∈ {α₀, α₁, β₀, β₁}, at the D1 §3 canonical fractional positions.
  - home(i) is the cage assigned to ion i by input.
  - margin(i) = min_{c' ≠ home(i)} d(i, c') − d(i, home(i)) is its
    membership margin.
- **Coordination count.**
  - For ions i and j of opposite charge, j is a first neighbour of i iff
    d(i, j) < r_min1(i, j).
  - N_i is the number of first neighbours of i, counted over all ions in the
    cell in any role and any cage.
- **Floor and cap.**
  - L_i is the per-ion coordination floor, and U_i = U_M for species M (D1
    §7).
  - For α `cage_salt` cations: L_i = max(1, ⌈a·C_M(c) − 10⁻⁹⌉). The −10⁻⁹
    guard stops float error from raising an exact integer.
  - For β `cage_salt` cations and for every Cl: L_i = 1 (anchoring), and
    there is no cap.

## 3. Pipeline (concretizes D1 §11)

| Step | Operation | Stream | D2 section |
|---|---|---|---|
| 1 | framework ring placement | legacy `random.Random(seed)` | unchanged |
| 2 | `occlusion_counterion` placement | v7.15 secondary stream | unchanged |
| 3 | cation density smoothing (smoothing entry point only) | existing | unchanged |
| 4 | charge rebalancing | R_q = `v7.16-charge-rebalance:{seed}` | §5 |
| 5 | V_target = f_V(n̄)·V_ref and cell scaling | none (deterministic) | §4.1 |
| 6 | V_eff, n_melt and floors L_i | none (deterministic) | §4.2 |
| 7 | cage-salt build with D-12 relocation | R_s = `v7.16-cage-salt:{seed}` | §6 |
| 8 | exact verification on the coordinates to be written | none | §7 |
| 9 | write | none | D3 |

- **Order of steps 4 and 5.** Step 5 depends only on the inputs, so it may be
  applied while the population is built, as v7.14/v7.15 already do. It must
  precede steps 6-8.
- **Derived streams.** Each derived stream is a NumPy `Generator(PCG64)`
  seeded with the first 8 bytes of SHA-256 of its tag string. This follows
  the v7.15 `v7.15-occlusion:{seed}` pattern. With seed = None, each stream
  is seeded from OS entropy and the run is not reproducible, as in v7.15.
- **Smoothing entry point.** That entry point tries several attempt seeds
  and selects one population. Steps 4-8 run once, on the selected population,
  using that population's attempt seed. This mirrors v7.15
  `finalize_occluded_chlorides`.

## 4. Deterministic quantities

### 4.1 Volume

- n̄ = N_Cl/2 and f_V(n̄) = 1 − 0.03·min(n̄/10.5, 1) (D1 §9).
- An explicit `--volume-factor` replaces f_V.
- V_target = f·V_ref, and isotropic scaling follows v7.14.

### 4.2 Effective volume V_eff,α(c; d) and n_melt

**Estimator.** Deterministic midpoint quadrature on the fractional grid:

```
f_{ijk} = ((i + 1/2)/G, (j + 1/2)/G, (k + 1/2)/G),   i, j, k = 0..G-1
V_eff(c; d) = (V_cell / G^3) * #{ f_ijk : nearest_centre(f_ijk) = c  and  min_{o in F_O} d(f_ijk, o) > d }
```

- **Grid.** G is the smallest integer with max_k |a_k|/G ≤ h, where |a_k| are
  the lattice-vector lengths and h = 0.10 Å.
- **Offset.** The half-step offset keeps grid points off the partition
  bisector planes of the canonical centres. On the fixture, review 3 found
  that an aligned grid biases V by −1.2 %.
- **Ties.** A grid point equidistant from two centres within 10⁻⁹ Å is
  assigned to the cage first in the order α₀, α₁, β₀, β₁.

**Error.** The indicator-function quadrature error scales with the surface
area times h. The required accuracy is |V(G) − V(2G)|/V ≤ 0.5 % on the
canonical fixture, checked at qualification (§9).
- A 0.5 % volume error moves a·C by at most 0.5 %. That can flip
  ⌈a·C⌉ only when a·C is within 0.5 % of an integer.
- Such cases are deterministic, so they are reproducible, and they are
  reported with their margin (§8).

**Dependent quantities.**
- d_wall(c), x_Li(c), v_pair, n_melt(c), the reported range (± 0.3 Å),
  φ(c), C_M(c) and L_i follow exactly as defined in D1 §7.
- They are computed once per structure, after step 5 and before step 7.
- They depend only on the inventory and the scaled cell, so the build cannot
  change them.

### 4.3 Cage charge observable

q_c^exc, q_c^fw, the weights w(j, c) and z_max(c) are as defined in D1 §4.
They are computed in exact rational arithmetic, as halves of integers, so
the comparisons |q_c| ≤ q_max(c) involve no floating-point rounding.

## 5. Charge rebalancing (step 4): exact, not heuristic

**Classes.** Ring sites fall into *window classes* by their weight pattern
w(·, c). On the canonical template these are:
- 4R@α₀ and 4R@α₁;
- 6R@(β_k, α_l), four classes;
- 8R@(α₀, α₁).

Every site in a class contributes identically to every q_c.

**Exact class-count search.** A counterion's ring family and species are
fixed by input. The decision is therefore an integer vector: for each
family, species and class, how many counterions sit there.
- The capacity of a class is its number of sites that are empty after step 3
  or already hold an `occlusion_counterion`.
- The search space is the product over families of the compositions of each
  species count into that family's classes: at most a few times 10⁴ vectors.
- It is enumerated exhaustively.

Among feasible vectors with every |q_c| ≤ q_max(c), choose lexicographically
the one that minimises:
1. the number of counterions that change class from their step-3 placement;
2. Σ_c q_c²;
3. the lexicographically smallest vector in a fixed canonical order.

If no vector is feasible, this is a *seed-level failure*: the run reports the
best achievable max_c(|q_c| − q_max(c)) and the blocking framework sites.

**Up-front check.** This is the same search with step-3 occupancy ignored
(every site of each family available) and q_c^fw ≡ 0. It is the deliberately
stricter input rule of D1 §4, run before any placement.

**Realization.**
- Counterions already in a class that keeps them stay where they are.
- Each counterion leaving an over-full class moves to a site of a
  deficient class of the same family. Both the counterion and the target
  site are chosen with R_q, uniformly among the empty sites of that class.
- Departing counterions are taken in canonical site-id order.
- This changes smoothed positions; the moves are reported (D1 §4).

## 6. Cage-salt build (step 7)

### 6.1 Variables and objective

- **Variables.** The Cartesian positions x ∈ ℝ^{3|S|} of the `cage_salt`
  ions. Ring cations and framework atoms are fixed except through discrete
  D-12 moves (§6.4).
- **Penalty function.** h(u) = max(0, u)².
- **Smooth count.** s_i(x) = Σ_{j opposite charge} σ((r_min1(i, j) − d_ij)/w),
  where σ is the logistic function and w = 0.12 Å.

The objective at weight W is:

```
E_W(x) = W * [ sum_{i<j in S} h(rho_cut(i,j) - d_ij)
             + sum_{i in S, k in R} h(rho_cut(i,k) - d_ik)
             + sum_{i in S, o in F_O} h(rho_cut(i,o) - d_io)
             + sum_{i in S, k in F} h(1.5 + delta_b - d_ik)
             + sum_{i in S} h(eps_m - margin(i)) ]
         + sum_{i in S, cation} [ h(L_i + mu - s_i) + h(s_i - (U_i - mu)) ]
         + sum_{i in S, Cl}      h(1 + mu - s_i)
```

- The hard terms are exact-zero at feasibility and once differentiable.
- The count terms push the smooth count beyond the discrete threshold by a
  margin μ, so the discrete count passes.
- Ring-cation and framework positions enter only as fixed parameters.

### 6.2 Initialization (one start)

- **Order.** Ions are inserted in an order drawn from R_s.
- **Candidates.** For each ion, positions are drawn uniformly in the cell from
  R_s. A candidate is accepted when:
  - it lies in home(i), with margin ≥ ε_m;
  - it is at least ρ_cut from every framework O, every ring cation and
    every salt ion already placed.
- **Fallback.** After M_ins = 2000 rejected candidates, the candidate with the
  smallest maximum violation is kept. Relaxation (§6.3) then removes the
  overlap.
- **Purpose.** This is random sequential insertion with a hard core. It is
  an initial condition, not a sampling law.

### 6.3 Relaxation and coordination shaping (continuation)

- **Minimizer.** L-BFGS-B minimises E_W from the start, with
  W = 10, 100 and 1000 in turn. Each stage starts from the previous result,
  with at most 6000 iterations per stage.
- **Gradients** are analytic.
- **No tolerance trust.** The minimiser's convergence flag is not used.
  Acceptance is decided only by the exact verification (§7).

### 6.4 D-12 relocation (discrete, within the start)

If verification fails on hard minima involving a ring cation k, one
relocation round runs. In descending order of violation, each conflicting
ring cation is offered the candidate sites allowed by D1 §8:
- empty sites of its own family;
- other variants of its own physical site (O1a).

A candidate is accepted only if all three hold:
1. it strictly reduces k's maximum violation against S;
2. it keeps every |q_c| ≤ q_max(c) (§4.3, recomputed exactly);
3. it is the best such candidate. Ties are broken by canonical
   (site_id, variant_id) order, then R_s.

After a round with at least one move, §6.3 is re-run from the current x. At
most R_reloc = 4 rounds are allowed per start.

### 6.5 Restarts and fail-closed

- Starts s = 1..K_max with K_max = 20 are run in order.
- The first start whose structure passes §7 is returned. This is
  *first-success*: later starts are never run, so the output does not depend
  on K_max beyond s.
- If no start passes, the run fails closed. It writes no structure and
  reports, per start, the best slack per constraint class and the CN
  deficits.

### 6.6 Fallback rule A (S-12, only if adopted at qualification)

The count terms of §6.1 are replaced by a smooth version of the per-cage,
per-species mean condition. Here μ is applied to the mean, the pooling of D1
§7 is used, and anchoring and U_M stay per-ion. Everything else is
unchanged.

## 7. Exact verification (step 8)

Verification runs on the *coordinates that will be written*.
- **Canonicalisation.** Fractional coordinates are rounded to 10 decimals,
  wrapped into [0, 1), and every writer emits these. Re-checking from a
  writer's output must reproduce the same result (D4 test).
- **Pairs.** All pairwise distances are computed by brute force in float64
  under §2's minimum image, with no neighbour-list approximation. With at
  most a few hundred atoms this costs O(10⁵) pairs.

| # | Check (all must pass) | Tolerance |
|---|---|---|
| V1 | Per-cage, per-role, per-species inventory, by strict nearest centre (margin > 0) | exact |
| V2 | d_ij ≥ r_cut(i, j) for every pair with at least one `cage_salt` ion, including ring cations and framework O; d ≥ 1.5 Å to every atom | none below r_cut |
| V3 | Anchoring and L_i ≤ N_i ≤ U_i (or rule A if adopted) | exact integer counts |
| V4 | Every |q_c| ≤ q_max(c) | exact rational |
| V5 | Q_fw = 24 and Q_exc = N_Cl | exact |

- **Precision.** Distance values near r_cut are compared in float64.
  Rounding to 10 decimals moves a position by at most 10⁻¹⁰ in fractional
  units, about 2×10⁻⁹ Å. That is far below δ_b, so the build margin, not
  verification, absorbs it.
- **Coordination counts** use the strict inequality d < r_min1 on the same
  canonical coordinates.

## 8. Stochastic, error and reproducibility semantics

- **What is random.** The insertion order and positions (R_s), the
  relocation tie-breaks (R_s), and the choice of site within a class during
  rebalancing (R_q). The output is a deterministic function of the inputs,
  the seed and the numerical stack.
- **No sampling law is claimed.** The distribution of outputs over seeds is
  whatever this procedure induces. It is biased toward configurations with
  pairs near ρ_cut and counts near L_i + μ, because penalty relaxation
  settles on constraint boundaries. That is accepted (D-9: MD melts the
  structure) and made visible (§8.1).
- **Reproducibility envelope.**
  - Same inputs, seed, and pinned NumPy/SciPy/BLAS versions on the same
    platform give byte-identical output.
  - Across numerical stacks, L-BFGS-B trajectories may diverge, so outputs
    may differ. Every output still satisfies §7, and the qualification grid
    (§9) is re-run when the stack changes.
  - The v7.15 legacy model and the dry path do not use this machinery, and
    keep their existing byte-identity guarantees.
- **Numerical error sources.**

  | Source | Bound | Effect |
  |---|---|---|
  | V_eff quadrature | ≤ 0.5 % | can shift a floor L_i only near integer thresholds; reported |
  | float64 distances | negligible against δ_b | none |
  | Minimiser non-convergence | — | caught by §7: it can only cause a restart or a failure, never a wrong output |

  The D1 parameter uncertainties (a, d_wall, f_V, r_cut) are D1's and are not
  propagated here.

### 8.1 Numerical inspectability (O1)

**Material realized record per run**, written to the run report that D1 §14
accepts:
- **Rebalancing:** the class-count vector chosen in step 4, the moves made,
  and the alternatives' objective values.
- **Effective volume:** G and h, and V_eff at d_wall and ± 0.3 Å.
- **Per start:**
  - E_W after each stage;
  - the best slack per constraint class (salt–salt, salt–ring, salt–O,
    guard, membership);
  - the CN deficit and excess per ion class;
  - the relocation moves;
  - pass/fail with the failing check.
- **Final structure:**
  - the V1-V5 table;
  - histograms of d − r_cut per pair class over the first 0.5 Å;
  - histograms of N_i − L_i per species;
  - the starts used and K_max.

**Reader and questions.** The stakeholder and qualifier ask:
- Did it pass, and on which start?
- How close to the bounds did it sit?
- Was the solution boundary-hugging?
- Did relocation or rebalancing move ring cations?
- How sensitive are the floors to V_eff?

**Product surface.** This is the D1 §14 per-run report, already accepted at
D1 ratification. D2 adds only the fields above to it.

## 9. Evidence and oracles (to be produced in D4; defined here first)

**Oracles independent of the production code.**
1. **Volume quadrature.** On a manufactured cubic cell with a single "O" at
   the origin, V_eff of the complement of a sphere is known in closed form.
   The error must be ≤ 0.5 % at h = 0.10 Å, and refinement h → h/2 must
   converge.
2. **Brute-force verifier.** An independent O(N²) verifier, written
   separately from the production checks, recomputes V1-V5 on written files
   (CIF, POSCAR and LAMMPS). The two verdicts must agree.
3. **Manufactured feasibility.** A hand-built salt configuration known to
   satisfy V1-V5 must pass verification. Perturbing it to violate each check
   alone must be rejected by exactly that check (one test per check).
4. **Rebalancing.** On small synthetic class problems, the exhaustive search
   is compared with brute force over individual site assignments.
   Σ_c q_c = 0 always holds.
5. **Streams.** For a fixed seed, steps 1-3 are byte-identical with and
   without the cage-salt model. Changing R_s leaves steps 1-6 unchanged.
6. **Determinism.** Two runs with the same seed and stack give identical
   bytes.

**S-12 qualification protocol (pre-specified).**
- **Grid.** Na-A and K-A rings; seeds 0-4; eutectic loadings of 2, 4, 6, 8,
  10 and 12 Cl per α, each β holding 1 Cl and 1 Li; plus the three P5 dopant
  compositions (Li5K4Cs1, Li5K3Sr1, Li6K4Ba1 per α) at seeds 0-4. That is
  2 × 5 × (6 + 3) = 90 requests.
- **Per request.** Pass means a structure that passes §7 within K_max = 20.
  Each request is deterministic, so the result is a census of the grid, not
  a statistical estimate.
- **Decision rule, fixed now.**
  - Rule B is retained iff ≤ 4 of the 90 requests fail (pass rate ≥ 95 %).
  - Otherwise the release default becomes rule A with a = 0.7, which must
    then pass the same criterion. If it does not, that reopens D1.
  - The record states both pass rates, and the Wilson 95 % interval for
    information.
- **V_eff refinement check.** Run at G and 2G on the fixture, as in oracle 1.

**Prior exploratory evidence.** Probes P5-P8 (D1 §12) used an earlier
variant of this method: insertion, L-BFGS penalty with CN continuation, and
greedy relocation. Shaped runs reached the per-ion floors in P7 (18/18) and
P8 (26/30), and all four P8 misses were anchoring failures. These probes
motivate the method but are not acceptance evidence. The production method
differs in four ways:
- explicit anchoring and membership terms;
- the build margin δ_b;
- charge-constrained relocation;
- exact verification.

## 10. Parameter register (agent-proposed D2 defaults)

| Parameter | Value | Origin | Basis |
|---|---|---|---|
| δ_b build margin | 0.01 Å | agent | ≫ float error, ≪ MD-relevant length |
| ε_m membership margin | 0.05 Å | agent | keeps membership strict after rounding |
| w, smooth-count width | 0.12 Å | agent (P6-P8) | sharp enough to count, smooth enough for L-BFGS |
| μ, count margin | 0.3 | agent | pushes smooth count past the integer threshold |
| W continuation | 10, 100, 1000 | agent (P6-P8) | penalty continuation |
| L-BFGS-B iterations per stage | 6000 | agent (P6-P8) | — |
| M_ins insertion tries | 2000 | agent | — |
| K_max starts | 20 | agent | 4 to 5 × the P6-P8 start counts |
| R_reloc rounds per start | 4 | agent (P5) | — |
| h, quadrature spacing | 0.10 Å | agent | error target 0.5 % (§4.2) |
| Coordinate canonicalisation | 10 decimals | agent | §7 |
| Stream tags | `v7.16-charge-rebalance:{seed}`, `v7.16-cage-salt:{seed}` | agent | follows v7.15 |
| S-12 grid and threshold | §9 | stakeholder (fallback and 95 % idea) and agent (grid) | D1 §7 |

These are governed defaults. Changing them is a D2 change unless D3 exposes
one as a user option with this default and a recorded value per run. Only
K_max and the seed are proposed as user options.

## 11. D2 → D3 hand-off

D3 must preserve:
- the step order and stream derivation of §3, including that steps 4-8 run
  once on the population the smoothing entry point selected;
- exact rational charge arithmetic and the exhaustive class-count search (§5);
- the V_eff estimator definition, including the offset grid and tie order
  (§4.2);
- the objective and continuation (§6.1-6.3), relocation acceptance (§6.4) and
  first-success restarts (§6.5);
- verification on the canonical written coordinates by brute force (§7),
  with all writers emitting those coordinates;
- the run-report fields of §8.1;
- SciPy as a declared runtime dependency, for L-BFGS-B. NumPy `Generator`
  PCG64 is used for the derived streams;
- the reproducibility envelope (§8) recorded in qualification (stack
  versions).

D3 owns modules, data structures, neighbour handling outside verification,
vectorisation, the CLI and schema, error types and report format.

## 12. Acceptance request and open points

- **Requested:** independent falsification review of this draft. Then human
  acknowledgement of the method-level choices that can change outputs or
  failure rates: first-success restarts, K_max = 20, the boundary-hugging
  bias (§8), and the S-12 decision rule (§9).
- **Tension search.** Searched the D1 spec §13 (T1-T22), the workplan and
  the repository issues (none). Tensions bearing on D2:
  - T8, unshaped packing isolates ions: addressed by the anchoring terms;
  - T17 and T21, V_eff uncertainty: handled by the deterministic estimator,
    with the range reported;
  - T22, resolved by S-11/S-12: the fallback is implemented in §6.6.
  None is closed by this document.
- **Known risks.**
  - The per-ion floor may be hard to meet for dopants (not probed; D1 §12).
  - First-success makes the output depend on the order of starts. That is
    deterministic, but it means a later start is never "better".
