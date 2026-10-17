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

# v7.16 variable-loading molten-salt occlusion extension

**Lifecycle: PROPOSED (D1 draft, not accepted).** Drafted 2026-10-17 by an AI
agent (Claude Code session
https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi) under workplan
`LTA-V7.16-MOLTEN-SALT-OCCLUSION`, SSDP 7.2.1. Nothing here binds D2-D4 until
independent falsification review and stakeholder ratification are complete.
Revision r1 (same day) folds in independent falsification review 1 (§17): it
raised two SERIOUS CHALLENGES against this draft, both now carried as
stakeholder decisions or fixes. Revision r2 (same day) records the
stakeholder's answers. Where the stakeholder delegated, it records the agent
recommendations adopted (§15). Revision r3 records stakeholder decisions S-4
to S-6:
- S-4: "dense" is defined by melt density at 800 K;
- S-5: a charge-balance criterion with rebalancing;
- S-6: the contraction form.
Revision r4 folds in independent falsification review 2 (§17). It changes
the calibration route for the charge threshold, adds framework cations to the
charge observable, and corrects the rationale for the contraction form. Open
items are marked **[CONFIRM]**; they are listed in §15. The document stays
PROPOSED until the stakeholder ratifies it as a whole.
Items marked **[DECIDE]** are open stakeholder decisions; items marked
**(proposed)** are agent-chosen defaults that become binding only if ratified.
Stakeholder decisions D-1 to D-13 are cited from the workplan §2; they are
instructions, and this document is where they would become ratified D1.

**Predecessor.** `doc/specs/v7_15_chloride_occlusion_extension.md`
(accepted-current, human-ratified) at base commit
`4228470e274afb3d55ca321f29c1ce9a0504b34d`. Three of its own D1 reopen
triggers (archived workplan §11) fire: dense loading / multiple guests per cage;
counterions as free cage guests; and evidence that "dry volume + explicit
factor" is unsuitable (salt-loaded LTA contracts [R18]).

## 1. Scope and authority map

On ratification this document:

- **adds** a new default occlusion model, the *cage-salt model* (§3-§9);
- **retains** the v7.15 sparse model unchanged as an explicitly selected
  *legacy model* (§10);
- **supersedes, for the cage-salt model only,** these v7.15 sections:
  "Sparse chloride domain and canonical setting" (Cl domain and siting),
  "Collision guard", "Dry-volume compatibility boundary", "Smoothing and
  assembly" (step order) and "Epistemic scope of the literature" (extended);
- **generalises** the v7.15 "Independent charge invariants" (§4) and keeps the
  v7.15 Al24 scope, canonical-template gating, maximum-symmetry and
  unique-batch rejection, RNG isolation and single assembly path.

Retained authorities, with boundary notes needed on ratification:
- `v7_13_group_i_ii_cation_extension.md` keeps the supported cations, formal
  charges and Q_fw = +24. Boundary: when the cage-salt model is *active* (§3),
  §5 restricts the cations in the structure to S.
- `v7_14_group_i_ii_volume_reference_extension.md` keeps the dry reference
  volume V_ref. Boundary: when the cage-salt model is active, the default
  volume factor becomes f_V (§9).

The dry (no-occlusion) generator path is unaffected and must remain
output-identical.

Non-authoritative evidence bound by this document (stakeholder-supplied
package rev 3, 2026-10-17, `doc/evidence/`):

| Tag | File | Role here |
|---|---|---|
| E-occ | `lta_salt_occlusion_evidence.md` (rev 2.1) | loading, siting, cation effects [Rn] |
| E-pair | `molten_salt_pair_distances.csv` / `.md` (rev 2.1, 91 rows) [Pn] | hard minima, first-shell boundaries |
| E-cn | `lta_coordination_numbers.csv` / `.md` (rev 3.0, 215 rows) [Qn] | coordination targets |
| E-vol | `lta_salt_loading_lattice.csv` / `lta_salt_loading_lattice_contraction.md` (rev 3.0) | volume contraction |

Bound identities (SHA-256, as supplied 2026-10-17; re-confirmed at
ratification):

| File | SHA-256 |
|---|---|
| `molten_salt_pair_distances.csv` | `a59f08eb7192caa7051c34e6cebfcacdd2245df94bc7eef20b8f6901398a27bd` |
| `lta_coordination_numbers.csv` | `8ac1d074b41366075450088c0b92496e6889b297ca8b557a2a2ee26d1c36cfbc` |
| `lta_salt_loading_lattice.csv` | `a95c54b8e2f22a41b971121af05217a314de419f8d742a0cdea4dd63a9a21496` |

Both generator scripts regenerate their CSVs byte-identically (agent check,
2026-10-17, run in an isolated directory).

## 2. Scientific question and epistemic scope

**Question.** Produce *initial configurations* for AIMD/MLFF molecular dynamics
of dehydrated Al24 zeolite A carrying an occluded chloride-salt phase at a
user-specified per-cage loading and composition, such that MD can melt and
equilibrate them. The generator's product is the configuration; it estimates
no thermodynamic quantity.

**What a configuration claims.** Exact inventory and charge (§4); no ion pair
closer than the bound melt onset distance (§6); every ion anchored to an
opposite-charge first neighbour and per-species first-shell coordination
approximately at bound melt values where the loading makes that meaningful
(§7); a cell volume from the bound contraction model (§9). Ring cations are
not coordination centres. "Amorphous" means only: positions are random subject
to these constraints (D-9).

**What it does not claim.** Equilibrium siting, g(r), occupancy statistics,
uniform sampling of configuration space, energetic stability, or any
equilibrium loading. β-cage Cl remains a training hypothesis (as in v7.15;
the β default of one halide is inference [E-occ §4]). Literature values below
are melt or crystal data transferred to a confined phase by assumption.

**Validity regime.** Al24Si24O96 primitive cell of the canonical Fm-3c template
(v7.15 gating); base salt LiCl–KCl (the only base salt the loading evidence
supports, 773–823 K [R1, R2, R15]); dopants from the species set S (§5);
0 ≤ Cl ≤ 12 per α-cage, 0–1 per β-cage. Pair data refer to 950–1300 K melts;
CN data to 1100 K or T_m. The contraction working value rests on 294 K
crystals at 4 Br per pseudo-cell; only its bound comes from in-situ
650–800 K data. The
generator does not take a temperature; these are interpretation limits.

## 3. Definitions

- **Cell, cages.** The primitive cell contains α-cages α₀, α₁ at fractional
  (1/4,1/4,1/4), (3/4,3/4,3/4) and β-cages β₀, β₁ at (0,0,0), (1/2,1/2,1/2)
  (v7.15 canonical centres). One pseudo-cell = one α + one β; literature
  loadings are per pseudo-cell (D-1).
- **Ring sites.** Physical sites of families R ∈ {4, 6, 8} with capacities
  24 / 16 / 6 and their template positional variants (v7.15).
- **Roles.** `framework` (ring site, compensates the framework);
  `occlusion_counterion` (ring site, compensates Cl; v7.15 role); `cage_salt`
  (new: an ion of the salt phase assigned to a cage). Cl⁻ is always
  `cage_salt`.
- **Active.** The cage-salt model is *active* iff any cage-salt input is
  nonzero: some N_Cl(c) > 0 or some N_i^salt(c) > 0. Otherwise every
  cage-salt rule (§5-§9) is inert and the dry or legacy path applies
  unchanged.
- **Cage assignment (adopted, O8, §15).** An ion belongs to the cage whose centre
  is nearest under the minimum-image Cartesian metric (nearest-centre
  partition). Ties are a construction failure.
  - Geometry note: the α–β centre distance is 10.63 Å, so the bisector lies
    5.32 Å from the α centre, about 1.1 Å on the α side of the 6R sites. An
    ion on the α face of a 6R is therefore counted as β. This affects only
    accounting, not distances.
  - Not adopted: a sphere of radius 5.7 Å about each α centre (11.4 Å free
    diameter [R31]). It lacks a β region and a rule for the 5.7–6.1 Å window
    shell.
  - Feasibility in P5 was insensitive to the choice, but only 24 points
    tested the partition (§12).
- **First neighbour.** For ions i, j of opposite charge, j is a first
  neighbour of i iff d_ij < r_min1(i, j) (E-pair column `r_min1`), minimum
  image, both ions in any role.
- **Coordination number.** N_i = number of first neighbours of i. Only
  opposite-charge ions count (D-13); framework O is not an ion of the table and
  does not count.

## 4. Inventory and charge invariants

**Input.** Per cage c ∈ {α₀, α₁, β₀, β₁}: N_Cl(c) and N_i^salt(c) for each
cation species i; ring-family counts for `framework` and `occlusion_counterion`
roles (exact integers, v7.13 typing rules).

**Domain.**

```
0 <= N_Cl(alpha_k) <= 12                  (D-2, hard cap; never lowered)
N_Cl(beta_k) in {0, 1}                    (D-2)
sum_i N_i^salt(beta_k) <= N_Cl(beta_k)    (adopted; D-6 "MCl unit")
```

**Invariants (independent, never collapsed).**

```
Q_fw  = sum_i z_i N_i^fw = 24                                    (v7.13)
N_Cl  = sum_c N_Cl(c)
Q_exc = sum_i z_i ( N_i^occ + sum_c N_i^salt(c) ) = N_Cl          (generalises v7.15)
```

A β-cage Cl may be compensated by a cation in the β-cage or on a ring site
(D-6).

**Cage charge balance (stakeholder decision S-5, 2026-10-17).** Charge
separation must not be too severe. The code computes the charge distribution
and rebalances regions where charge is unreasonably concentrated. The
threshold "could be inspired by" studies of charge fluctuations in molten
salts.

- **Observable (proposed).** The local charge deviation of cage c is
  q_c = q_c^exc + q_c^fw:

  ```
  q_c^exc = sum_{i in cage_salt(c)} z_i - N_Cl(c) + sum_{j in occlusion_counterion} z_j w(j, c)
  q_c^fw  = sum_{j in framework} z_j w(j, c) - <same quantity averaged over cages of c's type>
  ```

  - **Sharing weights.** These follow the window topology on the canonical
    template (verified independently in review 2):
    - a 4R site belongs to one α: w = 1 (5.0 Å from that α centre, 7.16 Å
      from the nearest β);
    - a 6R site separates one β and one α: w = ½ each (4.24 Å and 6.40 Å);
    - an 8R site separates α₀ and α₁: w = ½ each (6.27 Å from both).
  - **Framework term.** q_c^fw measures how unevenly the framework-compensating
    cations are spread relative to the same-type average. The framework anion
    shares are equal for same-type cages by symmetry. Example: on the fixture
    (Na:2/16/6, seed 0), both 4R framework Na sit in α₀, so q^fw = ±1 e
    between α₀ and α₁ (review 2).
  - Σ_c q_c = 0 for every structure.
- **Criterion.** For every cage, |q_c| ≤ q_max. It is checked on the final
  written structure (§11).
- **Rebalancing (proposed mechanism; §11 step 4).**
  - After smoothing and before the salt build, `occlusion_counterion` cations
    are re-sited to empty sites of their own ring family, to bring every
    |q_c| ≤ q_max. This is the D-12 move set.
  - Later D-12 relocations of any ring cation are accepted only if they keep
    every |q_c| ≤ q_max.
  - Per-cage salt inventories are input (D-1) and are never changed
    **[CONFIRM]**. Moving salt ions between cages would require relaxing D-1.
- **Failure handling.**
  - *Up-front rejection* is a capacity bound. The input is rejected, naming
    the cages, if no assignment of the requested counterions to any sites of
    their requested ring families (best case) reaches |q_c^exc| ≤ q_max. Note
    that 8R moves cannot shift charge between α₀ and α₁, because every 8R
    site joins both.
  - *Seed-level failure* is a construction failure for that seed. It occurs
    when framework occupancy, which depends on the seed through the legacy
    stream, blocks the needed sites. Example from review 2: α₀ with 12 Cl,
    no cage cations and 12 counterions on 4R needs all 12 of α₀'s 4R sites,
    but on seed 0 two of them hold framework Na. The report names the blocking
    sites.
- **Threshold q_max (stakeholder-set; provisional 1 e, [CONFIRM]).**
  - **Literature (context only; agent search, 2026-10-17; abstracts and
    snippets, full texts blocked by this environment's network policy).**
    Molten salts screen charge perfectly: the charge–charge structure factor
    satisfies S_ZZ(k) → 0 as k → 0 (zeroth moment), and the Stillinger–Lovett
    condition fixes its k² coefficient (second moment). Molten NaCl is
    discussed in [P4: Zeidler et al., *J. Chem. Phys.* 157, 094504 (2022)].
    For large sub-volumes, the variance of the net charge then grows with the
    boundary area, not the volume (charge hyperuniformity):
    - Ph. A. Martin, T. Yalcin, *J. Stat. Phys.* 22, 435 (1980);
    - J. L. Lebowitz, *Phys. Rev. A* 27, 1491 (1983);
    - T. Hoang Ngoc Minh, B. Rotenberg, S. Marbach, *Faraday Discuss.* 246,
      225 (2023), by simulation.

    No study gives a number at cage scale.
  - **Why a bulk-melt sphere statistic is not used (review 2, SERIOUS
    CHALLENGE SC-A).** Counting ion centres in a cage-sized bulk sphere
    (R ≈ 5.9 Å) is dominated by ions cut by the arbitrary boundary. Even a
    perfect rock-salt lattice at melt density gives σ_Q ≈ 1.35 e and |Q| up
    to 4 (review-2 probe). A rule ⌈3σ_Q⌉ would therefore allow about 5 e and
    fail S-5. A cage boundary is instead mostly an impenetrable wall with
    fixed framework charge, so the two statistics are different estimands.
  - **Calibration route (proposed, optional).** Run MD of salt-loaded LTA at
    800 K with the stakeholder's MLFF. Compute the distribution of q_c, with
    exactly this observable, partition and weights. Set q_max to a stated
    upper percentile of |q_c|, for example the 99th.
  - **Provisional value.** Until calibrated, q_max = 1 e is an agent
    placeholder for the stakeholder to accept or replace. It allows a single
    β-Cl counterion on a 4R site (q_β = −1) and the fixture's framework
    imbalance of ±1 e. The E-occ statement "each Cl⁻ arrives with a
    charge-balancing cation" concerns inventory, not locality, and is not
    cited as support.
- **Reported:** every q_c^exc and q_c^fw, the rebalancing and relocation moves
  with their charge effect, and the distribution of q_c against q_max.

Inputs that are cation-only in a cage (N_i^salt(c) > 0, N_Cl(c) = 0) are
rejected up front (adopted), because those cations could never be anchored.
Ring-site capacity sum_i (N_ri^fw + N_ri^occ) ≤ C_r holds as in v7.15.

The final structure is a construction failure unless every role, cage and
species count equals its input exactly, with cage membership evaluated on the
written coordinates.

## 5. Supported salt species

S = {Li, Na, K, Rb, Cs, Sr, Ba}: the v7.13 cations for which E-pair supplies
every needed pair. When the cage-salt model is active, **every cation in the
structure, in any role, must be in S**, because every cation is a potential
neighbour of a cage Cl. Any other species (H, Fr, Be, Mg, Ca, Ra) is rejected
naming the missing pair; there is no radius fallback (D-3).

Consequence (adopted, §15): Ca-containing frameworks (zeolite 5A) and Mg
forms cannot be salt-loaded until a pair table for them is supplied. E-pair
§2.2 holds measured Mg–Cl 2.42 Å and Ca–Cl 2.78 Å peaks, but no recommended
rows, and D-3 forbids deriving them here.

Trivalent cations (U³⁺, Ln³⁺, Y³⁺) are outside v7.13 and out of scope (T5).
Free oxide O²⁻ is not a supported species; the O-bridged and O²⁻ rows of E-pair
are not used.

## 6. Hard minimum distances

For every pair of distinct ions (a, b) in the written structure, with at least
one in role `cage_salt`, the minimum-image distance must satisfy
d_ab ≥ r_cut(a, b), with r_cut from E-pair:

| Pair class | Row used |
|---|---|
| cation–Cl | (M, Cl) |
| cation–cation | (M, M′) "Cl-bridged" |
| Cl–Cl | "around M" rows, mole-fraction weighted: r_cut = Σ_M x_M r_cut(Cl–Cl around M) |

- **Weights x_M (adopted, §15).** Mole fractions over all `cage_salt`
  cations of the cell; if there are none, over all cations. The values range
  only from 2.97 Å (Li/Na/Sr) to 3.49 Å (Cs), and MD relaxes the result (D-9).
- **Ring cations get the same table values as any other ion** (D-12). No
  reduced or scaled cut-off exists in this model. The v7.16 draft's earlier
  scaling factor s (workplan O2) is retired, because the table's r_cut column
  is already a closest-approach onset (E-pair §0, §4.3).
- **Ring–ring pairs** (both on ring sites) keep the template geometry and are
  not constrained by this document. Their slack against the Cl-bridged r_cut is
  reported.
- **Ion–framework minima (stakeholder decision S-1, 2026-10-17: option (b)).**
  E-pair §6.4 defers cation–framework contacts to framework-oxygen values,
  which the evidence does not supply. The evidence holds only crystal
  Ag–O 2.49–2.52 Å [R27, R28] and Sr–O3 2.21–2.31 Å [Q25]. The stakeholder
  accepts the E-pair free-oxide rows as an explicit grade-E approximation:
  - cation–O(framework) uses the (M, O) r_cut;
  - Cl–O(framework) uses the "Cl–O around M" r_cut, weighted as for Cl–Cl
    below.
  Check: Sr–O 2.09 Å ≤ crystal Sr–O3 2.21 Å. The rejected alternatives were a
  supplied table and guard-only. Guard-only would have been scientifically
  inadmissible, since it permits Cl⁻ or Li⁺ 1.5 Å from framework O.
  Framework Si/Al: no data. Proposed: no separate T-atom minimum beyond the
  catastrophic guard below. (The probes used 0.85 × Shannon sums for Si/Al.
  That is a radius formula outside D-3 and is disclosed in §12.)
- **Reading of D-3.** "Same-sign cut-off larger than opposite-sign" holds per
  pair member in the table: every like pair's r_cut exceeds both members' r_cut
  against Cl. It does not hold globally (for example Li–Sr 2.65 < Cs–Cl
  2.78 Å). This document adopts the per-member reading.
- **Catastrophic guard (retained from v7.15).** Every `cage_salt` ion is at
  least 1.5 Å from every atom. This is not a contact distance.

A structure violating any hard minimum is never written (fail closed).

## 7. Anchoring and coordination targets

**Anchoring (hard, proposed).** Every Cl has ≥ 1 first-neighbour cation in any
role, and every `cage_salt` cation has ≥ 1 first-neighbour Cl. Basis: halides
in LTA are anchored to cations and never seen isolated [E-occ §4]; D-5 relies
on coordination in place of a binding rule. Probe P5 showed that unconstrained
random packing leaves isolated cations (CN 0), so this needs to be stated.

**Coordination target (cation-centred, α only).** For each species M with
n_M ≥ 1 `cage_salt` centres *in α-cages*, the mean of N_i(Cl) over those centres
must satisfy the coordination rule below (density-reduced, per S-2/S-4).
β-cage cations are excluded. Under D-2 they
have at most one Cl in their own cage (CN of about 1, as measured in review 1).
Their CN is reported separately and subject only to anchoring. Ring cations
count as neighbours of Cl (O4c, proposed yes), but are not CN centres. Only
cation-centred targets are imposed (D-13). Cl-centred counts are reported, not
targeted, because:
- by number conservation they follow from the cation targets plus the ring
  contribution;
- with ring cations counted they exceed bulk values (finding T7).

**Density reference (stakeholder decisions S-2 and S-4, 2026-10-17).**
- S-2: "Coordination is only for the dense case. The expected coordination
  should be reduced in relationship to the density, though the exact
  relationship I am unsure." 
- S-4: "dense" corresponds to the density of the molten salt at 800 K, the
  typical operating temperature of LiCl–KCl with LTA.

**Melt number density at 800 K (grade C/D).**
- **Correlations.** Janz pure-salt correlations, as cited in later papers:
  ρ_LiCl = 1.8842 − 4.33×10⁻⁴ T and ρ_KCl = 2.1359 − 5.83×10⁻⁴ T
  (g cm⁻³, T in K).
- **Mixing.** Ideal mixing of molar volumes reproduces the Janz eutectic
  point (58:42, 773 K): calculated 1.6227 against 1.6226 g cm⁻³ as cited,
  for example by *J. Nucl. Fuel Cycle Waste Technol.* 15, 117 (2017)
  (e-sciencecentral.org/upload/jnfcwt/pdf/JNFCWT-15-117.pdf) and by the
  density study of eutectic LiCl–KCl with UCl₃/CeCl₃/LaCl₃ on ResearchGate
  (publication 320063307). The originals were not opened.
- **At 800 K.** The eutectic gives ρ = 1.609 g cm⁻³, about 57.7 Å³ per ion
  pair (57.4 Å³ at x_Li = 0.592). For a cage, x_Li is the Li:(Li+K) ratio
  of its `cage_salt` cations. If the cage has none (only dopant cations, or Cl
  compensated by ring counterions), the [R1] eutectic x_Li = 0.592 is used
  (proposed).
- **Dopants.** Their volume is not modelled (dilute assumption), because no
  sourced density is bound for them. If dopants exceed 20 % of a cage's
  `cage_salt` cations, the report flags that cage's density reference as
  outside validity (proposed). CsCl has about 100 Å³ per pair, so the
  assumption matters for Cs-rich cages.
- **Extrapolation.** Both pure-salt correlations are used far below the pure
  melting points (LiCl 883 K, KCl 1044 K). The 10⁻⁴ agreement with the cited
  eutectic value may be circular if that value was itself computed by
  additivity (T19).

**Cage capacity at that density.** n_melt(c) = V_free,α(c) / v_pair(800 K,
x_Li(c)). The proposed V_free,α is the cage's nearest-centre cell volume lying
more than 1.40 Å from every framework O (Shannon O²⁻ [P7]) and more than the
Si/Al radii from T atoms. It is computed on the realized, scaled cell with a
deterministic estimator of stated tolerance (D2). The exclusion radius is an
agent choice and strongly affects the result (review-2 sensitivity, Na-A
fixture):

| O exclusion radius (Å) | 1.2 | 1.35 | 1.40 (proposed) | 1.6 | 1.8 | 2.0 |
|---|---|---|---|---|---|---|
| V_free,α (Å³) | 893 | 875 | 867 | 816 | 744 | 667 |
| n_melt (x_Li 0.58) | 15.5 | 15.2 | 15.0 | 14.1 | 12.9 | 11.6 |

The 1.35–1.40 Å radii follow the usual O-radius convention behind
free-diameter figures such as the 11.4 Å of [R31]. Radii of 1.4–2.0 Å are
"peak-consistent" (Cl–O peak minus half the Cl–Cl peak). The choice moves the
12/α Li target from 3.4 to 4.3, more than δ.

| α volume definition | V (Å³) | Cl at 800 K melt density |
|---|---|---|
| Nearest-centre free volume (adopted; Na-A, V_ref, agent Monte Carlo) | 867 | **15.0** (14.7 also excluding ring Na) |
| 11.4 Å free-diameter sphere [R31] | 776 | 13.4 |
| 11.05 Å IZA maximum included sphere (idealized LTA) | 707 | 12.2 |
| Adopted definition with f_V = 0.97 applied | about 841 | about 14.6 |

**Finding (T17).** At the proposed radius and for x_Li ≳ 0.07, the 800 K
melt capacity of an α-cage (12.2–15.0 across the volume definitions above) is
at or above the D-2 cap of 12. Such cages are never "dense". K-rich cages
(pure KCl: n_melt ≈ 11.7, or 11.3 with f_V) and radii of 2.0 Å or more can be
dense. The measured capacity of about 10.5 Cl per pseudo-cell [R1,
R2], about 9.5 per α with 1 Cl in β, is about 63 % of melt density under the
adopted definition.

**Coordination rule under S-2/S-4 ([CONFIRM]; three readings).** Notation:
φ(c) = min(1, N_Cl(c)/n_melt(c)) and T_M(c) = max(1, T_M^bulk φ(c)). Basis
(grade D): N(M–Cl) = 4π ρ_Cl ∫₀^{r_min1} g_MCl(r) r² dr, which is linear in
the Cl number density *if g is unchanged*. That assumes the salt is spread
uniformly over the cage.

| Reading | Rule in an α-cage | Consequence |
|---|---|---|
| (i) literal | Anchoring and per-ion range everywhere. The band T^bulk ± δ only if φ ≥ 1 (dense). | No mean target in nearly all cages (T17). Unshaped packing then gives Li–Cl about 2 at 8–12 Cl/α (P5). |
| (ii) two-sided linear (r3) | mean ∈ T_M(c) ± δ | At ≤ 4 Cl/α it forces CN of about 1, i.e. isolated ion pairs (P8 at n = 2, 4). That contradicts the cluster evidence below. |
| **(iii) one-sided linear floor (recommended)** | mean ≥ T_M(c) − δ when φ < 1. Two-sided T^bulk ± δ when φ ≥ 1. Per-ion 1 ≤ N_i ≤ U_M always. | Requires enough coordination at high loading and allows clustered salt (higher local CN) at low loading. U_M caps extremes. |

- Example floors T_M(c) at n_melt = 15:

  | Cl per α | 4 | 8 | 9.5 | 12 |
  |---|---|---|---|---|
  | Li | 1.15 | 2.29 | 2.72 | 3.44 |
  | K | 1.39 | 2.77 | 3.29 | 4.16 |

- A physical reason for low Cl coordination of wall cations, if one is
  needed, is coordination by framework O (for example Sr–O3 × 3 at 6R
  [Q25]). This is not modelled.

**Literature on CN vs density (agent search, 2026-10-17).** No study was
found that reports the coordination of a confined chloride melt as a function
of loading or density. Indirect support for the direction only:
- molten NaCl under compression: CN rises from 3.5 at 0.1 MPa to 4.5 at
  5 GPa. Melt density is unknown at pressure, per the authors (Urakawa et al.,
  *Am. Mineral.* 84, 341 (1999), doi:10.2138/am-1999-0317; snippet);
- chloride vapour-phase ions and clusters: CN 2–3 (Rong et al., *J. Energy
  Storage* (2024), doi:10.1016/j.est.2024.110712; abstract);
- an ionic liquid in carbon nanotubes shows reduced coordination under
  confinement (PubMed 38660932; snippet);
- confined α-cage halides in LTA crystals have 3–4 cation neighbours (E-cn
  §4).
- *Counterevidence to a uniform linear reduction:* small chloride clusters
  and LTA halide units keep CN 2–4 at very low overall density. At low
  loading the salt is therefore likely to cluster, not disperse. This is why
  reading (iii) bounds the mean only from below.

The linear form is an assumption (T16).

**Evaluation.** The rule is evaluated per α-cage c and per species M with ≥ 1
centres in c. T_M(c) uses only the cage's own Cl, but first neighbours may
include Cl in adjacent cages through the windows. This is accepted as an
approximation.

**Bound bulk target values T_M^bulk (Option A adopted, §15).**

| M | Option A (adopted) | Option B (not adopted): base-salt E-cn §1c |
|---|---|---|
| Li | 4.3 ± 0.5 (§1a, inside r_min1) | 4.2 ± δ_B |
| Na | 5.2 ± 0.5 (§1a, inside r_min1) | unavailable → A |
| K | 5.2 ± 0.5 (§1a, inside r_min1) | 7.0 ± δ_B (definition-dependent) |
| Rb | 5.4 ± 0.5 (§1a; r_min1 here is an E ratio estimate) | 5.8 ± δ_B (dilute, XAFS, no cutoff) |
| Cs | 5.5 ± 0.5 (§1a; r_min1 here is an E ratio estimate) | unavailable → A |
| Sr | 6 ± 1 (E working value, dilute; not defined inside r_min1) | same |
| Ba | 7 ± 1 (E working value; 6.4 vs 7.7 conflict open) | same |

Option A is adopted because, for Li, Na and K, it is the only set whose
integration boundary is the same `r_min1` the generator counts within (E-cn
§0). For Rb, Cs, Sr and Ba no supplied value shares that boundary, and the
targets are E-grade transfers. Option B matches the base-salt composition but
uses an unstated boundary. Its source ranges (Li 4.1–4.3, K 6.8–7.2) are run
spreads, not tolerances.

**Per-ion admissible range (adopted, §15).** A mean-only target would admit
a bimodal population, for example Li at CN 8 and CN 1. In every α-cage,
every `cage_salt` cation must therefore have 1 ≤ N_i(Cl) ≤ U_M:

| M | Li | Na | K | Rb | Cs | Sr | Ba |
|---|---|---|---|---|---|---|---|
| U_M | 6 | 6 | 7 | 7 | 6 | 8 | 9 |

U_M is the integer ceiling of the upper "spread across all sources" in E-cn
§1a, or T + 2 for Sr and Ba. The lower bound is the anchoring rule, not a
spread bound, because the density-reduced floor goes down to 1.

- **Distribution form.** The evidence supplies means ± spread, not P(CN)
  histograms (O4b answered by data availability). The per-cage mean and the
  per-ion range are targeted. The spread inside that range arises from randomness,
  satisfying D-9, and is reported.
- **Failure policy (adopted, §15).** Violating the coordination rule or the per-ion range
  after the D2 restart budget is a construction failure (fail closed), as is a
  hard-minimum or anchoring violation.

## 8. Ring cations during the salt build (D-12)

A ring cation (role `framework` or `occlusion_counterion`, O1b adopted)
that conflicts with the cage salt may move only to:

- an empty physical site of its own ring family; or
- another positional variant of its own physical site (O1a, adopted).

If no such site or variant resolves the conflict, the cation stays. Species
counts per ring family and role never change, and no cut-off is relaxed for
ring cations. Each relocation is recorded with its site identifiers.
Relocation happens after cation density smoothing, so a relocated ring cation
is no longer at its smoothed site. That is accepted and reported.

## 9. Volume model (D-11)

Let n̄ = N_Cl / 2 (Cl per pseudo-cell) and V_ref the v7.14 dry reference
volume of the framework population. In the cage-salt model the default volume
factor is

```
f_V(n) = 1 - (1 - f_ref) * min(n / n_ref, 1)
V_target = f_V(n_bar) * V_ref
```

**Bindings (all grade E).**

- f_ref = 0.97 ± 0.03, with hard bounds 0.919 ≤ f_ref < 1 (E-vol §5).
- n_ref = 10.5 Cl per pseudo-cell, the calibration loading of 4A + LiCl–KCl
  at 823 K [R1, R2] (E-occ §6.4).

**Form (S-6, 2026-10-17: "go with your recommendation"; recommendation
revised after review 2 from (c) to (a)).**

```
f_V(n) = 1 - (1 - f_ref) * min(n / n_ref, 1),   n_ref = 10.5 Cl per pseudo-cell
```

- **Correction.** The r3 rationale for (c), n_ref = 4, was wrong in two ways.
  - The K/Br datum K₉(K₄Br)·0.75Ag₄Br₄ holds 3 of its 4 Br per pseudo-cell
    in β-cage Ag₄Br₄ clusters (E-vol §3). The model caps β at 1 Cl, so the
    datum is no loading equivalent of n̄ = 4 mostly-α Cl.
  - [R18] states that the partly loaded phase "continues to decrease but never
    reaches the lattice parameter of the new zeolite phase". Contraction thus
    keeps growing up to full loading, which the flat-from-4 form (c)
    contradicts.
- **Adoption.** Form (a) assigns f_ref = 0.97 to the fully loaded state, as
  E-vol §5b does, and is the minimal monotone form consistent with [R18]. It
  is adopted under the stakeholder's delegation, and the change is reported
  to them.
- **Rejected:** (b), a step at n > 0, and (c).
- **Thermal caveat.** At n̄ ≲ 3.5, 1 − f_V < 1 %. In that regime E-vol §4
  warns that the thermal term (about 0.3 %, sign unknown) is not negligible.
- Confidence is low (T11).

- **Fixed framework population.** f_V is an occlusion-only factor at a fixed
  framework population. Ion exchange (Na → Li) is represented only through the
  framework population in V_ref, never through f_V (tension T14).
- No dependence on composition or cation form is modelled, because there are
  no data for one.

**Override (adopted, §15).** An explicit `--volume-factor` *replaces*
f_V(n̄), so V_target = volume_factor · V_ref, exactly the v7.14 meaning
(D-11 "override").

The legacy model and the dry path keep the v7.14/v7.15 behaviour (default
factor 1.0). The temperature mismatch (in-situ 650–800 K contraction applied
to a dry reference of another temperature) is below the stated uncertainty:
the thermal term is about 0.3 % in V (E-vol §4).

## 10. Legacy sparse model (D-10)

The v7.15 sparse model is reachable only through an explicit legacy selector.
For the same inputs and seed its output must be byte-identical to v7.15.0. The
default occlusion model is the cage-salt model. v7.15 flags keep their v7.15
per-cell meaning.

**Legacy flags (adopted, O7, §15).** The v7.15 flags (`--alpha-cl`, `--beta-cl`,
`--occlusion-n4/6/8`, `--occlusion-min-distance`) are accepted only with the
legacy selector, and are rejected otherwise with a message naming the
selector. Old commands then fail loudly instead of being reinterpreted. New
per-cage inputs use new names (D3).

## 11. Construction order and determinism (D2/D3 constraints)

1. framework ring placement (legacy stream);
2. `occlusion_counterion` ring placement (v7.15 secondary stream);
3. cation density smoothing (ring cations only; Cl and `cage_salt` excluded).
   Its greedy repair moves are finished before step 4;
4. charge rebalancing (§4): re-site `occlusion_counterion` cations within
   their ring family to satisfy |q_c| ≤ q_max, on its own stream derived from
   the seed. Seed-level failure is reported here;
5. volume model (§9): V_target depends only on the inputs, so the cell is
   scaled to V_target *before* the salt is built;
6. cage-salt build in the final cell, including D-12 relocation, on its own
   stream derived from the seed. A relocation is accepted only if it keeps
   every |q_c| ≤ q_max;
7. exact verification on the final scaled coordinates. It covers:
   - inventory (§4);
   - cage charge, |q_c| ≤ q_max (§4);
   - hard minima (§6), with no tolerance below r_cut;
   - anchoring, the CN rule and the per-ion range (§7).

   Verification of an unscaled cell does not count, because isotropic scaling
   by f_V^{1/3} would shrink verified distances by up to 1 %;
8. write.

Enabling the cage-salt model must not perturb the streams of steps 1-3. The D2
method (proposed: insertion plus penalty relaxation plus light reverse-Monte-
Carlo CN shaping) is delegated and claims no uniform sampling. Some modes
remain rejected when any occlusion is active (v7.15):
- maximum-symmetry mode;
- unique batches larger than one;
- ratio allocation;
- the v7.15 `--auto-rebalance` of framework allocations (distinct from the
  §4 charge rebalancing);
- candidate-variant output.

## 12. Feasibility evidence (agent probes; exploratory, not acceptance evidence)

Scripts: `workplans/v716-probes/`. P5 and P6 used the E-pair r_cut and
r_min1 columns exactly, the D-12 relocation, the fixture
`tests/real_lta_primitive.cif`, a base salt at x_Li = 0.592, and 1 Cl + 1 Li
in each β (B1a) unless stated. The framework proxy F2 is option (b) of §6.
Si/Al used 0.85 × Shannon sums. Except in the sensitivity grid, the cage region
was C1: a 5.7 Å sphere about each α centre, plus a 1.0 Å confinement of β
guests that neither O8 option contains. Builds took place in the scaled cell.
"Feasible" means the exact final check passed with slack ≥ −1e-3 Å, a probe
tolerance that §11 does not allow for production. Reported CN is α-only.

| Probe | Variants | Result |
|---|---|---|
| P5 main | ring Na/K/Li-A × seeds 0-2 × n = 8, 10, 11, 12 per α × relocation on/off; 8 starts | Relocation on: 24/24 at every point. Fixed rings: 19-23/24. |
| P5 sensitivity (Na-A, seed 0) | F1/F2 × C1/C2 × β B0/B1a/B1b × n 10, 12 × f_V 1.0/0.97, plus O1a off | 8/8 at all 50 points |
| P5 dopants | {Li5 K4 Cs1}, {Li5 K3 Sr1}, {Li6 K4 Ba1} × 3 seeds | 9/9 points 8/8 |
| P5 Ca ring | Ca-A | rejected: missing pair (Ca, Cl), as §5 requires |
| P5 unshaped CN | all feasible points | Group means Li–Cl 1.9-2.1 and K–Cl 3.5-3.7. Per-point ranges Li 1.2-2.8, K 2.5-5.1. Sr–Cl 2.0-2.5, Ba–Cl 3.0-4.5. Min 0 (isolated ions). Ring–ring slack ≥ +0.27 Å |
| P7 density-reduced CN (Option A, superseded reference N_Cl/12) | Na/K-A × seeds 0-2 × n 8, 10, 12 × nearest-centre partition (C2); 4 starts | 18/18 points: per-α-pooled means within ±0.5 of T·n/12, 1 ≤ N_i ≤ U_M for every Li/K, hard minima met. Achieved Li 2.6-4.4, K 3.0-5.5 against targets 2.87-4.30 and 3.47-5.20 |
| P8 melt-density CN (S-4 reading, Option A, n_melt = 15) | Na/K-A × seeds 0-2 × n 2, 4, 6, 8, 12 × C2; 4 starts; targets max(1, T·n/15) | Hard minima met at all 30 points. Pooled α means lie within ±0.5, counting band edges as inside. Under spec rules 26/30 comply: 4 points (n = 2 three times, n = 6 once) leave one cation with CN 0, because the probe's smooth objective does not enforce anchoring, so D2 must enforce it explicitly. The probe pools both α-cages and uses n_melt = 15 throughout, whereas §7 is per cage (x_Li 0.5–0.67 would give 14.5–15.9). All points also satisfy reading (iii). |
| P6 CN-shaped | Na/K-A × seeds 0-2 × n 8, 10, 12 × options A/B; 4 starts | All hard minima met. α-only, Option A: Li 3.9-4.5, K 5.1-5.5, all 18 points within ±0.5. Option B: Li 3.6-4.3, K 6.2-7.1; 3 of 18 points outside ±0.5, and 17 of 18 outside the source ranges. Cl-centred (cage + ring) 5.4-7.1. Anchoring min ≥ 1 (α). If β-Li (CN 1) were included in the mean, 6 of 18 Option A points would miss (review 1). |

**Variant disclosure (cumulative v7.16).**
- P1-P4 (previous session): about 90 points (workplan §5).
- P5: 132 points.
- P6: 36 points.
- P7: 18 points (one form, chosen after the stakeholder's S-2 instruction,
  not selected from alternatives).
- P8: 30 points (the S-4 reading, one form).
- The α free-volume estimate: one Monte Carlo run (4×10⁵ points, Na-A,
  seed 0, O radius 1.35 and 1.40 Å), with the two sphere definitions as
  comparison.
- Every evaluated point is reported. Nothing was selected as a default from
  these results. All probe parameters were chosen by the agent, as were the F1
  proxy, the 5.7 Å sphere, the smooth CN surrogate width of 0.12 Å and the
  restart counts.
- Review 1 re-ran 2 probe cases with parameters unchanged and reproduced the
  stored values. It found 0 nearest-centre misassignments in a C1 build at
  12/α, and a minimum ion–framework distance of 1.72 Å.
- Not examined:
  - mixed or Group-II ring populations;
  - more than 3 seeds;
  - dopant-rich salts;
  - n < 2;
  - the §4 charge rebalancing and the q_c criterion;
  - CN shaping with Cs/Sr/Ba;
  - the per-ion range rule, which was not a shaping objective. Recounted after
    the fact, 2 of 18 Option A P6 points have at least one α ion outside it;
  - β anchoring;
  - option (c) of §6;
  - the spec's own cage partition beyond 24 points;
  - realism or uniformity of the configurations;
  - CN-band misses under a production restart budget.

## 13. Tension assessment (agent-asserted, proposed until ratification)

Searched: `doc/specs`, `doc/archs`, `workplans/` (including archive), the
`qualification/` records, the repository's GitHub issues (none exist) and the
evidence package. Recorded items found: the v7.15 archived reopen triggers
(three fire, see Predecessor) and the v7.14 §9 limitations (consistent). No
recorded tension against v7.13 was found.

| ID | Tension | Assessment for this revision |
|---|---|---|
| T1 | Crowding at 10-12 Cl/α under rock-salt cut-offs | **Inapplicable to this revision.** The rock-salt bases are replaced by the E-pair onset r_cut, and P5 is feasible at 12/α. Reopen if r_cut is recalibrated upward (E-pair §6.1, Zeidler NaCl partials). |
| T2 | Unsourced melt densities | **Superseded by T19.** The E-occ §5 handbook densities are not used. The §7 density reference uses the cited Janz correlations. |
| T3 | No ion–framework or dopant pair data | **Resolved by S-1 with an approximation.** Rb, Sr and Ba pairs are supplied, mostly grade E. Ion–framework minima use the free-oxide rows (grade E). Reopen if framework-O data arrive. |
| T4 | R10/R14 author lists vs v7.15 citations | **Resolved.** E-occ rev 2.1 lists R14 as Simpson & Gougar and R10 as Lewis, Fischer, Smith, matching v7.15. |
| T5 | Trivalent dopants outside v7.13 | **Applicable; out of scope** (§5). v7.13 is not reopened. |
| T6 | K–Cl CN is definition-dependent (5.2 vs 7.0) | **Applicable; Option A adopted** (consistent with the r_min1 boundary). |
| T7 | Cl-centred CN with ring neighbours (6-7) exceeds the bulk eutectic value (5.1-5.3) and confined-crystal values (3-4) | **Applicable.** Reported, not targeted. |
| T8 | Unshaped packing isolates ions | **Applicable** → anchoring rule (§7). D2 must shape CN. |
| T9 | BaCl₂ CN conflict (6.4 vs 7.7) | **Applicable.** E working value 7 ± 1. |
| T10 | Ca/Mg forms (5A) excluded | **Applicable** → §5 decision. |
| T11 | Contraction rests on one E-grade analogue point (K/Br-LTA at 4 Br per pseudo-cell, 294 K), with 3 of the 4 Br in β-cage Ag₄Br₄, so it is no loading equivalent for this model. The loading dependence is assumed: form (a) is linear to 10.5. The r3 form (c) was withdrawn after review 2 | **Applicable.** Reopen when matched 4A data or the user's DFT/MLFF f_V(n) arrive (E-vol §5b route ii). |
| T13 | Bulk melt CN targets are undefined at low loading (review 1 SC1) | **Resolved by S-2 and S-4:** density-reduced target max(1, T·φ). Interpretation per §7 [CONFIRM]. |
| T14 | The f_V bound (f_a ≥ 0.972) equals the v7.14 Li-A/Na-A ratio 0.9724, and 4A exchanges Na → Li in LiCl–KCl [R22, R32]. The R18 small-a phase may be partly an exchange effect, so f_V on a Li-rich V_ref could double-count (hypothesis, review 1) | **Applicable, open.** f_V is defined as occlusion-only (§9). Reopen with matched data. |
| T15 | Charge separation between α-cages is permitted without a per-cage rule (review 1) | **Resolved by S-5:** criterion abs(q_c) ≤ q_max with rebalancing (§4). Threshold per T18. |
| T16 | No study gives CN vs loading or density for a confined chloride melt. The linear-in-density form (§7) assumes uniform dispersal. Cluster evidence (CN 2–4 at low overall density) is counterevidence to applying it two-sided | **Applicable.** Reopen with an MD-computed CN(n) in LTA, for example from the user's MLFF (cf. E-cn §6.3). |
| T17 | The 800 K melt-density capacity of an α-cage (11.6–15.5 Cl across volume definitions and radii, 15.0 proposed) is mostly at or above the D-2 cap of 12. Only K-rich cages or large exclusion radii can be dense | **Applicable** → §7 [CONFIRM] readings and the proposed radius. |
| T18 | The q_max threshold is not calibrated, and no study quantifies cage-scale charge fluctuations. A bulk-sphere statistic is the wrong estimand (review 2, SC-A) | **Applicable** → stakeholder-set provisional 1 e; calibrate on salt-loaded LTA MD (§4). |
| T20 | The r3 q_c excluded framework cations, hiding a ±1 e local imbalance (review 2) | **Resolved:** q_c^fw added (§4). |
| T19 | Melt density rests on Janz correlations read via citing papers (grade C), ideal mixing, and no dopant volumes | **Applicable, minor.** The check reproduces the cited eutectic value to 10⁻⁴ g cm⁻³. |
| T12 | Committed evidence file was rev 2.0. Rev 2.1 corrects the packing range (79-87 %, not 70-80 %) and R33/R32 conversions | **Resolved by committing rev 2.1.** The workplan T1 wording is superseded. |

No SERIOUS CHALLENGE against v7.13, v7.14 or the retained parts of v7.15 is
raised. v7.15's sparse semantics are superseded only by explicit model
selection, not found false.

## 14. Inspectability (O1)

**Material realized record.**
- Per-cage inventory by role and species.
- Ring relocations (ion, from-site, to-site).
- Per pair class: achieved minimum distance and slack against r_cut.
- The catastrophic-guard minimum.
- Ring–ring slack.
- Per-species CN mean, histogram and band, plus Cl-centred CN.
- Anchoring check.
- Restarts and failure reasons.
- Seed and derived stream identities.
- Bound evidence-file SHA-256s.
- Volume inputs: n̄, f_V, V_ref, and whether the factor was overridden.
- Per-cage excess charge q_c, the rebalancing moves, and q_max with its
  provenance (provisional or calibrated).
- n_melt(c), V_free,α(c), φ(c) and T_M(c) per α-cage.
- α and β CN reported separately, with the CN option, band and per-ion range
  in force.
- The x_M weights and the resulting Cl–Cl r_cut.
- Ion–framework slack per class.
- The minimum cage-assignment margin.
- A statement that verification ran on the final scaled cell.
- Ring-cation density after relocation (nearest-neighbour statistics).

**Reader.** The stakeholder preparing MD runs.

**Routine questions.**
- Did my requested loading go where I asked?
- How tight was the packing?
- Were ring cations moved?
- How far is the local coordination from the melt targets?
- Which volume was applied, and why?

**Falsifiers and reopen triggers.**
- MD shows systematic expulsion or crystallization traceable to the initial
  minima or CN.
- Recalibrated r_cut values.
- Matched f_V data.
- A CN definition resolved against the chosen option.

**Product surface (marked).** A per-run machine-readable report holding the
record above, written beside each structure. It goes beyond the stakeholder's
original request. It was **adopted as an agent recommendation under the
stakeholder's 2026-10-17 delegation (§15)**, and it is confirmed by
ratification of this document.

## 15. Decision record and open decisions

**Stakeholder decisions (human, via chat, 2026-10-17).**
- **S-1.** Ion–framework minima use the free-oxide rows as an approximation
  (§6).
- **S-2.** Coordination applies only in the dense case, and the expected CN
  is reduced with density; the exact relation is unknown and the agent was
  asked to check the literature (§7).
- **S-3.** "If you can handle the rest of the questions with certainty, go
  with your recommendation."

**Adopted under S-3** (agent recommendations; asserter: AI agent, this
session; confirmed only by ratification of the whole document):
- Option A bulk targets.
- The per-ion range 1 ≤ N_i ≤ U_M.
- Fail-closed on a CN miss.
- Exclusion of species without pair rows.
- Cell-wide Cl–Cl weights.
- `--volume-factor` replaces f_V.
- O1a and O1b: yes.
- O8: nearest-centre partition.
- β salt cations ≤ N_Cl(β); cation-only cages rejected.
- O7: legacy flags only with the legacy selector.
- The per-run report as a product surface.
- O6: trivalent cations out of scope.

**Stakeholder decisions, second round (human, via chat, 2026-10-17).**
- **S-4.** "Dense" means the density of the molten salt at 800 K (§7).
- **S-5.** Charge separation should not be too severe. The code computes the
  charge distribution and rebalances concentrated regions; the threshold could
  be inspired by molten-salt charge-fluctuation studies (§4).
- **S-6.** Contraction form: "go with your recommendation". The
  recommendation was revised after review 2 from (c) to (a) (§9). The change
  is reported to the stakeholder.

**Superseded.** The provisional dense threshold n_CN = 8, the reference of 12
(by S-4), and the bulk-sphere q_max calibration (by review 2).

**Agent confidence in items adopted under S-3.** S-3 was conditional on
certainty.
- *High:* the exclusion of species without pair rows (D-3); legacy-flag
  handling; replace-semantics of `--volume-factor` (D-11 "override"); O1a/O1b
  (the D-12 move set); β salt cations ≤ N_Cl(β); trivalent out of scope.
- *Moderate:*
  - O8 nearest-centre (only 24 points tested; the sphere alternative is
    incomplete);
  - E-grade targets for Rb/Cs/Sr/Ba;
  - U_M values;
  - fail-closed on a CN miss;
  - cell-wide Cl–Cl weights.

  The stakeholder may reverse any of these at ratification.

**Open [CONFIRM] at ratification.**
1. §7: coordination reading (i), (ii) or (iii). Reading (iii), the one-sided
   floor, is recommended and provisionally written in.
2. §7: V_free,α definition and O exclusion radius 1.40 Å (proposed), and the
   x_Li fallback and dopant-validity flag.
3. §4: the q_c observable including q_c^fw (proposed); q_max = 1 e
   (provisional placeholder, stakeholder-set); calibration from salt-loaded
   LTA MD (optional).
4. §4: rebalancing never changes per-cage salt inventory (D-1).

## 16. D1 → D2 hand-off (on ratification)

D2 must preserve:
- §4 exact inventory and invariants;
- §6 hard minima with exact final verification (no tolerance below r_cut);
- §7 anchoring (enforced explicitly; see P8), the density-reduced CN band per
  α-cage, the per-ion range, and V_free,α computed on the realized cell;
- §4 the cage charge criterion (q_c^exc + q_c^fw), the rebalancing move set,
  relocation acceptance under the criterion, and seed-level failure
  reporting;
- a deterministic V_free,α estimator with a stated tolerance, so output stays
  byte-deterministic;
- §8 relocation limits;
- §11 stream isolation;
- fail-closed behaviour.

D2 owns:
- the packing and relaxation method;
- CN shaping;
- the restart budget, which must be calibrated so feasible requests rarely
  fail;
- tie handling;
- how the CN miss probability is bounded.

D2 must not:
- add a radius or scaling fallback;
- reduce cut-offs against ring cations;
- lower the 12/α cap;
- claim uniform or equilibrium sampling.

## References

Evidence tags [Rn], [Pn], [Qn] resolve in the `doc/evidence/` files listed in
§1. The primary sources most used here:
- [R1] D. Lexa, I. Johnson, *Metall. Mater. Trans. B* 32, 429 (2001).
- [R2] D. Lexa, *Metall. Mater. Trans. B* 34, 201 (2003).
- [R18] D. Lexa, L. Leibowitz, A. J. Kropf, APS Activity Report 2001.
- [R31] N. C. Hyatt et al., *MRS Proc.* 807 (2004).
- [P1] M.-M. Walz, D. van der Spoel, *PCCP* 21, 18516 (2019).
- [P2] S. Roy et al., *PCCP* 22, 22900 (2020).
- [P22] A.-L. Rollet, M. Salanne, *Annu. Rep. Prog. Chem. C* 107, 88 (2011).
- [P23] D. K. Belashchenko, O. I. Ostrovski, SAIMM (2004).
- [Q14] J. Fuller et al., *PCCP* 24, 13102 (2022).
- [Q1] S. Turner et al., *Microporous Mesoporous Mater.* 107, 252 (2008).

All v7.13, v7.14 and v7.15 references remain applicable within their scopes.

## 17. Review record

**Independent falsification review 1** (2026-10-17). AI reviewer: a fresh
Claude Code subagent context of this session, with the same model family,
source corpus and toolchain as the author. It is therefore process-independent
but has common-mode risk.

**SERIOUS CHALLENGES.**
- **SC1 (low loading).** Resolved by S-2/S-4 and the §7 readings.
- **SC2 (β cations in the CN mean).** Fixed: the target is now α-only and the
  P6 claim is corrected.

**Findings and dispositions.**
- **Blocking: verify-then-scale.** Fixed in §11.
- **Major:**
  - per-cage charge separation → §4 decision;
  - option (c) inadmissible → §6;
  - Option A rationale overstated for Rb/Cs/Sr/Ba → §7;
  - Option B "±0.5" claim false → §7, §12;
  - mean-only CN → per-ion range decision;
  - volume datum at 4 Br contradicts the linear form → §9 decision, T11;
  - exchange double-count hypothesis → T14;
  - "active" undefined → §3;
  - probe and spec mismatches → §12.
- **Minor:**
  - authority boundary notes → §1;
  - partition geometry → §3;
  - β cap ≤ N_Cl(β) → §4;
  - §2 wording → fixed;
  - relocation diagnostics → §14;
  - record omissions → §14;
  - unapplied evidence README corrections → below.

The review verified as correct:
- the CSV hashes and byte-identical regeneration;
- the spot-checked r_cut and r_min1 values;
- the CN table values;
- f_V 0.97 ± 0.03, its bounds and n_ref;
- the reopen triggers;
- the charge-invariant generalisation;
- T4;
- every §12 count.

**Evidence corrections not yet applied.** The evidence README lists
corrections to the rev 2.1 pair file (drop P14; P5 SI provenance;
bibliographic completions). None changes a bound value. The CSVs are bound by
hash as supplied.

**Independent falsification review 2** (2026-10-17) covered the r2/r3 delta.
AI reviewer: a fresh subagent context of this session, with the same
common-mode caveat as review 1.

**SERIOUS CHALLENGE.**
- **SC-A.** The bulk-sphere calibration of q_max measures the wrong quantity
  and would allow about 5 e. Fixed: q_max is stakeholder-set, and
  calibration, if any, uses salt-loaded LTA and q_c itself.

**Findings and dispositions.**
- **Major:**
  - framework cations excluded from q_c → q_c^fw added;
  - rebalancing missing from the step order, unverified, ill-defined
    reachability → §11 steps 4, 6 and 7, and the §4 failure split;
  - V_free sensitivity, composition and T17 overstatement → §7 table, x_Li
    fallback, T17 qualified, marked proposed;
  - stale r3 text → fixed;
  - the [CONFIRM] framing (dropped "only", "adopted" label, a single reading,
    self-confirmation) → three readings, S-2 quoted verbatim, labels
    corrected;
  - S-6 basis wrong → recommendation revised to (a);
  - the linear form's physics → counterevidence recorded and reading (iii)
    recommended.
- **Minor:**
  - P8 row precision → fixed;
  - correlation extrapolation and circularity → T19 text;
  - literature wording (zeroth vs second moment; the unverifiable J. Phys.:
    Condens. Matter snippet dropped; E-occ not cited as support) → fixed;
  - deterministic V_free estimator and cross-window neighbours → §16, §7.

Review 2 verified as correct:
- the Janz ideal-mixing densities;
- 867 Å³ and n_melt 15.0;
- the example targets;
- the window topology;
- Σ q_c = 0;
- the P7/P8 counts;
- the f_V arithmetic;
- the E-vol rows;
- the CSV hashes.

A third independent review is required before ratification if the stakeholder's
[CONFIRM] answers change content beyond the listed options.
