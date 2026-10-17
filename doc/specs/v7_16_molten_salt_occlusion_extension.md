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
One interpretation item remains, marked **[CONFIRM]**, and one threshold is
marked **[CALIBRATE]**. The document stays PROPOSED until the stakeholder
ratifies it as a whole.
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
threshold is to be informed by studies of charge fluctuations in molten salts.

- **Observable.** The excess charge of cage c:

  ```
  q_c = sum_{i in cage_salt(c)} z_i - N_Cl(c) + sum_{j in occlusion_counterion} z_j w(j, c)
  ```

  The sharing weights w follow the window topology, checked on the canonical
  template:
  - a 4R site belongs to one α: w = 1 (5.0 Å from that α centre, 7.2 Å from
    the nearest β);
  - a 6R site separates one β and one α: w = ½ each (4.24 Å and 6.4 Å);
  - an 8R site separates two α: w = ½ each (6.27 Å from both).

  `framework` cations are excluded, because they compensate the framework
  charge (+24 invariant, v7.13). Their local distribution remains a
  v7.13–v7.15 matter.
- **Criterion.** For every cage, |q_c| ≤ q_max.
- **Rebalancing (proposed mechanism).**
  - After ring placement and before the salt build, `occlusion_counterion`
    cations may be re-sited to empty sites of their own ring family. This is
    the D-12 move set and is needed to bring every |q_c| ≤ q_max.
  - Per-cage salt inventories are input (D-1) and are never changed.
  - If the inventory alone makes the criterion unreachable, the input is
    rejected up front, naming the cages. "Unreachable" means that even the
    best assignment of counterions to available sites fails.
  - Changing per-cage salt inventories would require relaxing D-1.
- **Threshold (provisional q_max = 1 e; [CALIBRATE]).**
  - **Literature (agent search, 2026-10-17).** In a molten salt, perfect
    screening (the Stillinger–Lovett conditions) forces the charge–charge
    structure factor S_ZZ(k) → 0 as k → 0, as shown for molten NaCl
    [P4: Zeidler et al., *J. Chem. Phys.* 157, 094504 (2022)]. The net charge
    in a sub-volume therefore fluctuates far less than the ion count. For large
    volumes its variance grows with the boundary area rather than the volume
    (charge hyperuniformity):
    - Ph. A. Martin, T. Yalcin, *J. Stat. Phys.* 22, 435 (1980);
    - J. L. Lebowitz, *Phys. Rev. A* 27, 1491 (1983);
    - confirmed by simulation in T. Hoang Ngoc Minh, B. Rotenberg,
      S. Marbach, *Faraday Discuss.* 246, 225 (2023).

    MD of molten NaCl puts S_ZZ near Q ≈ 9 nm⁻¹ at about 10 % of S_NN (cited
    in *J. Phys.: Condens. Matter*, doi:10.1088/1361-648X/ac101c). These were
    read via search abstracts and snippets only; this environment's network
    policy blocks the full texts.
  - **Gap.** No study found gives the distribution of net charge in a
    cage-sized volume (about 870 Å³) of LiCl–KCl at 800 K. A cage holds about
    1.5–2 ion spacings, too small for the asymptotic area law.
  - **Calibration route (bound for D2 or evidence).**
    - Run bulk LiCl–KCl eutectic MD at 800 K, with the stakeholder's MLFF or
      a published polarizable model.
    - Sample the net formal charge Q in probe regions of volume
      V_free,α ≈ 870 Å³ (equal-volume sphere radius ≈ 5.9 Å).
    - Set q_max = ⌈3 σ_Q⌉.
  - **Provisional value.** Until that run, q_max = 1 e. This is an
    agent-chosen conservative placeholder, not sourced. It matches "each
    occluded Cl⁻ arrives with a charge-balancing cation" (E-occ §0, §6.3).
    The ½ weights let a single β counterion on a 6R site satisfy it.
- **Reported:** every q_c, the rebalancing moves, and the distribution of q_c
  against q_max.

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
must lie in the band T_M ± δ_M. β-cage cations are excluded. Under D-2 they
have at most one Cl in their own cage (CN of about 1, as measured in review 1).
Their CN is reported separately and subject only to anchoring. Ring cations
count as neighbours of Cl (O4c, proposed yes), but are not CN centres. Only
cation-centred targets are imposed (D-13). Cl-centred counts are reported, not
targeted, because:
- by number conservation they follow from the cation targets plus the ring
  contribution;
- with ring cations counted they exceed bulk values (finding T7).

**Density reference (stakeholder decisions S-2 and S-4, 2026-10-17).**
- S-2: the bulk coordination applies in the dense case, and expected
  coordination is reduced with density.
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
  pair (57.4 Å³ at x_Li = 0.592). For a cage, the Li:K ratio of its
  `cage_salt` cations is used.
- **Dopants.** Their volume is not modelled (dilute assumption), because no
  sourced density is bound for them.

**Cage capacity at that density.** n_melt(c) = V_free,α(c) / v_pair(800 K,
x_Li(c)). The adopted V_free,α is the cage's nearest-centre cell volume lying
more than 1.40 Å from every framework O (Shannon O²⁻ [P7]) and more than the
Si/Al radii from T atoms. It is computed on the realized, scaled cell.

| α volume definition | V (Å³) | Cl at 800 K melt density |
|---|---|---|
| Nearest-centre free volume (adopted; Na-A, V_ref, agent Monte Carlo) | 867 | **15.0** (14.7 also excluding ring Na) |
| 11.4 Å free-diameter sphere [R31] | 776 | 13.4 |
| 11.05 Å IZA maximum included sphere (idealized LTA) | 707 | 12.2 |
| Adopted definition with f_V = 0.97 applied | about 841 | about 14.6 |

**Finding (T17).** Under every definition the 800 K melt capacity (12.2–15.0)
is at or above the D-2 cap of 12, so a "dense" α-cage cannot occur in this
model's domain. The measured capacity of about 10.5 Cl per pseudo-cell [R1,
R2], about 9.5 per α with 1 Cl in β, is about 63 % of melt density under the
adopted definition.

**Density-reduced target (adopted reading, [CONFIRM]).** Read literally, S-2
and S-4 together would mean no α-cage ever gets a coordination target, since
the cap stays below melt density. The adopted reading, which needs the
stakeholder's confirmation, is instead:

```
phi(c)  = min(1, N_Cl(c) / n_melt(c))
T_M(c)  = max(1, T_M^bulk * phi(c))      for every alpha-cage with N_Cl(c) >= 1
```

- T_M^bulk applies only at or above melt density.
- The floor of 1 is the anchoring rule.
- Basis (grade D): by definition, N(M–Cl) = 4π ρ_Cl ∫₀^{r_min1}
  g_MCl(r) r² dr, so the coordination number scales linearly with the Cl
  number density if g is unchanged.
- Example targets at n_melt = 15:

  | Cl per α | 4 | 8 | 9.5 | 12 |
  |---|---|---|---|---|
  | Li | 1.15 | 2.29 | 2.72 | 3.44 |
  | K | 1.39 | 2.77 | 3.29 | 4.16 |

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
The linear form is an assumption (T16).

**Evaluation.** The band T_M(c) ± δ_M is evaluated per α-cage c and per
species M with ≥ 1 centres in c.

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
a bimodal population, for example Li at CN 8 and CN 1. In dense α-cages
every `cage_salt` cation must therefore have 1 ≤ N_i(Cl) ≤ U_M:

| M | Li | Na | K | Rb | Cs | Sr | Ba |
|---|---|---|---|---|---|---|---|
| U_M | 6 | 6 | 7 | 7 | 6 | 8 | 9 |

U_M is the integer ceiling of the upper "spread across all sources" in E-cn
§1a, or T + 2 for Sr and Ba. The lower bound is the anchoring rule, not a
spread bound, because the density-reduced target goes down to about 2.9.

- **Distribution form.** The evidence supplies means ± spread, not P(CN)
  histograms (O4b answered by data availability). The per-cage mean and the
  per-ion range are targeted. The spread inside that range arises from randomness,
  satisfying D-9, and is reported.
- **Failure policy (adopted, §15).** Missing the band or the per-ion range
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
- n_ref = 4 Cl per pseudo-cell (S-6, below).

**Form (stakeholder decision S-6, 2026-10-17: agent recommendation (c)
adopted).**

```
f_V(n) = 1 - (1 - f_ref) * min(n / n_ref, 1),   n_ref = 4 Cl per pseudo-cell
```

This is the only considered form that agrees with the one quantitative datum
(the K/Br analogue: 4 Br per pseudo-cell, f_V = 0.968; E-vol rows 21/7 and CSV
row 28 note). It is consistent with [R18], which reports contraction already
in partially loaded material. The rejected forms were (a), linear to
n_ref = 10.5, which excludes the datum at its own loading, and (b), a step at
n > 0. Confidence is low, because the datum is a poor proxy (K/Ag/Br, solution
route) (T11).

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
3. cation density smoothing (ring cations only; Cl and `cage_salt` excluded);
4. volume model (§9): V_target depends only on the inputs, so the cell is
   scaled to V_target *before* the salt is built;
5. cage-salt build in the final cell, including D-12 relocation, on its own
   stream derived from the seed;
6. exact verification on the final scaled coordinates: inventory (§4), hard
   minima (§6), anchoring and CN (§7). No tolerance below r_cut. Verification
   of an unscaled cell does not count, because isotropic scaling by f_V^{1/3}
   would shrink verified distances by up to 1 %;
7. write.

Enabling the cage-salt model must not perturb the streams of steps 1-3. The D2
method (proposed: insertion plus penalty relaxation plus light reverse-Monte-
Carlo CN shaping) is delegated and claims no uniform sampling. Maximum-symmetry
mode, unique batches larger than one, ratio allocation, automatic rebalancing
and candidate-variant output remain rejected when any occlusion is active
(v7.15).

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
| P7 density-reduced CN (S-2 form, Option A) | Na/K-A × seeds 0-2 × n 8, 10, 12 × nearest-centre partition (C2); 4 starts | 18/18 points: per-α-pooled means within ±0.5 of T·n/12, 1 ≤ N_i ≤ U_M for every Li/K, hard minima met. Achieved Li 2.6-4.4, K 3.0-5.5 against targets 2.87-4.30 and 3.47-5.20 |
| P8 melt-density CN (S-4 reading, Option A, n_melt = 15) | Na/K-A × seeds 0-2 × n 2, 4, 6, 8, 12 × C2; 4 starts; targets max(1, T·n/15) | 30/30: means within ±0.5 and hard minima met. 4/30 (n = 2 and 6) leave one cation with CN 0, because the probe's smooth objective does not enforce anchoring. D2 must enforce anchoring explicitly. |
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
  - n < 8 (the SC1 threshold);
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
| T11 | Contraction rests on one E-grade analogue point (K/Br-LTA at 4 Br per pseudo-cell, 294 K), with the loading dependence assumed (S-6 form (c) reaches full contraction at that point) | **Applicable.** Reopen when matched 4A data or the user's DFT/MLFF f_V(n) arrive (E-vol §5b route ii). |
| T13 | Bulk melt CN targets are undefined at low loading (review 1 SC1) | **Resolved by S-2 and S-4:** density-reduced target max(1, T·φ). Interpretation per §7 [CONFIRM]. |
| T14 | The f_V bound (f_a ≥ 0.972) equals the v7.14 Li-A/Na-A ratio 0.9724, and 4A exchanges Na → Li in LiCl–KCl [R22, R32]. The R18 small-a phase may be partly an exchange effect, so f_V on a Li-rich V_ref could double-count (hypothesis, review 1) | **Applicable, open.** f_V is defined as occlusion-only (§9). Reopen with matched data. |
| T15 | Charge separation between α-cages is permitted without a per-cage rule (review 1) | **Resolved by S-5:** criterion abs(q_c) ≤ q_max with rebalancing (§4). Threshold per T18. |
| T16 | No study gives CN vs loading or density for a confined chloride melt. The linear-in-density form (§7) is an assumption; indirect literature supports only the direction | **Applicable.** Reopen with an MD-computed CN(n) in LTA, for example from the user's MLFF (cf. E-cn §6.3). |
| T17 | The 800 K melt-density capacity of an α-cage (12.2–15.0 Cl, by definition) is at or above the D-2 cap of 12, so the "dense" case never occurs | **Applicable** → the §7 [CONFIRM] reading. |
| T18 | The q_max threshold is not calibrated; no study quantifies cage-scale charge fluctuations in LiCl–KCl | **Applicable** → §4 calibration route; provisional 1 e. |
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
- **S-5.** Charge separation must not be too severe. The code computes the
  charge distribution and rebalances concentrated regions, with a threshold
  informed by molten-salt charge-fluctuation studies (§4).
- **S-6.** Contraction form: the agent recommendation (c) is adopted (§9).

**Superseded by S-4.** The provisional dense threshold n_CN = 8 and the
reference of 12.

**Open.**
1. **[CONFIRM]** §7 density-reduced reading: the target applies at all
   loadings as max(1, T·φ), because the cap stays below melt density (T17).
2. **[CALIBRATE]** §4 q_max (provisional 1 e) from bulk-melt MD at 800 K.
   This can be done after ratification, as a D2 or evidence task.
3. Implicit in S-5 and confirmed by the reading in §4: rebalancing re-sites
   only occlusion counterions and never changes per-cage salt inventory
   (D-1).

## 16. D1 → D2 hand-off (on ratification)

D2 must preserve:
- §4 exact inventory and invariants;
- §6 hard minima with exact final verification (no tolerance below r_cut);
- §7 anchoring (enforced explicitly; see P8), the density-reduced CN band per
  α-cage, the per-ion range, and V_free,α computed on the realized cell;
- §4 the cage excess-charge criterion and rebalancing move set;
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
- **SC1 (low loading).** Carried as decision 15.2.
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

A second review is required only if ratification changes the scientific
content beyond choosing among the options listed in §15.
