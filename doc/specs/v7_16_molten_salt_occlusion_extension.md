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

**Lifecycle: ACCEPTED-CURRENT D1 (human-ratified 2026-10-17; r8 amendment
2026-10-19).** The stakeholder ratified revision r7 in chat on 2026-10-17
("Sounds good. I ratify."). Revision r8 adds stakeholder decision S-13, a
salt-ion–T minimum (§6). It was chosen by the stakeholder on 2026-10-19 after
finding F-E1-1 (§13, T23). The authoring session (an AI agent) records the
ratification and the decision; it does not supply them.
- **r8 status.** The S-13 value is stakeholder-set. Keeping V_eff O-defined
  (§7) was an agent proposal; the stakeholder acknowledged it on 2026-10-19
  ("Keep V_eff O-defined, I acknowledge.", S-14, §15).
- **Census.** The S-12 census recorded under r7 is void for r8, because the
  frozen code changed. Rule B was re-qualified by census run 3 on the r8 code
  (90/90; `qualification/v716/census-record.md`).

Drafted under workplan `LTA-V7.16-MOLTEN-SALT-OCCLUSION`,
SSDP 7.2.1, by Claude Code session
https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi.

- **What it governs.** The scientific meaning of the v7.16 cage-salt model,
  for the D2 method, D3 design and D4 implementation that follow.
- **What it does not do.** It does not change shipped software behaviour.
  Until v7.16 is implemented and qualified, the released generator still
  behaves as v7.15.0. It does not self-qualify any implementation either:
  independent post-implementation review remains required.
- **History.** It was revised r1-r7 through three independent falsification
  reviews (§17), stakeholder decisions S-1 to S-12 (§15), and r8 for S-13.
- **Proposed items.** Items marked **(proposed)** were agent-chosen defaults.
  They are binding as ratified, and later changes go through this D1 owner.
- **Workplan decisions.** Stakeholder decisions D-1 to D-25 are recorded in
  the workplan §2; this document is their ratified D1 form.
- **Sole remaining gate.** The S-12 coordination-rule fallback is decided at
  qualification against the stated criterion.

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
closer than the bound melt onset distance (§6), and no salt ion within
T_min = 3.0 Å of framework Si/Al (§6, S-13); every ion anchored to an
opposite-charge first neighbour, and every α-cage salt cation's first-shell
Cl coordination between a conservative fraction a = 0.4 of a
density-projected melt value and a bound cap (§7); a cell volume from the bound contraction model (§9). Ring cations are
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
  - **Confirmed reading (D1 owner, 2026-10-18, via D2 acknowledgement).**
    The tie rule and nearest-centre membership apply to `cage_salt` ions.
    Ring cations are checked per (role, family, species) and site: every 8R
    site lies exactly between α₀ and α₁. Window weights are fixed per
    physical site. See D2 `v7_16_cage_salt_numerical_method.md` R-1 and R-2.
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
- **Criterion.** For every cage, |q_c| ≤ q_max(c). It is checked on the final
  written structure (§11).
- **Rebalancing (proposed mechanism; §11 step 4).**
  - After smoothing and before the salt build, `occlusion_counterion` cations
    are re-sited to empty sites of their own ring family, to bring every
    |q_c| ≤ q_max(c). This is the D-12 move set.
  - Later D-12 relocations of any ring cation are accepted only if they keep
    every |q_c| ≤ q_max(c).
  - Per-cage salt inventories are input (D-1) and are never changed.
    **No inter-cage moves (stakeholder decision S-10, 2026-10-17; paraphrase
    of "don't move an ion from an alpha to a beta cage or vice versa. You can
    move them around as long as it stays inside the cage"):** a `cage_salt`
    ion never moves from an α- to a β-cage or back, and never between cages.
    It may move freely inside its own cage. Applying S-10 to `cage_salt`
    ions only is an agent reading (proposed). Ring-site cations move only
    between ring sites under D-12 and step 4, which can change their sharing
    weights.
- **Failure handling.**
  - *Up-front rejection (proposed; deliberately stricter than necessary).*
    The input is rejected, naming the cages, if no assignment of the requested
    counterions to any sites of their requested ring families (best case)
    reaches |q_c^exc| ≤ q_max(c) on the excess charge alone. A favourable
    q^fw could offset the excess for some seeds, but the input is required to
    be balanceable without relying on that. Note
    that 8R moves cannot shift charge between α₀ and α₁, because every 8R
    site joins both.
  - *Seed-level failure* is a construction failure for that seed. It occurs
    when framework occupancy, which depends on the seed through the legacy
    stream, blocks the needed sites. Example from review 2: α₀ with 12 Cl,
    no cage cations and 12 counterions on 4R needs all 12 of α₀'s 4R sites,
    but on seed 0 two of them hold framework Na. The report names the blocking
    sites.
- **Threshold q_max (S-9, 2026-10-17: the stakeholder delegated the choice to
  the agent's reading of the charge-fluctuation literature).**

  ```
  q_max(c) = z_max(c) * e
  ```

  This gives 1 e for cages touched only by alkali cations, and 2 e for cages
  touched by Sr²⁺ or Ba²⁺.
  - **Basis.** Molten salts screen charge perfectly: S_ZZ(k) → 0 as k → 0
    (zeroth moment), and the Stillinger–Lovett condition fixes the k²
    coefficient (second moment). For large regions, the variance of the net
    charge then grows with the boundary area, not the volume (charge
    hyperuniformity). These are large-region sum rules and imply no bound at
    cage scale. For a region of cage size (radius about 1.5–2 ion spacings),
    whose boundary is mostly a wall the ions cannot cross, the agent's
    judgment is that a well-screened configuration should carry at most about
    one unscreened ion's charge. The rule allows at most one ion's worth of
    imbalance.
  - **Per cage (revised after review 3).** z_max(c) is the largest formal
    charge among the cations that contribute to cage c: its `cage_salt`
    cations and every ring cation with w(j, c) > 0. A single Sr or Ba
    therefore relaxes only the cages it touches. A cell-wide z_max would have
    doubled the tolerance everywhere.
  - Trade-off against S-5: framework Sr/Ba placement alone can produce ±2 e
    in the cages it touches. Only counterion re-siting can offset that,
    because q^fw is never rebalanced. The tolerance 2 e is accepted there for
    feasibility, not for screening.
  - The literature gives no number at cage scale. This threshold is an
    agent judgment under S-9, not an imported value.
  - **Literature (context; agent search, 2026-10-17; abstracts and snippets,
    full texts blocked by this environment's network policy):**
    - [P4: Zeidler et al., *J. Chem. Phys.* 157, 094504 (2022)] for molten
      NaCl;
    - Ph. A. Martin, T. Yalcin, *J. Stat. Phys.* 22, 435 (1980);
    - J. L. Lebowitz, *Phys. Rev. A* 27, 1491 (1983);
    - T. Hoang Ngoc Minh, B. Rotenberg, S. Marbach, *Faraday Discuss.* 246,
      225 (2023).
  - **Why a bulk-melt sphere statistic is not used (review 2, SC-A).**
    Counting ion centres in a cage-sized bulk sphere is dominated by ions
    cut by the arbitrary boundary. A perfect rock-salt lattice already gives
    σ_Q ≈ 1.35 e, and ⌈3σ_Q⌉ would allow about 5 e.
  - **Optional later calibration.** Run MD of salt-loaded LTA at 800 K with
    the stakeholder's MLFF. Compute the distribution of q_c with exactly this
    observable, and replace z_max by a stated upper percentile of |q_c| if
    the two differ materially (reopen trigger).
- **Reported:** every q_c^exc and q_c^fw, the rebalancing and relocation moves
  with their charge effect, and the distribution of q_c against q_max(c).

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
- **Framework Si/Al (stakeholder decision S-13, 2026-10-19; replaces the r7
  "no separate T-atom minimum").** Every `cage_salt` ion is at least
  **T_min = 3.0 Å** from every framework Si and Al. This interim value is the
  same for every salt species, Cl included.
  - **Basis (grade E).** The template's own crystallographic cation sites lie
    at least 3.05 Å from T (4R 3.05, 6R 3.23, 8R 3.05 Å). E-pair has no T rows.
    T_min is replaceable by bound per-cation M–T data when they are supplied.
  - **Why.** Under r7, small cations (mostly Li) passed every O minimum while
    sitting 1.5–2.6 Å from Si/Al, inside wall pockets. This happened in 28 of
    39 exploratory structures (F-E1-1, §13 T23).
  - Ring cations keep their template sites; T_min applies to `cage_salt` ions
    only.
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
n_M ≥ 1 `cage_salt` centres in an α-cage c, the coordination of those centres
must satisfy the rule below, evaluated per cage (density-reduced, per
S-2/S-4/S-7).
β-cage cations are excluded. Under D-2 they
have at most one Cl in their own cage (CN of about 1, as measured in review 1).
Their CN is reported separately and subject only to anchoring. Ring cations
count as neighbours of Cl (O4c, proposed yes), but are not CN centres. Only
cation-centred targets are imposed (D-13). Cl-centred counts are reported, not
targeted, because:
- they are coupled to the cation-centred counts by number conservation plus
  the ring contribution;
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

**Effective salt volume (stakeholder decision S-8, 2026-10-17).** The cage
volume is not well defined, because the interfacial structure between the
cage's inner surface and the salt is unknown. An interface of assumed
thickness is used, and the effective volume carries an uncertainty range.

```
d_wall(c)        = sum_M x_M(c) * [ r_peak(Cl–O around M) - r_peak(Cl–Cl around M) / 2 ]
V_eff,α(c; d)    = volume of the nearest-centre cell of alpha-cage c at distance > d from every framework O
n_melt(c)        = V_eff,α(c; d_wall(c)) / v_pair(800 K, x_Li(c))          (used by the §7 rule)
n_melt range(c)  = [ n_melt at d_wall(c) + 0.3 Å , n_melt at d_wall(c) - 0.3 Å ]   (reported)
```

- **Dividing surface (proposed, agent choice under S-8).** d_wall is the
  position of the dividing surface between salt and wall, measured from the
  O centres. It is not a closest-approach distance: salt ions may come
  closer, down to their r_cut.
- **Per-species values.** d_M comes from the E-pair r_peak columns, weighted
  by the cage's `cage_salt` cation mole fractions x_M(c). If there are none,
  the eutectic is used.

  | M | Li | Na | K | Rb | Cs | Sr | Ba |
  |---|---|---|---|---|---|---|---|
  | d_M (Å) | 1.39 | 1.87 | 1.97 | 2.12 | 2.30 | 1.72 | 1.85 |

  For the eutectic, d_wall = 1.62 Å, an interface about 0.22 Å thicker than
  the Shannon O²⁻ radius of 1.40 Å [P7].
- **Uncertainty ±0.3 Å.** This is half the Li–K spread of d_M (0.58 Å). It is
  an agent choice, and the range is reported, not used.
- **T-atom exclusion dropped from V_eff.** V_eff stays defined by the
  O-distance d_wall alone.
  - Under r7: Si/Al lie inside the O shell, so excluding them at a small T
    radius changes V_eff by 4 Å³ or less (review 3).
  - Under r8 (agent proposal, stakeholder-acknowledged 2026-10-19, S-14): the S-13
    T_min of §6 constrains ion *centres* only. V_eff measures the volume the
    melt occupies, ion bodies included.
  - Applying a 3.0 Å centre exclusion to V_eff would remove about 40 % of the
    α volume and 79 % of the β volume (Na-A fixture at f_V). It would cut
    n_melt to about 0.6×, which *raises* φ = min(1, N_Cl/n_melt) and so
    raises C_M and the floors. For example, at 12 Cl the K floor would go
    from 2 to 3. That would double-count the exclusion.
  - A T dividing surface built like the O one (T_min minus the Cl–O onset gap,
    about 1.8 Å from T) removes only about 1 % of the α volume. That is below
    the neglected ring-cation term (about 2 %) and the ±0.3 Å range (about
    ±8 %).
  - T_min removes no Cl-accessible α volume (319 Å³ either way). It removes
    only the small-cation wall pockets: Li-allowed α volume falls from 431 to
    302 Å³ (review of r8).
- **Ring cations not excluded.** Ring-cation volume is not subtracted (about
  −2 % if it were). This is conservative for the minimum.
- V_eff is computed on the realized, scaled cell with a deterministic
  estimator. Its tie rule on partition bisectors must avoid the −1 %
  aligned-grid bias found in review 3 (D2).

Na-A fixture at V_ref (agent Monte Carlo, 6×10⁵ points; three seeds agree
within 4 Å³):

| d (Å from O centres) | 1.20 | 1.40 | 1.62 | 1.80 | 2.00 |
|---|---|---|---|---|---|
| V_eff,α (Å³) | 895 | 868 | 810 | 745 | 667 |

| Composition | d_wall | n_melt (used) | Reported range (± 0.3 Å) |
|---|---|---|---|
| LiCl–KCl eutectic, x_Li 0.592 | 1.62 | 14.1 | about 12.5–15.3 |
| pure LiCl | 1.39 | 19.0 | about 17–20 |
| pure KCl | 1.97 | 9.2 | about 8–10.5 |

For comparison, the 11.4 Å [R31] and 11.05 Å IZA spheres give 13.5 and 12.3
Cl at x_Li 0.592. The d_M values for Rb and Cs (2.12 and 2.30 Å) lie above
the Li–K span, so cages rich in them carry an extrapolated d_wall.

**Finding (T17).** For the eutectic, the 800 K melt capacity of an α-cage is
14.1 Cl (range about 12.5–15.3). A fully loaded eutectic cage at the D-2 cap
of 12 is therefore at about 85 % of melt density. K-rich cages
(pure KCl: about 9) are dense at 10–12 Cl. The measured capacity of about
10.5 Cl per pseudo-cell [R1, R2], about 9.5 per α with 1 Cl in β, is
62–76 % of melt density (67 % at the central value).

**Coordination rule (stakeholder decision S-7, 2026-10-17).** "Coordination
is a distribution, not a strict number, so you should conservatively set a
minimum that is below the ideal. Say if the projected coordination is C(ρ),
then the minimum is set to a·C(ρ), where a is a number from 0 to 1."

```
C_M(c) = T_M^bulk * phi(c),   phi(c) = min(1, N_Cl(c) / n_melt(c))
for every cage_salt cation i of species M in alpha-cage c:
    max(1, ceil(a * C_M(c))) <= N_i(Cl) <= U_M          (a = 0.4)
```

- **Projection C(ρ) (grade D).** N(M–Cl) = 4π ρ_Cl ∫₀^{r_min1}
  g_MCl(r) r² dr, which is linear in the Cl number density if g is
  unchanged. C reaches T^bulk at melt density. n_melt uses the S-8 central
  dividing surface, and C over the S-8 range is reported.
- **Per ion (stakeholder decision S-11, 2026-10-17: reading B).** The minimum
  applies to every ion. Because a coordination number is an integer, the
  rule N_i ≥ a·C is the same as N_i ≥ ⌈a·C⌉. The anchoring floor of 1 always
  applies. This replaces the mean reading A, which review 3 had raised as an
  open choice.
- **a = 0.4 (stakeholder decision S-11).** The stakeholder chose a more
  conservative value than the agent's 0.7. Context for that 0.7: it was the
  median ratio of the lower edge of the experimental definition-dependent
  spread to the E-cn §1a value, and the simulated spreads would have given
  about 0.93. a remains a governed parameter, 0 < a ≤ 1.
- **Per-ion floors ⌈a·C⌉ at a = 0.4, eutectic (n_melt about 14):**

  | Cl per α | ≤ 6 | 8 | 10 | 12 | melt density |
  |---|---|---|---|---|---|
  | Li (a·C) | 1 (≤ 0.73) | 1 (0.98) | 2 (1.22) | 2 (1.46) | 2 (1.72) |
  | K (a·C) | 1 (≤ 0.89) | 2 (1.18) | 2 (1.48) | 2 (1.77) | 3 (2.08) |

  Sr and Ba at 12 Cl: a·C = 2.0 and 2.4, so each needs ≥ 3. A single dopant
  ion therefore faces the same floor as any other ion; no pooling is needed
  under B.
- **Flexibility.** Above the floor and below U_M the distribution is free.
  Clustered salt with higher local coordination remains allowed, consistent
  with the cluster counterevidence below.
- **Fallback (stakeholder decision S-12, 2026-10-17).** "If under this rule
  the structure generation is still too difficult to pass, then go back to A
  with a = 0.7."
  - Rule A is: mean_{i ∈ M, α-cage c} N_i(Cl) ≥ 0.7·C_M(c), plus anchoring
    and U_M, with small-population pooling. Pooling applies when a species
    has fewer than 3 centres in a cage: its mean is pooled over both α-cages,
    and if there are still fewer than 3 centres only anchoring and U_M apply.
  - **When the fallback is decided (proposed).** It is decided once, at D2/D4
    qualification. It is not a silent per-run switch, so every structure of
    a release has one coordination rule.
  - **Reference grid (proposed).**
    - Na-A and K-A rings;
    - seeds 0-4;
    - eutectic loadings of 2, 4, 6, 8, 10 and 12 Cl per α, with 1 Cl and 1 Li
      per β;
    - plus the three P5 dopant compositions.
  - **Pass criterion (proposed).** At least 95 % of grid requests build
    within the D2 restart budget, with every other rule met.
  - If rule B misses the criterion, the release default becomes rule A with
    a = 0.7, and the qualification record states the measured B pass rate.
  - The rule in force is written to every run report (§14).
- A physical reason for low Cl coordination of wall cations is coordination by
  framework O (for example Sr–O3 × 3 at 6R [Q25]). This is not modelled.

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
  the S-7 rule is one-sided (a lower bound only).

The linear form is an assumption (T16).

**Evaluation.** The rule is evaluated per α-cage c and per species M with ≥ 1
centres in c. C_M(c) uses only the cage's own Cl, but first neighbours may
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
spread bound, because the S-7 minimum goes below 1 at low loading.

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
   their ring family to satisfy |q_c| ≤ q_max(c), on its own stream derived from
   the seed. Seed-level failure is reported here;
5. volume model (§9): V_target depends only on the inputs, so the cell is
   scaled to V_target *before* the salt is built;
6. cage-salt build in the final cell, including D-12 relocation, on its own
   stream derived from the seed. A relocation is accepted only if it keeps
   every |q_c| ≤ q_max(c);
7. exact verification on the final scaled coordinates. It covers:
   - inventory (§4);
   - cage charge, |q_c| ≤ q_max(c) (§4);
   - hard minima (§6), with no tolerance below r_cut or T_min;
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

**Evidence predates S-13.** P5-P8 used 0.85 × Shannon sums for Si/Al, a
radius formula outside D-3. E1 and the r7 census used no T minimum. None of
them tests the r8 T_min; the re-frozen census does.

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
| P8 melt-density CN (S-4 reading, Option A, n_melt = 15) | Na/K-A × seeds 0-2 × n 2, 4, 6, 8, 12 × C2; 4 starts; targets max(1, T·n/15) | Hard minima met at all 30 points. Pooled α means lie within ±0.5, counting band edges as inside. Under spec rules 26/30 comply: 4 points (n = 2 three times, n = 6 once) leave one cation with CN 0, because the probe's smooth objective does not enforce anchoring, so D2 must enforce it explicitly. The probe pools both α-cages and uses n_melt = 15 throughout, whereas §7 is per cage (x_Li 0.5–0.67 would give 14.5–15.7 at d = 1.40 Å, and less at the r6 dividing surface). Recounted against the S-7 minimum, the mean condition holds at 30/30 points
for a = 0.6 and 0.7, and at 22/30 for a = 0.8. With all rules (anchoring, cap,
minima) it is 26/30 at a = 0.7. P8 shaped toward C itself, not toward the
minimum, so it is an existence check only. |
| Recount under S-11 (B, a = 0.4; per-cage x_Li and the r6 dividing surface; agent, after the fact) | stored P5-P8 results | Shaped toward density-reduced targets: P7 18/18 and P8 26/30; the 4 P8 misses are the isolated-ion anchoring cases. P6, shaped toward bulk targets: 26/36. Unshaped P5: 0/122. For comparison, a = 0.7 per ion gives P8 15/30 and P7 15/18. Per-ion floors are reachable with shaping, so D2 must shape CN and enforce anchoring. |
| P6 CN-shaped | Na/K-A × seeds 0-2 × n 8, 10, 12 × options A/B; 4 starts | All hard minima met. α-only, Option A: Li 3.9-4.5, K 5.1-5.5, all 18 points within ±0.5. Option B: Li 3.6-4.3, K 6.2-7.1; 3 of 18 points outside ±0.5, and 17 of 18 outside the source ranges. Cl-centred (cage + ring) 5.4-7.1. Anchoring min ≥ 1 (α). If β-Li (CN 1) were included in the mean, 6 of 18 Option A points would miss (review 1). |

**Variant disclosure (cumulative v7.16).**
- P1-P4 (previous session): about 90 points (workplan §5).
- P5: 132 points.
- P6: 36 points.
- P7: 18 points (one form, chosen after the stakeholder's S-2 instruction,
  not selected from alternatives).
- S-11 recount: stored P5-P8 results at a ∈ {0.4, 0.7}, per ion. a = 0.4
  was set by the stakeholder, not selected from the recount.
- P8: 30 points (the S-4 reading, one form). Recounted after the fact
  against a ∈ {0.6, 0.7, 0.8}. a = 0.7 was chosen from E-cn §1a data, not
  from this recount.
- V_eff(t): 5 thicknesses × 2 Monte Carlo seeds, then 5 distances × 1 seed
  plus 3 compositions for the r6 dividing surface. The dividing-surface rule
  was chosen after review 3 pointed to the E-pair derivation, not by
  comparing the n_melt it produces.
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
  - the §4 charge rebalancing and the q_c criterion with q_max(c) = z_max(c);
  - V_eff on ring forms other than Na-A, and t > 0 in any CN probe;
  - CN shaping with Cs/Sr/Ba;
  - the per-ion range rule, which was not a shaping objective. Recounted after
    the fact, 2 of 18 Option A P6 points have at least one α ion outside it;
  - β anchoring;
  - option (c) of §6;
  - the spec's own cage partition beyond 24 points;
  - realism or uniformity of the configurations;
  - misses of the CN minimum under a production restart budget;
  - a CN probe shaped directly toward the S-11 per-ion floors, and dopant
    (Cs/Sr/Ba) shaping;
  - P6 and P7 were recounted on the Li and K species only.

## 13. Tension assessment (agent-asserted, proposed until ratification)

Searched: `doc/specs`, `doc/archs`, `workplans/` (including archive), the
`qualification/` records, the repository's GitHub issues (none exist) and the
evidence package. Recorded items found: the v7.15 archived reopen triggers
(three fire, see Predecessor) and the v7.14 §9 limitations (consistent). No
recorded tension against v7.13 was found.

| ID | Tension | Assessment for this revision |
|---|---|---|
| T1 | Crowding at 10-12 Cl/α under rock-salt cut-offs | **Inapplicable to this revision.** The rock-salt bases are replaced by the E-pair onset r_cut, and P5 is feasible at 12/α. Under r8 (T_min), exploratory builds at 12/α still pass on start 1; the re-frozen S-12 census decides. Reopen if r_cut is recalibrated upward (E-pair §6.1, Zeidler NaCl partials). |
| T2 | Unsourced melt densities | **Superseded by T19.** The E-occ §5 handbook densities are not used. The §7 density reference uses the cited Janz correlations. |
| T3 | No ion–framework or dopant pair data | **Resolved by S-1 with an approximation.** Rb, Sr and Ba pairs are supplied, mostly grade E. Ion–framework minima use the free-oxide rows (grade E). The salt-ion–T minimum is S-13 (T23). Reopen if framework-O or M–T data arrive. |
| T4 | R10/R14 author lists vs v7.15 citations | **Resolved.** E-occ rev 2.1 lists R14 as Simpson & Gougar and R10 as Lewis, Fischer, Smith, matching v7.15. |
| T5 | Trivalent dopants outside v7.13 | **Applicable; out of scope** (§5). v7.13 is not reopened. |
| T6 | K–Cl CN is definition-dependent (5.2 vs 7.0) | **Applicable; Option A adopted** (consistent with the r_min1 boundary). |
| T7 | Cl-centred CN with ring neighbours (6-7) exceeds the bulk eutectic value (5.1-5.3) and confined-crystal values (3-4) | **Applicable.** Reported, not targeted. |
| T8 | Unshaped packing isolates ions | **Applicable** → anchoring rule (§7). D2 must shape CN. |
| T9 | BaCl₂ CN conflict (6.4 vs 7.7) | **Applicable.** E working value 7 ± 1. |
| T10 | Ca/Mg forms (5A) excluded | **Applicable** → §5 decision. |
| T11 | Contraction rests on one E-grade analogue point (K/Br-LTA at 4 Br per pseudo-cell, 294 K), with 3 of the 4 Br in β-cage Ag₄Br₄, so it is no loading equivalent for this model. The loading dependence is assumed: form (a) is linear to 10.5. The r3 form (c) was withdrawn after review 2 | **Applicable.** Reopen when matched 4A data or the user's DFT/MLFF f_V(n) arrive (E-vol §5b route ii). |
| T12 | Committed evidence file was rev 2.0. Rev 2.1 corrects the packing range (79-87 %, not 70-80 %) and R33/R32 conversions | **Resolved by committing rev 2.1.** The workplan T1 wording is superseded. |
| T13 | Bulk melt CN targets are undefined at low loading (review 1 SC1) | **Resolved by S-2, S-4 and S-7:** one-sided minimum a·C(ρ) plus anchoring. |
| T14 | The f_V bound (f_a ≥ 0.972) equals the v7.14 Li-A/Na-A ratio 0.9724, and 4A exchanges Na → Li in LiCl–KCl [R22, R32]. The R18 small-a phase may be partly an exchange effect, so f_V on a Li-rich V_ref could double-count (hypothesis, review 1) | **Applicable, open.** f_V is defined as occlusion-only (§9). Reopen with matched data. |
| T15 | Charge separation between α-cages is permitted without a per-cage rule (review 1) | **Resolved by S-5:** criterion abs(q_c) ≤ q_max with rebalancing (§4). Threshold per T18. |
| T16 | No study gives CN vs loading or density for a confined chloride melt. The linear-in-density form (§7) assumes uniform dispersal. Cluster evidence (CN 2–4 at low overall density) is counterevidence to applying it two-sided | **Applicable.** Reopen with an MD-computed CN(n) in LTA, for example from the user's MLFF (cf. E-cn §6.3). |
| T17 | The 800 K melt capacity of a eutectic α-cage is 14.1 Cl (range about 12.5–15.3), so cages at the D-2 cap are about 85 % dense and K-rich cages can be dense | **Applicable; handled** by the S-7 minimum at the S-8 central dividing surface, with the range reported. |
| T18 | No study quantifies cage-scale charge fluctuations, so q_max(c) = z_max(c) e is an agent judgment under S-9. A bulk-sphere statistic is the wrong estimand (review 2, SC-A) | **Applicable.** Reopen if salt-loaded LTA MD shows a materially different q_c distribution. |
| T19 | Melt density rests on Janz correlations read via citing papers (grade C), ideal mixing, and no dopant volumes | **Applicable, minor.** The check reproduces the cited eutectic value to 10⁻⁴ g cm⁻³. |
| T20 | The r3 q_c excluded framework cations, hiding a ±1 e local imbalance (review 2) | **Resolved:** q_c^fw added (§4). |
| T21 | The dividing surface d_wall (E-pair derived) and its ±0.3 Å range are agent choices. The interfacial structure is unknown (S-8). Rb- and Cs-rich cages are extrapolated | **Applicable.** Reopen with MD density profiles of salt in LTA. |
| T22 | Whether the S-7 minimum applies to the mean (A) or to each ion (B). The basis of a = 0.7 is the spread between definitions, a stand-in (review 3) | **Resolved by S-11:** per ion, a = 0.4, stakeholder-set. |
| T23 | F-E1-1 (AI agent, exploratory CLI campaign E1, 2026-10-18): with no T minimum, salt cations (mostly Li) sat 1.5–2.6 Å from Si/Al in 28 of 39 structures, satisfying every O minimum (`workplans/V7.16-D4-IMPLEMENTATION.md`, E1 findings) | **Resolved by S-13:** T_min = 3.0 Å for every `cage_salt` ion (§6), stakeholder-set, interim. |

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
- Per-species CN mean, histogram and the S-7 minimum, plus Cl-centred CN.
- Anchoring check.
- Restarts and failure reasons.
- Seed and derived stream identities.
- Bound evidence-file SHA-256s.
- Volume inputs: n̄, f_V, V_ref, and whether the factor was overridden.
- Per-cage q_c^exc and q_c^fw, z_max(c) and q_max(c), the rebalancing moves,
  and any blocking sites on seed-level failure.
- Per α-cage: x_Li(c) (or the fallback used), the dopant-validity flag,
  d_wall(c), V_eff,α and n_melt at d_wall and ± 0.3 Å, φ(c), C_M(c) over that
  range, the per-ion floor ⌈a·C_M(c)⌉, the achieved per-ion minimum and
  histogram, and whether the floor or anchoring was binding.
- α and β CN reported separately, with the CN option, a, and the per-ion
  floor and cap in force.
- The x_M weights and the resulting Cl–Cl r_cut.
- Ion–framework slack per class, including the salt-ion–T (Si/Al) minimum and its slack (S-13).
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
- Is charge locally balanced?
- Was the CN minimum binding, and how close was it?

**Falsifiers and reopen triggers.**
- MD shows systematic expulsion or crystallization traceable to the initial
  minima or CN.
- Bound per-cation M–T (Si/Al) minimum data supplied (replaces the interim
  S-13 T_min).
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

**Stakeholder decisions, third round (human, via chat, 2026-10-17).**
- **S-7.** The coordination minimum is a·C(ρ), 0 < a ≤ 1, applied below the
  ideal because coordination is a distribution (§7). The agent first set a = 0.7 (superseded by S-11: per ion, a = 0.4)
  (proposed default, from E-cn §1a).
- **S-8.** Cage volume is not well defined. Assume an interface thickness so
  that the effective salt volume carries an uncertainty range (§7: t ∈ [0,
  0.6] Å, proposed).
- **S-9.** q_max is to be picked by the agent from its charge-fluctuation
  study (§4: z_max e; made per cage, z_max(c), after review 3).
- **S-10.** Ions are not moved between α- and β-cages. Moves within a cage
  are allowed (§4).

**Superseded by S-7.** The three readings and the two-sided band.

**Revised after review 3** (agent choices under S-8 and S-9):
- the S-8 central dividing surface d_wall(c) from E-pair, with the range
  reported, replaces the t = 0 end;
- q_max per cage, z_max(c), replaces the cell-wide value;
- small-population pooling (retained only for the S-12 fallback to reading
  A);
- the up-front charge check is labelled as deliberately stricter;
- S-10 applies to `cage_salt` ions only.

**Stakeholder decision, fourth round (human, via chat, 2026-10-17).**
- **S-11.** "I was intending for B but your result also speaks against it.
  Try to go with B but with a more conservative a = 0.4." The per-ion reading
  (B) is adopted with a = 0.4 (§7).

- **S-12.** If generation under S-11 is still too difficult to pass, go back
  to reading A with a = 0.7. Recorded as a qualification-time fallback with a
  proposed reference grid and pass criterion (§7).

**Stakeholder decision, fifth round (human, via chat question, 2026-10-19).**
- **S-13.** Asked how salt ions should be kept off framework Si/Al after
  finding F-E1-1, the stakeholder chose "Interim 3.0 Å (Recommended)": one
  minimum of 3.0 Å from any Si/Al for every salt cation and Cl. Per-cation data
  may replace it later (§6).

- **S-14.** The stakeholder acknowledged the agent proposal to keep V_eff
  defined by the O distance alone under S-13 (§7): "Keep V_eff O-defined, I
  acknowledge." (chat, 2026-10-19).

**Re-qualification.** Census run 3 on the r8 code retained rule B (90/90).

**Open.** None.

## 16. D1 → D2 hand-off (on ratification)

D2 must preserve:
- §4 exact inventory and invariants;
- §6 hard minima with exact final verification (no tolerance below r_cut);
- §7 anchoring (enforced explicitly; see P8);
- the S-11 per-ion rule max(1, ⌈a·C(ρ)⌉) ≤ N_i ≤ U_M for every α-cage
  `cage_salt` cation, with a as a governed parameter (default 0.4);
- the S-12 fallback: measure B's pass rate on the §7 reference grid at
  qualification, and switch the release default to A with a = 0.7 if it falls
  below 95 %. A must therefore also be implementable, with pooling;
- n_melt from V_eff,α at d_wall(c) on the realized cell, with the x_Li
  fallback and the dopant flag;
- §4 S-10: packing and relaxation must preserve each `cage_salt` ion's
  nearest-centre cage membership;
- q_max(c) = z_max(c) e, and the step-4 moves that re-site counterions after
  smoothing (they partly undo it; reported);
- §4 the cage charge criterion (q_c^exc + q_c^fw), the rebalancing move set,
  relocation acceptance under the criterion, and seed-level failure
  reporting;
- a deterministic V_eff,α estimator with a stated tolerance and a tie rule
  on partition bisectors, so output stays
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
- how the CN miss probability is bounded. (Confirmed deviation, D1 owner,
  2026-10-18: D2 bounds it only by the S-12 census, not by a probability
  bound.)

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
- **SC1 (low loading).** Resolved by S-2/S-4/S-7 (one-sided minimum a·C(ρ) plus anchoring, §7).
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

**Independent falsification review 3** (2026-10-17) covered the r4→r5 delta
and a full-document pass. AI reviewer: a fresh subagent context, with the
same common-mode caveat.

**SERIOUS CHALLENGE.**
- **SC-R3-1.** Mean vs per-ion reading of S-7 was recorded as decided →
  returned to the stakeholder as §7 [CONFIRM].

**Findings and dispositions.**
- **Major:**
  - the basis of a = 0.7 (wrong estimand, experimental vs simulated spread)
    → disclosed, with the alternative value 0.93;
  - the t = 0 end left S-8 without effect, plus a closest-distance
    misdescription → central dividing surface d_wall(c), corrected wording;
  - cell-wide z_max → per cage;
  - stale §2 claim → restated;
  - single-centre strict floor → small-population pooling.
- **Minor:**
  - the up-front check is not a necessary condition → labelled as
    deliberately stricter;
  - S-10 scope and membership preservation → marked, and added to §16;
  - the T-atom term had no effect, and the ring-cation note was lost →
    dropped and restored;
  - numeric and wording errors → fixed (the 14.5–15.7 range, the sphere
    values at x 0.592, the anchoring thresholds, the density percentage, the
    O-radius convention, "radius" not "across", the sum-rule wording,
    number-conservation wording, per-cage wording);
  - stale text → fixed in the spec, workplan and hand-off;
  - §14 and §16 omissions → added.

Review 3 verified as correct:
- the V_eff table;
- the densities and n_melt;
- the a·C examples;
- the a ratios and U_M;
- the P8 recount;
- the fixture q^fw;
- the sharing distances;
- the edge cases;
- the attribution and lifecycle.

**Further review.** Changes after review 3 respond only to its findings. The
stakeholder's S-11 answer selected listed reading B and set the governed
parameter a = 0.4. The pooling rule now applies only within the S-12
fallback to reading A. A further independent review is therefore not required unless
ratification changes content beyond that.

**r8 amendment (2026-10-19).** It adds stakeholder decision S-13 (§6, §15) and
T23 (§13), and records V_eff's treatment of T_min as an agent proposal (§7).
Basis: the AI agent's exploratory finding F-E1-1. A focused independent review
of the r8 and D2 r3 amendments (AI agent, 2026-10-19) found no SERIOUS
CHALLENGE and no blocker, and returned "acceptable with fixes":
- **M1:** the census was not marked void; now stated in the header, D2 r3 and
  §15.
- **M2:** sign error in the §7 V_eff rationale; corrected, with its measured
  support.
- **Minors:** §2/§11/§12/§13/§14 wording, governance of the V_eff proposal, D2
  F_T definition, parent line, T23 search, salt–T histogram, margin
  acknowledgement; all folded in.
