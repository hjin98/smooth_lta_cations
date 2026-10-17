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

Unchanged authorities: `v7_13_group_i_ii_cation_extension.md` (supported
cations, formal charges, Q_fw = +24), `v7_14_group_i_ii_volume_reference_extension.md`
(dry reference volume V_ref), and the dry (no-occlusion) generator path, which
must remain output-identical.

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
approximately at bound melt values (§7); a cell volume from the bound
contraction model (§9). "Amorphous" means only: positions are random subject
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
CN data to 1100 K or T_m; the contraction value to 650–800 K in situ. The
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
- **Cage assignment (proposed, O8).** An ion belongs to the cage whose centre
  is nearest under the minimum-image Cartesian metric (nearest-centre
  partition). Ties are a construction failure. **[DECIDE]** alternative: a
  sphere of radius 5.7 Å about each α centre (11.4 Å free diameter [R31]).
  Probe P5 found feasibility insensitive to this choice.
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
sum_i N_i^salt(beta_k) <= 1               (proposed; D-6 "MCl unit")
```

**Invariants (independent, never collapsed).**

```
Q_fw  = sum_i z_i N_i^fw = 24                                    (v7.13)
N_Cl  = sum_c N_Cl(c)
Q_exc = sum_i z_i ( N_i^occ + sum_c N_i^salt(c) ) = N_Cl          (generalises v7.15)
```

There is no per-cage neutrality requirement (proposed; D-5, D-6): a β-cage Cl
may be compensated by a cation in the β-cage or on a ring site. Ring-site
capacity sum_i (N_ri^fw + N_ri^occ) ≤ C_r holds as in v7.15.

The final structure is a construction failure unless every role, cage and
species count equals its input exactly, with cage membership evaluated on the
written coordinates.

## 5. Supported salt species

S = {Li, Na, K, Rb, Cs, Sr, Ba}: the v7.13 cations for which E-pair supplies
every needed pair. When the cage-salt model is active, **every cation in the
structure, in any role, must be in S**, because every cation is a potential
neighbour of a cage Cl. Any other species (H, Fr, Be, Mg, Ca, Ra) is rejected
naming the missing pair; there is no radius fallback (D-3).

Consequence: Ca-containing frameworks (zeolite 5A) and Mg forms cannot be
salt-loaded until a pair table for them is supplied. **[DECIDE]** accept this
exclusion, or supply Mg/Ca rows (E-pair §2.2 holds measured Mg–Cl 2.42 Å and
Ca–Cl 2.78 Å peaks, but no recommended rows).

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

- **Weights x_M (proposed).** Mole fractions over all `cage_salt` cations of
  the cell. If there are none, over all cations. **[DECIDE]** alternatives:
  the Cl's own cage; or the conservative max over the cations present.
- **Ring cations get the same table values as any other ion** (D-12). No
  reduced or scaled cut-off exists in this model. The v7.16 draft's earlier
  scaling factor s (workplan O2) is retired, because the table's r_cut column
  is already a closest-approach onset (E-pair §0, §4.3).
- **Ring–ring pairs** (both on ring sites) keep the template geometry and are
  not constrained by this document. Their slack against the Cl-bridged r_cut is
  reported.
- **Ion–framework minima: [DECIDE, blocking].** E-pair §6.4 says to use
  framework-oxygen values for cation–framework contacts, not the melt table.
  The evidence gives only crystal Ag–O 2.49–2.52 Å [R27, R28] and Sr–O3
  2.21–2.31 Å [Q25]: nothing for Li/Na/K/Rb/Cs/Ba–O(framework) or for
  Cl–O(framework). D-3 forbids a radius fallback. Options:
  (a) the stakeholder supplies an ion–framework table;
  (b) the stakeholder ratifies, as an explicit grade-E proxy, the E-pair free
      oxide rows: M–O r_cut for cation–O(framework) and the weighted "Cl–O
      around M" r_cut for Cl–O(framework). Check: Sr–O 2.09 Å ≤ crystal
      Sr–O3 2.21 Å;
  (c) keep only the v7.15 catastrophic guard.
  Framework Si/Al: no data. Proposed: no separate T-atom minimum beyond the
  catastrophic guard below.
- **Catastrophic guard (retained from v7.15).** Every `cage_salt` ion is at
  least 1.5 Å from every atom. This is not a contact distance.

A structure violating any hard minimum is never written (fail closed).

## 7. Anchoring and coordination targets

**Anchoring (hard, proposed).** Every Cl has ≥ 1 first-neighbour cation in any
role, and every `cage_salt` cation has ≥ 1 first-neighbour Cl. Basis: halides
in LTA are anchored to cations and never seen isolated [E-occ §4]; D-5 relies
on coordination in place of a binding rule. Probe P5 showed that unconstrained
random packing leaves isolated cations (CN 0), so this needs to be stated.

**Coordination target (cation-centred).** For each `cage_salt` cation species
M with n_M ≥ 1 centres, the mean of N_i(Cl) over those centres must lie in the
band T_M ± δ_M. Ring cations count as neighbours of Cl (O4c, proposed yes), but
ring cations are not CN centres. Only cation-centred targets are imposed
(D-13). Cl-centred counts are reported, not targeted, because:

- by number conservation they follow from the cation targets plus the ring
  contribution;
- with ring cations counted they exceed bulk values (finding T7).

**Bound target values [DECIDE between two options].**

| M | Option A (proposed): E-cn §1a, defined inside `r_min1` | Option B: eutectic E-cn §1c |
|---|---|---|
| Li | 4.3 ± 0.5 | 4.2 (4.1–4.3), R_min unstated |
| Na | 5.2 ± 0.5 | unavailable → A |
| K | 5.2 ± 0.5 | 7.0 (6.8–7.2), definition-dependent |
| Rb | 5.4 ± 0.5 | 5.8 (dilute, XAFS) |
| Cs | 5.5 ± 0.5 | unavailable → A |
| Sr | 6 ± 1 (E, dilute working value) | same |
| Ba | 7 ± 1 (E; 6.4 vs 7.7 conflict open) | same |

Option A is proposed because it is the only set whose integration boundary is
the same `r_min1` the generator counts within (E-cn §0: "use §1a values as
targets inside r_min1"). Option B matches the base-salt composition but uses
an unstated boundary. Probe P6 reached both options within about ±0.5 (§12).

- **Distribution form.** The evidence supplies means ± spread, not P(CN)
  histograms (O4b answered by data availability). Only the mean is targeted.
  The spread arises from randomness, satisfying D-9, and is reported.
- **Failure policy (proposed).** Missing the band after the D2 restart budget
  is a construction failure (fail closed), as is a hard-minimum violation.
  **[DECIDE]** alternatively write the structure and flag the CN miss in the
  report.

## 8. Ring cations during the salt build (D-12)

A ring cation (role `framework` or `occlusion_counterion`, O1b proposed yes)
that conflicts with the cage salt may move only to:

- an empty physical site of its own ring family; or
- another positional variant of its own physical site (O1a, proposed yes).

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

**Bindings (proposed, all grade E).**

- f_ref = 0.97 ± 0.03, with hard bounds 0.919 ≤ f_ref < 1 (E-vol §5).
- n_ref = 10.5 Cl per pseudo-cell, the calibration loading of 4A + LiCl–KCl
  at 823 K [R1, R2] (E-occ §6.4).

**Form.** The evidence has one loaded point and no loading dependence. It
does show the contraction grows with occlusion and that the fully loaded phase
is the smallest [R18]. Linear interpolation from the dry endpoint, constant
above n_ref, is therefore the minimal monotone model. No dependence on
composition or cation form is modelled, because there are no data for one.

**Override [DECIDE].** Proposed: an explicit `--volume-factor` *replaces*
f_V(n̄), so V_target = volume_factor · V_ref, exactly the v7.14 meaning. The
alternative is that it multiplies f_V(n̄).

The legacy model and the dry path keep the v7.14/v7.15 behaviour (default
factor 1.0). The temperature mismatch (in-situ 650–800 K contraction applied
to a dry reference of another temperature) is below the stated uncertainty:
the thermal term is about 0.3 % in V (E-vol §4).

## 10. Legacy sparse model (D-10)

The v7.15 sparse model is reachable only through an explicit legacy selector.
For the same inputs and seed its output must be byte-identical to v7.15.0. The
default occlusion model is the cage-salt model. v7.15 flags keep their v7.15
per-cell meaning.

**[DECIDE] (O7).** Proposed: the v7.15 flags (`--alpha-cl`, `--beta-cl`,
`--occlusion-n4/6/8`, `--occlusion-min-distance`) are accepted only with the
legacy selector, and are rejected otherwise with a message naming the
selector. Old commands then fail loudly instead of being reinterpreted. New
per-cage inputs use new names (D3).

## 11. Construction order and determinism (D2/D3 constraints)

1. framework ring placement (legacy stream);
2. `occlusion_counterion` ring placement (v7.15 secondary stream);
3. cation density smoothing (ring cations only; Cl and `cage_salt` excluded);
4. cage-salt build, including D-12 relocation, on its own stream derived from
   the seed;
5. exact verification on the final coordinates: inventory (§4), hard minima
   (§6), anchoring and CN (§7);
6. volume model (§9);
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
"Feasible" means the exact final check passed with slack ≥ −1e-3 Å.

| Probe | Variants | Result |
|---|---|---|
| P5 main | ring Na/K/Li-A × seeds 0-2 × n = 8, 10, 11, 12 per α × relocation on/off; 8 starts | Relocation on: 24/24 at every point. Fixed rings: 19-23/24. |
| P5 sensitivity (Na-A, seed 0) | F1/F2 × C1/C2 × β B0/B1a/B1b × n 10, 12 × f_V 1.0/0.97, plus O1a off | 8/8 at all 50 points |
| P5 dopants | {Li5 K4 Cs1}, {Li5 K3 Sr1}, {Li6 K4 Ba1} × 3 seeds | 9/9 points 8/8 |
| P5 Ca ring | Ca-A | rejected: missing pair (Ca, Cl), as §5 requires |
| P5 unshaped CN | all feasible points | Li–Cl mean 1.9-2.1, K–Cl 3.5-3.7, min 0 (isolated ions); ring–ring slack ≥ +0.27 Å |
| P6 CN-shaped | Na/K-A × seeds 0-2 × n 8, 10, 12 × options A/B; 4 starts | all hard minima met. A: Li 3.9-4.5, K 5.1-5.5. B: Li 3.6-4.3, K 6.2-7.1. Cl-centred (cage + ring) 5.4-7.1; anchoring min ≥ 1 |

**Variant disclosure (cumulative v7.16).**
- P1-P4 (previous session): about 90 points (workplan §5).
- P5: 132 points.
- P6: 36 points.
- Every evaluated point is reported. Nothing was selected as a default from
  these results. All probe parameters were chosen by the agent, as were the F1
  proxy, the 5.7 Å sphere, the smooth CN surrogate width of 0.12 Å and the
  restart counts.
- Not examined: mixed or Group-II ring populations, more than 3 seeds,
  dopant-rich salts, n < 8, realism or uniformity of the configurations, and
  CN-band misses under a production restart budget.

## 13. Tension assessment (agent-asserted, proposed until ratification)

Searched: `doc/specs`, `doc/archs`, `workplans/` (including archive), the
`qualification/` records, the repository's GitHub issues (none exist) and the
evidence package. Recorded items found: the v7.15 archived reopen triggers
(three fire, see Predecessor) and the v7.14 §9 limitations (consistent). No
recorded tension against v7.13 was found.

| ID | Tension | Assessment for this revision |
|---|---|---|
| T1 | Crowding at 10-12 Cl/α under rock-salt cut-offs | **Inapplicable to this revision.** The rock-salt bases are replaced by the E-pair onset r_cut, and P5 is feasible at 12/α. Reopen if r_cut is recalibrated upward (E-pair §6.1, Zeidler NaCl partials). |
| T2 | Unsourced melt densities | **Inapplicable.** Densities no longer enter any binding; they remain only in the E-occ §5 order-of-magnitude packing argument. |
| T3 | No ion–framework or dopant pair data | **Partly resolved.** Rb, Sr and Ba pairs are supplied, mostly grade E. Ion–framework data are still absent: blocking decision in §6. |
| T4 | R10/R14 author lists vs v7.15 citations | **Resolved.** E-occ rev 2.1 lists R14 as Simpson & Gougar and R10 as Lewis, Fischer, Smith, matching v7.15. |
| T5 | Trivalent dopants outside v7.13 | **Applicable; out of scope** (§5). v7.13 is not reopened. |
| T6 | K–Cl CN is definition-dependent (5.2 vs 7.0) | **Applicable** → §7 decision. |
| T7 | Cl-centred CN with ring neighbours (6-7) exceeds the bulk eutectic value (5.1-5.3) and confined-crystal values (3-4) | **Applicable.** Reported, not targeted. |
| T8 | Unshaped packing isolates ions | **Applicable** → anchoring rule (§7). D2 must shape CN. |
| T9 | BaCl₂ CN conflict (6.4 vs 7.7) | **Applicable.** E working value 7 ± 1. |
| T10 | Ca/Mg forms (5A) excluded | **Applicable** → §5 decision. |
| T11 | Contraction rests on one E-grade analogue point (K/Br-LTA), linear in loading by assumption | **Applicable.** Reopen when matched 4A data or the user's DFT/MLFF f_V(n) arrive (E-vol §5b route ii). |
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
record above, written beside each structure. **Beyond the stakeholder's
explicit request; binds only if the stakeholder accepts it at ratification.**
Otherwise these questions are answered through agent reports only.

## 15. Decisions requested at ratification

1. §6 ion–framework minima: (a) supply data / (b) ratify the free-oxide proxy / (c) guard only. **Blocking.**
2. §7 CN target set: A (proposed) or B; and the failure policy for a CN-band miss.
3. §5 exclusion of Ca/Mg/H/Be/Fr/Ra from salt-mode structures, or supply rows.
4. §6 Cl–Cl weighting population.
5. §9 n_ref = 10.5, and whether `--volume-factor` replaces (proposed) or multiplies.
6. O1a and O1b (§8): proposed yes to both.
7. O8 cage assignment: nearest-centre (proposed) or 5.7 Å sphere.
8. §4 the β cage-salt cation cap of 1, and no per-cage neutrality.
9. O7 legacy flag handling (§10).
10. §14 the per-run report as a product surface.
11. O6: confirm trivalent cations are out of scope.

## 16. D1 → D2 hand-off (on ratification)

D2 must preserve:
- §4 exact inventory and invariants;
- §6 hard minima with exact final verification (no tolerance below r_cut);
- §7 anchoring and the CN band on the mean;
- §8 relocation limits;
- §11 stream isolation;
- fail-closed behaviour.

D2 owns:
- the packing and relaxation method;
- CN shaping;
- the restart budget;
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
