# Chloride-salt occlusion in zeolite A (LTA): evidence file for the occlusion model

**Revision 2.1, 2026-10-17.** This revision follows two independent reviews; §7 lists the changes. The file is the literature basis for revising the occlusion model in `smooth_lta_cations`.

Every number below carries a reference tag `[Rn]` (full list in §8) and a verification grade:

| Grade | Meaning |
|---|---|
| **A** | Read directly from the paper's abstract or publisher/OSTI record |
| **B** | Read from the full text of a patent or report |
| **C** | Secondary: the number is reported by another paper or patent citing the original; original full text not checked |
| **D** | Our own derived calculation (formula in §5, reproducible) |

All values come from the cited sources, except the melt densities in §5, which are approximate handbook values flagged there. Where a bibliographic detail could not be confirmed, it is marked *(to confirm)*.

---

## 0. Conventions

- **Pseudo-cell** = (AlSiO₄)₁₂ = M₁₂Al₁₂Si₁₂O₄₈ for monovalent M.
  - Cubic *Pm-3m*, a ≈ 12.3 Å.
  - Contains **one α-cage and one β-cage**.
  - Molar mass of Na₁₂Al₁₂Si₁₂O₄₈ = 1704.6 g/mol (D).
- **True cell of Na-A** is *Fm-3c*, a = 24.55 Å, ≈ Na₉₆Al₉₆Si₉₆O₃₈₄ = 8 pseudo-cells [R25].
  - The F-centred **primitive cell = 2 pseudo-cells**: 24 Al, 24 Si, 2 α-cages and 2 β-cages.
- **Loading units.** Experimental loadings in the literature are quoted **per pseudo-cell** ("per zeolite unit cell (AlSiO₄)₁₂"). Multiply by 2 for the primitive cell.
- **Charge balance.** Each occluded Cl⁻ arrives with a charge-balancing cation. A pseudo-cell at *n* Cl therefore holds 12 + *n* monovalent cations (cf. the Na₂₂ phase for 10 NaNO₃ [R11]).
- **wt% ↔ Cl/cell conversions depend on the mass basis** (§5):
  - Salt mass / (zeolite + salt): 12 Cl/cell = **28.1 wt%** LiCl–KCl eutectic.
  - Cl mass only / (zeolite + Cl): 12 Cl/cell = **20.0 wt%**.
  - INL's wt% is salt mass / (salt + zeolite). Their 20 wt% batch is 20.0 g LiCl–KCl + 80.0 g dehydrated 4A [R33, Table 1].
    - On that basis INL 20 wt% ≈ **7.7 Cl/cell** and 10 wt% ≈ **3.4 Cl/cell**.
    - [R33] nevertheless equates 20 wt% with the "theoretical occlusion limit of 12 Cl⁻ per pseudo unit cell". That equivalence holds numerically only on the Cl-only basis, so it is **an inconsistency in the source**.
    - **Do not convert INL loadings on the Cl-only basis.**
  - **State the basis whenever converting.**

---

## 1. Equilibrium loading data (calibration targets)

| # | System | T | Loading, Cl⁻ (or anion) per pseudo-cell | Grade | Ref |
|---|---|---|---|---|---|
| 1 | Zeolite 4A + LiCl–KCl–MeCl, Me = Na, Rb, Cs (x_MeCl 0–0.5, x_KCl/x_LiCl = 0.69) | 823 K | **11 ± 1** at infinite MeCl dilution; varies with MeCl content of the zeolite; no chemical reaction | A | [R1] |
| 2 | Zeolite 4A + LiCl–KCl–MeClₙ, Me = Ca, Sr, U (x up to 0.45) | 823 K | **10 ± 1** at infinite MeClₙ dilution, i.e. effectively 4A in LiCl–KCl (no Ca/Sr/U in the zeolite); irreversible reaction with UCl₃ | A | [R2] |
| 3 | ANL statement of zeolite capacity (LiCl–KCl waste salt) | not stated | **12** ("12 Cl⁻ per unit cell for zeolite") | B | [R5] |
| 4 | NaA + simulated IFR salt (92.5 wt% LiCl–KCl eutectic, 5 % NaCl, 1 % CsCl, 1 % BaCl₂, 0.5 % SrCl₂), washed | 673 ± 25 K | See below | B (data), D (conversion) | [R6] |
| 5 | Zeolite 4A in molten **pure LiCl** | not stated in abstract | **"between 10 and 11.5 salt molecules"**; basis (per pseudo-cell?) not stated in the abstract; conflicts with the rapid-decomposition reports in §3 | A | [R9] |
| 6 | Li-A (per [R6]) + LiCl–KCl eutectic | — | 3.5 meq Cl/g; LiCl occluded, KCl not. Mass basis not stated → **≈ 5.3** (per g LiA), **≈ 6.2** (per g product), or 6.0 / 7.0 if per g NaA | C (data), D (conversion) | [R7] via [R6] |
| 7 | Nitrate analogues, single-cation forms | nitrate melts (companion papers at 330–450 °C [R12, R13]) | LiNO₃ **11.7**, NaNO₃ **10**, AgNO₃ **10**, KNO₃ **≈ 0** per formula weight M₁₂[(AlO₂)₁₂(SiO₂)₁₂] | A | [R11] |
| 8 | Patent background statement: zeolite A (Na, K or Li form) + molten chloride waste salt | — | "up to about 25 % by weight" salt → **≈ 10.2** Cl/cell on the salt-mass basis | B (data), D (conversion) | [R8] |

**Row 4 details [R6]:**

- The patent reports 16.9–19.2 wt% Cl, which it calls equivalent to 4.8–5.3 meq Cl/g of product. That equivalence is the source's own arithmetic, and it does not quite hold: 19.2 wt% Cl is 5.42 meq/g.
- Cation and chloride analyses agree on **5.0 ± 0.3 meq salt/g**.
- From the Si analysis, zeolite is 66 wt% of the product.
- Converted on the Si basis: **≈ 12.9 Cl/cell** (12.1–13.7 for 5.0 ± 0.3 meq/g; 12.3–14.0 for 4.77–5.42 meq/g).
- The source has an internal mass-balance gap. 5.0 meq/g of IFR salt implies about 28 wt% salt, i.e. about 72 wt% zeolite, not 66 %.

**Rows 4 and 8 are not independent.** Row 8's statement in [R8] is background describing the same ANL work. The difference between 10.2 and 12.9 Cl/cell comes from the conversion basis and the mass-balance gap, not from different data.

**Rows 1 and 2 measure the same quantity.** Both are infinite-dilution limits for 4A in LiCl–KCl at 823 K, from two independent datasets, and they agree within error: **≈ 10.5 ± 1 Cl per pseudo-cell**. Neither shows that capacity depends on the cation. Both papers report that loading varies with MeCl/MeClₙ content in the zeolite, but the abstracts do not quantify it.

**Do not confuse with process loading.** The ceramic-waste-form process deliberately loads only **3.4–3.9 Cl⁻ per pseudo-cell (target 3.9)** [R4], or 3.8 [R5]. That figure is set by the sodalite reaction stoichiometry Na₁₂(AlSiO₄)₁₂ + 4 NaCl → 2 Na₈(AlSiO₄)₆Cl₂ [R4], not by the occlusion capacity. In the mixing study, keeping chloride below **80 % of the full salt capacity** reduced free salt [R17].

**Equilibration conditions in INL work:** 500 °C and 24 h contact, "to allow for equilibrium to be reached" [R15, ANS abstract; R33]. Their fitted model has separate ion-exchange and occlusion sites, with a fitted total occlusion capacity [R15].

At 10 and 20 wt% eutectic loading, 4A occluded > 99 % of the salt (99.4 % and 99.3 %) [R33]. **This validates near-complete uptake only up to about 7.7 Cl/cell.** It does not validate uptake at 12 Cl/cell.

**Above "capacity" ([R32]):** [R32] calls loadings above 20 wt% "super salt occluded zeolite". It says such loadings exceed "the reported capacity", from which the authors infer the zeolite is at capacity and the excess salt is free. The free salt is mainly NaCl, formed by Na⁺ exchange.

**Do not use 20 wt% as a capacity or free-salt onset in the model.** On the salt-mass basis 20 wt% is only about 7.7 Cl/cell, well below the 10–12 Cl/cell equilibrium data in rows 1–4. R32 inherits the INL 12 Cl ↔ 20 wt% equivalence noted in §0.

---

## 2. Cation and salt dependence

1. **Large cations reduce or block occlusion.**
   - KNO₃ is essentially not occluded in K-A [R11].
   - In K-A in nitrate melts, partial Na⁺NO₃⁻ occlusion starts only after 7 K⁺ have been exchanged for Na⁺, and small amounts of K⁺ block the channels [R13].
   - In [R13], occlusion starts at much lower Na content at 450 °C than at 350 °C. Blocking therefore weakens with temperature.
   - On Li-A in LiCl–KCl, LiCl but not KCl was occluded [R7 via R6].
2. **Nitrate occlusion is in the α-cage.**
   - Liquornik & Marcus place occluded nitrate in the α-cage. Where Li, Na or Ag nitrate is occluded, all cations can be exchanged; with other cations exchange is limited [R12].
   - Petranović & Krstanović describe NaA·NaNO₃ with nitrate in the "central α cage" [R26].
3. **Occluded and framework cation selectivities differ.**
   - Two-site model fitted to Cs–Rb–Na–Li–K data: Cs is preferred in the occluded (cage) sites, Li in the framework sites [R14].
   - Multivalent extension: Sr²⁺ is preferred in the ion-exchange site and Cs⁺ in the occlusion site; Li⁺ has the highest combined selectivity [R16].
4. **Capacity depends on the third component's concentration.**
   - Occlusion extent is a function of the MeCl (or MeClₙ) concentration in the zeolite. The values in §1 are infinite-dilution limits [R1, R2].
   - Separation factors relative to Li at infinite dilution:

     | Salt | Separation factor | Ref |
     |---|---|---|
     | NaCl | 0.84 | [R1] |
     | RbCl | 0.87 | [R1] |
     | CsCl | 2.31 | [R1] |
     | CaCl₂ | 0.9 | [R2] |
     | SrCl₂ | 11.9 | [R2] |
     | UCl₃ | 13 | [R2] |

**Gap:** no single-salt chloride melt equilibrium, such as NaCl in NaA or KCl in KA, was found. Pure NaCl and KCl melt at 801 °C and 771 °C. [R6] advises that "temperatures over about 600 °C are to be avoided since they may result in destruction of the zeolite". In practice, per-salt capacities therefore exist only as components of a LiCl–KCl solvent. That last point is inference.

---

## 3. Validity regime / framework stability

- **Occlusion pathway (LiCl–KCl + 4A).** Solid-state occlusion runs from about 450 K to 625 K. Abrupt molten-salt occlusion follows at 628 K (the eutectic melting point), accompanied by lattice contraction [R18].
- **UCl₃–LiCl–KCl + 4A.** LiCl/KCl begins diffusing into the lattice at 450 K. UCl₃ "does not appear to react with the zeolite 4A" up to the 685 K melting/occlusion step [R19]. In contrast, [R2] reports an irreversible reaction with UCl₃ at 823 K.
- **Temperature.** Mixing studies found occlusion improved above 775 K [R17]. Lexa's equilibria were at 823 K [R1, R2], the INL equilibria at 773 K [R15], and the [R6] data at 673 K.
- **4A → sodalite.** Salt-loaded 4A tends to transform to sodalite. Na/Ca 5A retains LTA more readily [R20].
- **Molten pure LiCl: conflicting evidence, unresolved.**
  - Seo et al. report 10–11.5 salt molecules occluded, and Cs/Sr uptake that plateaus after 2–4 h [R9]. Temperature is not stated in the abstract.
  - At 620–680 °C, "LTA structure in LiCl salt was decomposed to LiAlSiO₄ in 10 min"; Cs selectivity disappeared in 20–30 min. In LiCl–KCl, by contrast, a fraction of LTA remained after 60 min [R23].
  - 4A changed to Li-A when mixed with molten LiCl at 920 K [R22].
  - A patent states that in LiCl "the zeolite structure completely breaks down" (no temperature given) [R24].
  - **Treat pure-LiCl loadings as unverified and possibly transient.**

---

## 4. Structural evidence on halide siting

Diffraction on genuinely melt-equilibrated LTA is limited.

- **Neutron diffraction on salt-occluded 4A/5A (Li, K, Cs, Sr, Ba, Y chlorides) [R20].** Cations act both as charge compensators and as members of occluded clusters. The clusters roughly fill the 6–8 Å pores. Abstract only; site tables not checked.

Single-crystal X-ray studies give halide motifs inside LTA. **These come from vapour or solution routes, not melt equilibrium** (e.g. KCl or KI in methanol at room temperature; WCl₆ vapour at 453 K).

| Structure | β-cage content | α-cage content | Total halide / pseudo-cell (D, from formula) | Ref |
|---|---|---|---|---|
| \|K₂.₃₅Ag₁.₁(Ag₄Cl₄)₀.₄₅(C₃H₃O₃)₁.₁(K₃Cl)₃.₄₅(K₃(OH)₂)₀.₅₅(H₂O)≥₃\|[Si₁₂Al₁₂O₄₈]-LTA | Ag₄Cl₄ in about 45 % of sodalite cavities; each Cl bonds a K⁺ through a 6-ring (Cl–K 2.70 Å); Ag–O(framework) 2.493 Å | K₃Cl²⁺ units predominate | 0.45×4 + 3.45 = **5.25 Cl** | [R27] |
| K₉(K₄Br)·0.75 Ag₄Br₄-LTA | Ag₄Br₄ (T_d, about 7.9 Å) in about 75 % of sodalite units | near square-planar K₄Br³⁺ (see note) | 3 + 1 = **4 Br** | [R28] |
| \|K₉(K₄I)(Ag₄I₄)₀.₅\|-LTA | Ag₄I₄ in half of the sodalite cavities | "near square planar K₄I³⁺" (one extra KI per large cavity) | 2 + 1 = **3 I** | [R34] (A) |
| \|(WCl)₀.₆(Tl₇Cl₆)₀.₅Tl₈.₅\|-LTA | — | Half of the α-cages: Tl₇Cl₆⁺, a central Tl⁺ octahedrally coordinated by 6 Cl⁻, each Cl also bonded to an 8-ring Tl⁺. Other half: WCl⁵⁺ opposite 6-rings | 0.6 + 3 = **3.6 Cl** (**6 Cl in one α-cage**) | [R30] |

Note on the [R28] row: one Br⁻ per cell lies **opposite a 4-ring**, "held there by two eight-ring and two six-ring K⁺ ions".

Take-aways for siting:

- In the α-cage, halide is anchored to cations and is not seen as an isolated anion. Observed positions:
  - on the 4-fold axes, between a cage-centre cation and an 8-ring cation [R30];
  - **opposite a 4-ring**, bridging 8-ring and 6-ring cations [R28, R34];
  - on a 3-fold axis opposite a 6-ring, but only as the Cl of WCl⁵⁺ [R30]. This is indirect evidence.
- In [R27], the Cl bonded to K⁺ "through a 6-ring" belongs to the β-cage Ag₄Cl₄ cluster. It is not α-cage evidence.
- The β-cage can hold **more than one** halide when Ag⁺ is present (M₄X₄ clusters). For Na/K/Li, the sodalite-like [M₄Cl]³⁺ unit (1 Cl per β) remains the natural default, consistent with the sodalite reaction stoichiometry [R4]. This default is inference.

---

## 5. Derived numbers (Grade D), reproducible

Atomic masses: Na 22.990, Li 6.941, K 39.098, Al 26.982, Si 28.086, O 15.999, Cl 35.453, Cs 132.905.

LiCl–KCl eutectic composition follows [R1]: x_KCl/x_LiCl = 0.69, i.e. 59.2 : 40.8 mol%. INL uses 58 : 42 [R33]. The difference changes the results below by ≤ 0.7 % (relative).

| Quantity | Formula | Result |
|---|---|---|
| M(Na₁₂Al₁₂Si₁₂O₄₈) | — | 1704.6 g/mol (Si 19.77 wt%, matching 19.8 % in [R6]) |
| M(Li₁₂Al₁₂Si₁₂O₄₈) | — | 1512.1 g/mol |
| [R6] Cl per cell | q / (0.66 / 1704.6), q = 4.7 / 5.0 / 5.3 meq/g | 12.1 / **12.9** / 13.7 |
| [R7] Cl per cell, per g LiA | 3.5e-3 × 1512.1 | 5.3 |
| [R7] Cl per cell, per g product | 3.5e-3·M_LiA / (1 − 3.5e-3·M_LiCl) | 6.2 |
| [R7] Cl per cell, per g NaA / product (NaA basis) | same with M_NaA | 6.0 / 7.0 |
| Eutectic mean M (59.2 : 40.8) | 0.592·42.39 + 0.408·74.55 | 55.5 g/mol |
| wt% salt at 10 / 11 / 12 Cl per cell (salt-mass basis) | 55.5n / (1704.6 + 55.5n) | 24.6 / 26.4 / 28.1 % |
| wt% at 10 / 11 / 12 Cl per cell (Cl-only basis) | 35.45n / (1704.6 + 35.45n) | 17.2 / 18.6 / **20.0** % |
| Cl per cell at 25 wt% salt [R8] | w·M_NaA / ((1−w)·55.5) | 10.2 |
| Cl per cell at 20 wt% salt [R33] | same | 7.7 (salt-mass basis) |
| α-cage volume, 11.4 Å free diameter [R31] | (4/3)π(5.7 Å)³ | ≈ 776 Å³ |
| Ion pairs per α-cage at bulk-melt packing | V_α / (M / ρ_melt / N_A); ρ ≈ 1.50 (LiCl), 1.55 (NaCl), 1.52 (KCl), 2.79 (CsCl) g/cm³ | LiCl 16.5, NaCl 12.4, KCl 9.5, CsCl 7.7; eutectic (ideal mixing) 12.7 |

**About the melt densities.** They are approximate handbook values near the melting point, consistent with Janz-type correlations to about 0.01 g/cm³. They were **not** retrieved from a source for this file.

**About the packing estimate.** It is an order-of-magnitude size argument. The measured infinite-dilution totals of 10–11 Cl/cell [R1, R2] correspond to the following α-cage packing:

- 71–79 % of bulk-melt packing if 1 Cl sits in the β-cage (the §4 default);
- 79–87 % if all Cl sits in the α-cage.

- Packing alone orders the salts LiCl > NaCl > KCl > CsCl.
- The near-zero KNO₃ uptake in K-A is a window-blocking effect [R13], not a packing effect.

Geometry check for candidate α-cage Cl positions at radius r = 3.5 Å from the cage centre:

| Site pair | Nearest Cl–Cl (Å) |
|---|---|
| 8 sites on the 3-fold axes ([111]) | 4.04 |
| 6 sites on the 4-fold axes ([100]) | 4.95 |
| 3-fold site to 4-fold site | 3.22 |
| 4-ring-facing site ([110], 12 sites) to 3-fold site | 2.12 |
| 4-ring-facing site to 4-fold site | 2.68 |

The site sets are strongly mutually exclusive, so an exclusion rule is essential.

---

## 6. Implications for the model (recommendations, not evidence)

1. **β-cage.** Keep one Cl site at the cage centre as the default for Li/Na/K (sodalite-like [M₄Cl]³⁺). Make the capacity a per-cation parameter so Ag-type M₄X₄ clusters (up to 4 X per β) can be represented [R27, R28, R34].
2. **α-cage.** Replace the single cage-centre Cl site with a multi-site set. Candidate **Cl** sites:
   - 4-fold-axis sites near the 8-ring windows (6) [R30: Cl⁻ around a central Tl⁺, each also bonded to an 8-ring Tl⁺];
   - sites facing the 4-rings (12) [R28, R34];
   - 3-fold-axis sites opposite the 6-rings (8). Support is indirect only: the Cl of WCl⁵⁺ bound opposite a 6-ring [R30].

   **The α-cage centre is a cation site, not a Cl site.** In [R30] it holds a Tl⁺ coordinated by 6 Cl⁻. No cited structure has a halide at the α-cage centre. If a centre Cl site is kept anyway, mark it "no structural support" and exclude it from all shell sites within r (= 3.5 Å at the r used in §5).

   Enforce a Cl–Cl exclusion distance (§5 geometry), or use off-lattice / grand-canonical sampling, so that the capacity emerges rather than being hard-coded.
3. **Charge balance.** Insert or remove Cl⁻ together with a cation. Track the occluded-site cation identity separately from the framework sites [R14, R16].
4. **Calibration targets (4A in molten LiCl–KCl at 823 K):**
   - total ≈ **10.5 ± 1 Cl per pseudo-cell** at infinite dilution of the third component, from two independent datasets [R1: 11 ± 1; R2: 10 ± 1];
   - upper bound ≈ 12–13 [R5; R6 at 673 K].

   No source quantifies a cation-specific total capacity, so do not build one in. Use the separation factors in §2 and the INL two-site fits [R14–R16] to constrain the cation **partition** between framework and occluded sites. Do not use them to set the total.
5. **Cation-size penalty.** K⁺ (and larger) occupancy of the 8-ring windows can suppress α-cage loading: KNO₃ ≈ 0 in K-A [R11], and a few K⁺ block occlusion [R13]. The penalty should weaken with temperature [R13].

   **Calibrate it so it does not break target 4.** 4A in a 41 mol% KCl melt still occludes 10–11 Cl/cell at 823 K. The nitrate data [R11, R13] are at 330–450 °C.

---

## 7. Changes in revision 2 (from the independent review)

1. **K₄I³⁺ geometry restored.** Revision 1 called "near square planar" a snippet artefact; that was wrong. The original paper [R34] states it, and it is now cited directly at Grade A. The bromide analogue [R28] adds the observed 4-ring-facing halide site, now in §4 and §6.
2. **R14 authors corrected** to M. F. Simpson and M. L. D. Gougar (two authors).
3. **Pure-LiCl evidence presented as unresolved.** Revision 1 overstated the case against breakdown. [R23] reports decomposition to LiAlSiO₄ in 10 min at 620–680 °C. [R9]'s 10–11.5 figure has no stated basis or temperature. "Short-contact" was removed.
4. **Conversion basis made explicit.** INL equates 12 Cl/cell with 20 wt%, which is the Cl-only basis. [R8] is not independent of [R6]. The [R6] source arithmetic and its mass-balance gap are noted.
5. **R32 now quoted precisely, and R33 added** (INL loading, efficiency and eutectic composition).
6. **Eutectic composition** now taken from [R1] (59.2 : 40.8). The packing fractions are recomputed: 79–87 %, not "70–80 %". The over-reaching ordering claim has been removed.
7. **Smaller corrections:**
   - "process T" changed to "not stated" (row 3).
   - The temperature claims in §2 and §6 were corrected.
   - The UCl₃ reactivity contrast was added (§3).
   - "Aqueous" changed to "solution".
   - The R13 temperature effect was added.
   - The 4-ring site distances were added.
8. **Bibliographic details completed** throughout §8.

**Revision 2.1 (blocking-error review):**

- R1 and R2 are merged into one calibration target, ≈ 10.5 ± 1 Cl/cell. R2's 10 ± 1 is an infinite-dilution value for LiCl–KCl, not a Ca/Sr/U capacity.
- INL wt% is now stated as salt mass / (salt + zeolite), so INL 20 wt% ≈ 7.7 Cl/cell. Their "12 Cl ↔ 20 wt%" equivalence is flagged as a source inconsistency. 20 wt% must not be used as a capacity, and the > 99 % occlusion result validates uptake only up to about 7.7 Cl/cell.
- The α-cage centre is relabelled as a cation (Tl⁺) site in [R30]. The α-cage 3-fold Cl site is marked as indirect evidence. [R27]'s 6-ring Cl is reassigned to the β-cage cluster.
- Smaller fixes:
  - packing range now quoted with and without the β-cage Cl;
  - eutectic-composition effect corrected to ≤ 0.7 %;
  - K⁺-penalty calibration caveat added;
  - R15 source of the 500 °C / 24 h conditions given;
  - R32 details added, marked "to confirm".

---

## 8. References

The Status column records what was checked.

| Tag | Reference | Supports | Status |
|---|---|---|---|
| R1 | D. Lexa, I. Johnson, "Occlusion and ion exchange in the molten (lithium chloride–potassium chloride–alkali metal chloride) salt + zeolite 4A system with alkali metal chlorides of sodium, rubidium, and cesium," *Metall. Mater. Trans. B* **32**(3), 429–435 (2001), doi:10.1007/s11663-001-0028-4; OSTI 943218 (ANL/CMT/JA-38390) | 11 ± 1 Cl/cell, 823 K; separation factors | Abstract |
| R2 | D. Lexa, "Occlusion and ion exchange in the molten (LiCl + KCl + alkaline-earth chloride) salt + zeolite 4A system … and … with UCl₃," *Metall. Mater. Trans. B* **34**(2), 201–208 (2003), doi:10.1007/s11663-003-0007-z; OSTI 949686 | 10 ± 1 Cl/cell; separation factors; UCl₃ reaction | Abstract |
| R3 | D. Lexa, *J. Chem. Thermodyn.* (1999) *(details to confirm; cited in R1)* | Binary LiCl–KCl occlusion thermodynamics | Not found |
| R4 | US Patent 10,600,527, "Ceramic waste form production from used nuclear fuel." https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10600527 | Sodalite stoichiometry; 3.4–3.9 Cl/cell process loading | Full text (excerpts) |
| R5 | T. J. Battisti, K. M. Goff, K. J. Bateman, M. F. Simpson, J. P. Lind, "Ceramic waste form production and development at ANL-West," ANS 5th Topical Meeting, Charleston SC, 17–20 Sept 2002. https://corpora.tika.apache.org/base/docs/govdocs1/230/230873.pdf | 12 Cl/cell capacity vs 3.8 process loading | Full text (excerpts) |
| R6 | M. A. Lewis, T. R. Johnson, US Statutory Invention Registration H1227, filed 14 Aug 1991, issued 7 Sept 1993, "Method for immobilizing mixed waste chloride salts containing radionuclides and other hazardous wastes." https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/H1227 | NaA + IFR salt Cl analysis; >600 °C caution; prior-art summary | Full text |
| R7 | M. V. Šušić, N. A. Petranović, D. A. Mioč, "The properties of zeolite 4A treated in molten salts," *J. Inorg. Nucl. Chem.* **33**(8), 2667–2675 (1971), doi:10.1016/0022-1902(71)80244-0 | Li-A + LiCl–KCl: 3.5 meq Cl/g, LiCl only | Number secondary via R6; bibliography confirmed |
| R8 | M. A. Lewis, C. Pereira, US Patent 5,613,240 (1997), "Method of preparing sodalite from chloride salt occluded zeolite." https://patents.justia.com/patent/5613240 | "up to about 25 % by weight" (background, same work as R6) | Full text (excerpts) |
| R9 | C. S. Seo, B. H. Park, S. B. Park, K.-J. Jung, S. W. Park, S. H. Kim, "Study on the characteristics of the ion exchange of zeolite 4A in a molten LiCl system," *J. Chem. Eng. Japan* **39**(1), 27–33 (2006), doi:10.1252/jcej.39.27 | 10–11.5 salt molecules in molten LiCl | Abstract |
| R10 | M. A. Lewis, D. F. Fischer, L. J. Smith, "Salt-occluded zeolites as an immobilization matrix for chloride waste salt," *J. Am. Ceram. Soc.* **76**(11), 2826–2832 (1993), doi:10.1111/j.1151-2916.1993.tb04023.x | Foundational ANL paper | Record only |
| R11 | M. Liquornik, Y. Marcus, "Ion exchange in molten salts. II. The occlusion of lithium, sodium, potassium and silver nitrates in the respective forms of zeolite A," *Isr. J. Chem.* **6**(2), 115–121 (1968), doi:10.1002/ijch.196800017 | 11.7 / 10 / 10 / ≈ 0; Na₂₂ model | Abstract |
| R12 | M. Liquornik, Y. Marcus, "Ion exchange in molten salts. I. The ion-exchange properties of sodium zeolite A in molten NaNO₃; exchange reactions with alkali metal, thallium, and silver cations," *J. Phys. Chem.* **72**(8), 2885–2889 (1968), doi:10.1021/j100854a034 | α-cage occlusion; 330 °C (350 °C for K) | Abstract |
| R13 | M. Liquornik, Y. Marcus, "Ion exchange in molten salts. V. Potassium zeolite A as an ion exchanger in nitrate melts," *J. Phys. Chem.* **75**(16), 2523–2525 (1971), doi:10.1021/j100685a019 | K⁺ blocking; occlusion after 7 K → Na; lower threshold at 450 °C | Abstract |
| R14 | M. F. Simpson, M. L. D. Gougar, "Two-site equilibrium model for ion exchange between monovalent cations and zeolite-A in a molten salt," *Ind. Eng. Chem. Res.* **42**(18), 4208–4212 (2003), doi:10.1021/ie020896b; OSTI 961100 (ANL/ENT/JA-44603) | Cs prefers cage sites; Li prefers framework | Abstract and record |
| R15 | T.-S. Yoo, S. M. Frank, M. F. Simpson, P. A. Hahn, T. J. Battisti, S. Phongikaroon, "Salt-zeolite ion-exchange equilibrium studies for a complete set of fission products in molten LiCl-KCl," *Nucl. Technol.* **171**(3), 306–315 (2010), doi:10.13182/NT10-A10866. ANS abstract: https://wx1.ans.org/pubs/journals/nt/a_10866 | 500 °C / 24 h "to allow for equilibrium to be reached" (ANS abstract; not in the OSTI abstract); fitted total occlusion capacity | Abstract |
| R16 | S. Phongikaroon, M. F. Simpson, "Equilibrium model for ion exchange between multivalent cations and zeolite-A in a molten salt," *AIChE J.* **52**(5), 1736–1743 (2006), doi:10.1002/aic.10784; conf. version INL/CON-05-00460 | Sr prefers IX site; Cs prefers occlusion site | Abstract |
| R17 | C. Pereira, V. N. Zyryanov, M. A. Lewis, J. P. Ackerman, "Mixing of zeolite powders and molten salt" (ANL, 1996), abstract via https://worldwidescience.org/topicpages/s/salt+occluded+zeolite.html | > 775 K improves occlusion; < 80 % of capacity reduces free salt | Abstract |
| R18 | D. Lexa, L. Leibowitz, A. J. Kropf, "Occlusion of the LiCl-KCl eutectic salt by zeolite 4A by simultaneous differential scanning calorimetry/X-ray diffraction," APS Activity Report 2001. https://www-stg.aps.anl.gov/sites/www.aps.anl.gov/files/APS-sync/activity_reports/apsar2001/LEXAD1.PDF | 450–625 K solid-state, 628 K molten occlusion; contraction | Report (excerpts) |
| R19 | D. Lexa, L. Leibowitz, A. J. Kropf, "The interaction between UCl₃ and zeolite 4A by in situ temperature-resolved powder x-ray diffraction," APS report *(year not printed; filed under 1999)*. https://www.aps.anl.gov/sites/www.aps.anl.gov/files/APS-Uploads/APS-Science/1999/lexad1.pdf | Diffusion from 450 K; UCl₃ unreactive to 685 K | Report (excerpts) |
| R20 | J. W. Richardson Jr., "Salt-occluded zeolite waste forms: crystal structures and transformability," *MRS Proc.* **465**, 395–400 (symposium 1996; © 1997), doi:10.1557/PROC-465-395; OSTI 541788 | Neutron: cation clusters; 4A → sodalite, 5A resists | Abstract |
| R21 | M. A. Lewis, D. F. Fischer, J. J. Laidler, "Leach resistance properties and release processes for salt-occluded zeolite A," *MRS Proc.* **294** (1993), doi:10.1557/PROC-294-95 | Cites R7 | Record |
| R22 | J.-G. Kim, J.-H. Lee, J.-H. Yoo, J.-H. Kim, "Immobilisation of molten LiCl waste using zeolite A," 8th OECD-NEA Information Exchange Meeting, Las Vegas 2004, poster PS1-03. https://test.oecd-nea.org/pt/docs/iem/lasvegas04/09_Poster_Session_I/PS1_03.pdf | 4A → Li-A in molten LiCl at 920 K | Poster (excerpts) |
| R23 | T. Yang, H. Park, I. Kim, "Thermal reaction of zeolite and Cs removal in molten salt," Trans. Korean Nuclear Society Spring Meeting, Jeju, 10–11 May 2007. https://kns.org/files/pre_paper/14/596양태훈.pdf | LTA → LiAlSiO₄ in 10 min in LiCl (620–680 °C); LTA fraction survives 60 min in LiCl–KCl | Read |
| R24 | US Patent 8,535,636 (KAERI/KHNP; Eun, Cho, Park, Lee, Kim, Park), granted 17 Sept 2013. https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8535636 | "zeolite structure completely breaks down" in LiCl | Full text (excerpts) |
| R25 | R. X. Fischer, M. Šehović, W. H. Baur, C. Paulmann, T. M. Gesing, "Crystal structure and morphology of fully hydrated zeolite Na-A," *Z. Kristallogr.* **227**(7), 438–445 (2012), doi:10.1524/zkri.2012.1501 | *Fm-3c*, a = 24.550 Å | Abstract |
| R26 | N. A. Petranović, M. V. Šušić, I. Krstanović, "Structural change in the zeolite inclusion complex through high temperature treatment. Denitrated zeolite," *J. Inorg. Nucl. Chem.* **38**(2), 335–337 (1976), doi:10.1016/0022-1902(76)80423-X | NaA·NaNO₃: nitrate in the central α-cage | Abstract |
| R27 | S. H. Kim, T. N. Nguyen Thi, N. H. Heo, G. H. Kim, S. B. Hong, J. D. Head, K. Seff, "Single crystal structure of zeolite A (LTA) containing Ag₄Cl₄ nanoclusters and reduced 1,3,5-tripyrylium dimers…," *J. Phys. Chem. C* **112**(30), 11181–11193 (2008), doi:10.1021/jp801717e | Ag₄Cl₄ in β; K₃Cl²⁺ in α; Cl–K 2.70 Å; Ag–O 2.493 Å | Abstract |
| R28 | W.-T. Lim, S.-Y. Choi, B.-J. Kim, C.-M. Kim, I.-S. Lee, S.-H. Kim, N.-H. Heo, "Synthesis and crystal structure of Ag₄Br₄ nanoclusters in the sodalite cavities of fully K⁺-exchanged zeolite A (LTA)," *Bull. Korean Chem. Soc.* **26**(7), 1090–1096 (2005), doi:10.5012/bkcs.2005.26.7.1090 | Ag₄Br₄ in β; near square-planar K₄Br³⁺; Br⁻ opposite a 4-ring | Abstract |
| R29 | S. H. Kim, M. Park, Y. J. Son, H. J. Lee, G. C. Jeong, M. N. Bae, W. T. Lim, "Synthesis and crystal structure of zinc iodide in the sodalite cavities of zeolite A (LTA)," *Bull. Korean Chem. Soc.* **28**(4), 567–573 (2007), doi:10.5012/bkcs.2007.28.4.567 | Cites R34 (no longer used for numbers) | Full text (excerpts) |
| R30 | H. S. Lim, J. Y. Kim, N. H. Heo, K. Seff, "Crystal structure of zeolite LTA containing extraframework tungsten(VI) ions," *J. Phys. Chem. C* **122**(12), 6661–6668 (2018), doi:10.1021/acs.jpcc.7b12686 | Tl₇Cl₆⁺: 6 Cl in one α-cage; WCl⁵⁺ | Abstract |
| R31 | N. C. Hyatt, J. A. Hriljac, A. Choudhry, L. Malpass, G. P. Sheppard, E. R. Maddrell, "Zeolite–salt occlusion: a potential route for the immobilisation of iodine-129?" *MRS Symp. Proc.* **807** (2004), doi:10.1557/PROC-807-359 | α-cage effective diameter 11.4 Å (from the introduction) | Abstract/intro |
| R32 | A. M. Harward et al. (MS&T22 abstract lists Harward, Carlson, Yoo, Fredrickson, Patterson, Simpson), "Synthesis and characterization of super occluded LiCl-KCl in zeolite-4A as a chloride salt waste form intermediate," *J. Nucl. Mater.* (2023), PII S0022311523005834. Reported as vol. 588, art. 154816, doi:10.1016/j.jnucmat.2023.154816, with co-authors Elliott and Shaltry *(to confirm)* | Authors' statement that > 20 wt% "exceeds the reported capacity" (do not use as capacity, §1); free salt mainly NaCl; loadings up to 75 wt% | Abstract (publisher snippet) and conference abstract |
| R33 | A. Harward, L. Gardner, C. M. Decker Oldham, K. Carlson, T.-S. Yoo, G. Fredrickson, M. Patterson, M. F. Simpson, "Water sorption/desorption characteristics of eutectic LiCl-KCl salt-occluded zeolites," *J. Nucl. Fuel Cycle Waste Technol.* **20**(3), 259–268 (2022), doi:10.7733/jnfcwt.2022.021. https://inldigitallibrary.inl.gov/sites/sti/sti/Sort_61124.pdf | 20 wt% = "theoretical occlusion limit of 12 Cl⁻ per pseudo unit cell"; > 99 % occlusion; 500 °C / 24 h; 58 : 42 eutectic | Full text read |
| R34 | N. H. Heo, H. S. Kim, W. T. Lim, K. Seff, "Synthesis and crystal structure of Ag₄I₄ nanoclusters in the sodalite cavities of fully K⁺-exchanged zeolite A," *J. Phys. Chem. B* **108**(10), 3168–3173 (2004), doi:10.1021/jp031137p | Ag₄I₄ in half the β-cages; near square planar K₄I³⁺ | Abstract (OpenAlex) |

### Background (no occlusion numbers taken from these)

- R. M. Barrer, W. M. Meier, *J. Chem. Soc.* 299 (1958): early molten-salt occlusion in zeolites (cited by R11 and US 4,670,412).
- R. M. Barrer, W. M. Meier, *Trans. Faraday Soc.* **54**, 1074 (1958), doi:10.1039/TF9585401074: structural and ion-sieve properties of zeolite A.
- R. M. Barrer, A. J. Walker, *Trans. Faraday Soc.* **60**, 171 (1964), doi:10.1039/TF9646000171: salt imbibition from aqueous solution.

---

## 9. Full texts worth obtaining

| Ref | What it would add |
|---|---|
| R1, R2, R3 | Occlusion vs. composition curves |
| R7 | The only multi-halide, multi-cation-form melt dataset; would confirm the mass basis of 3.5 meq/g |
| R9 | Per-cell basis and temperature; reconciliation with R23 |
| R20 | Neutron site tables |
| R14–R16 | Two-site model parameters, including total occlusion capacity |
