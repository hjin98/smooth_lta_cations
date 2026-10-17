# Chloride-salt occlusion in zeolite A (LTA): evidence file for the occlusion model

Prepared 2026-10-16 as the literature basis for revising the occlusion model in `smooth_lta_cations`.
Every number below carries a reference tag `[Rn]` (full list in §8) and a verification grade:

| Grade | Meaning |
|---|---|
| **A** | Read directly from the paper's abstract or publisher/OSTI record |
| **B** | Read from the full text of a patent / report |
| **C** | Secondary: the number is reported by another paper or patent citing the original; original full text not checked |
| **D** | Our own derived calculation (formula given in §5 so it can be reproduced) |

Nothing here was taken from memory. Where a bibliographic detail could not be confirmed it is marked *(to confirm)*.

---

## 0. Conventions

- **Pseudo-cell** = (AlSiO₄)₁₂ = M₁₂Al₁₂Si₁₂O₄₈ for monovalent M, cubic *Pm-3m*, a ≈ 12.3 Å. Contains **one α-cage + one β-cage**. Molar mass of Na₁₂Al₁₂Si₁₂O₄₈ = 1704.6 g/mol (D).
- The true cell of Na-A is *Fm-3c*, a ≈ 24.55 Å, ≈ Na₉₆Al₉₆Si₉₆O₃₈₄ = 8 pseudo-cells [R25]. The F-centred **primitive cell = 2 pseudo-cells** (24 Al, 24 Si, 2 α + 2 β cages).
- All experimental loadings in the literature are quoted **per pseudo-cell** ("per zeolite unit cell (AlSiO₄)₁₂"). Multiply by 2 for the primitive cell.
- Each occluded Cl⁻ arrives with a charge-balancing cation, so a pseudo-cell at *n* Cl holds 12 + *n* monovalent cations (cf. Na₂₂ phase for 10 NaNO₃ [R11]).

---

## 1. Equilibrium loading data (calibration targets)

| # | System | T | Loading, Cl⁻ (or anion) per pseudo-cell | Grade | Ref |
|---|---|---|---|---|---|
| 1 | Zeolite 4A + LiCl–KCl–MeCl, Me = Na, Rb, Cs (x_MeCl 0–0.5, x_KCl/x_LiCl = 0.69) | 823 K | **11 ± 1** at infinite MeCl dilution; varies with MeCl content of zeolite | A | [R1] |
| 2 | Zeolite 4A + LiCl–KCl–MeClₙ, Me = Ca, Sr, U (x up to 0.45) | 823 K | **10 ± 1** at infinite dilution; irreversible reaction seen with UCl₃ | A | [R2] |
| 3 | ANL statement of zeolite capacity (LiCl–KCl waste salt) | process T | **12** ("12 Cl⁻ per unit cell for zeolite") | B | [R5] |
| 4 | NaA + simulated IFR salt (92.5 wt% LiCl–KCl eutectic, 5 % NaCl, 1 % CsCl, 1 % BaCl₂, 0.5 % SrCl₂), washed | 673 ± 25 K | 16.9–19.2 wt% Cl = 4.8–5.3 meq Cl/g product; cation and Cl analyses agree on 5.0 ± 0.3 meq salt/g; zeolite = 66 wt% of product (Si basis) → **≈ 12.9 (12.1–13.7)** | B (data), D (conversion) | [R6] |
| 5 | Zeolite 4A in molten **pure LiCl** | molten LiCl (T not stated in abstract) | **10–11.5 "salt molecules"** (basis presumed per pseudo-cell; *to confirm in full text*) | A | [R9] |
| 6 | Li-A + LiCl–KCl eutectic (Šušić, Petranović, Mioč 1971) | — | 3.5 meq Cl/g; only LiCl occluded, no KCl → **≈ 5.3** (per g zeolite) or **≈ 6.2** (per g product) | C (data), D (conversion) | [R7] via [R6] |
| 7 | Nitrate analogues, single-cation forms (Liquornik & Marcus) | nitrate melts (companion papers: 330–450 °C [R12, R13]) | LiNO₃ **11.7**, NaNO₃ **10**, AgNO₃ **10**, KNO₃ **≈ 0** per formula weight M₁₂[(AlO₂)₁₂(SiO₂)₁₂] | A | [R11] |
| 8 | Patent statement, zeolite A (Na, K or Li form) + molten chloride waste salt | — | "up to about 25 % by weight" salt → **≈ 10.2** Cl/cell for LiCl–KCl eutectic on NaA basis | B (data), D (conversion) | [R8] |

**Do not confuse with process loading.** The ceramic-waste-form process deliberately loads only **3.4–3.9 Cl⁻ per pseudo-cell (target 3.9)** [R4], or 3.8 [R5]. That is set by the sodalite reaction stoichiometry Na₁₂(AlSiO₄)₁₂ + 4 NaCl → 2 Na₈(AlSiO₄)₆Cl₂ [R4], not by the occlusion capacity. In the mixing study, keeping chloride below **80 % of the full salt capacity** reduced free salt [R17].

**Equilibration conditions used in the INL data set:** 500 °C, 24 h salt–zeolite contact, assumed to reach equilibrium; the fitted model has separate ion-exchange and occlusion sites and a fitted total occlusion capacity [R15].

---

## 2. Cation and salt dependence

1. **Large cations reduce or block occlusion.** KNO₃ is essentially not occluded in K-A [R11]. In K-A in nitrate melts, partial Na⁺NO₃⁻ occlusion starts only after 7 K⁺ have been exchanged for Na⁺; small amounts of K⁺ block the channels [R13]. On Li-A in LiCl–KCl, LiCl but not KCl was occluded [R7 via R6].
2. **Nitrate occlusion is in the α-cage.** Liquornik & Marcus place occluded nitrate in the α-cage; with Li/Na/Ag occlusion all cations exchange, with others exchange is limited [R12]. Petranović & Krstanović likewise describe NaA·NaNO₃ with nitrate in the central α cage [R26].
3. **Occluded vs framework cation selectivity differs.**
   - Monovalent two-site model fit to Cs–Rb–Na–Li–K data: Cs preferred in the occluded (cage) sites, Li preferred in framework sites [R14].
   - Multivalent extension: Sr²⁺ preferred in the ion-exchange site, Cs⁺ in the occlusion site; Li⁺ highest combined selectivity [R16].
4. **Capacity depends on the third component's concentration.** Lexa: occlusion extent is a function of MeCl (or MeClₙ) concentration in the zeolite; values above are the infinite-dilution limits [R1, R2]. Separation factors vs Li at infinite dilution: NaCl 0.84, RbCl 0.87, CsCl 2.31 [R1]; CaCl₂ 0.9, SrCl₂ 11.9, UCl₃ 13 [R2].

**Gap:** no single-salt NaCl-in-NaA or KCl-in-KA melt equilibrium in chloride was found. Pure NaCl/KCl melts lie above the temperature range where 4A survives (see §3), so per-salt capacities effectively exist only as components of LiCl–KCl solvent.

---

## 3. Validity regime / framework stability

- **Occlusion pathway.** LiCl–KCl + 4A: solid-state occlusion from ~450 K to ~625 K, abrupt molten-salt occlusion at 628 K (eutectic melting), accompanied by lattice contraction [R18]. UCl₃–LiCl–KCl: LiCl/KCl begins diffusing into the lattice at 450 K [R19].
- **Temperature.** Mixing studies: occlusion improved above 775 K [R17]. Lexa's equilibria were at 823 K [R1, R2].
- **4A → sodalite.** Salt-loaded 4A tends to transform to sodalite; Na/Ca 5A retains LTA more readily [R20].
- **Molten LiCl.** Short-contact occlusion of 10–11.5 per cell is reported [R9]. However:
  - 4A changed to Li-A when mixed with molten LiCl at 920 K [R22];
  - in molten LiCl at 620–680 °C, 4A converted to eucryptite (LiAlSiO₄), whereas in LiCl–KCl it behaved differently [R23];
  - a later patent states that in LiCl the zeolite structure completely breaks down [R24].
  Treat pure-LiCl loadings as metastable/short-time values.

---

## 4. Structural evidence on Cl⁻ siting

Diffraction on genuinely melt-equilibrated LTA is limited.

- **Neutron diffraction on salt-occluded 4A/5A (Li, K, Cs, Sr, Ba, Y chlorides) [R20].** Cations act both as charge compensators and as members of occluded clusters; clusters roughly fill the 6–8 Å pores. Abstract only — site tables not checked.

Single-crystal X-ray studies (vapour or aqueous routes, **not melt equilibrium**) give halide motifs inside LTA:

| Structure | β-cage content | α-cage content | Total halide / pseudo-cell (D, from formula) | Ref |
|---|---|---|---|---|
| \|K₂.₃₅Ag₁.₁(Ag₄Cl₄)₀.₄₅(C₃H₃O₃)₁.₁(K₃Cl)₃.₄₅(K₃(OH)₂)₀.₅₅(H₂O)≥₃\|[Si₁₂Al₁₂O₄₈]-LTA | Ag₄Cl₄ in ~45 % of sodalite cavities; each Cl bonds a K⁺ through a 6-ring (Cl–K 2.70 Å) | K₃Cl²⁺ units predominate | 0.45×4 + 3.45 = **5.25 Cl** | [R27] |
| K₉(K₄Br)·0.75 Ag₄Br₄-LTA | Ag₄Br₄ (T_d, ~7.9 Å) in ~75 % of sodalite units | K₄Br³⁺ | 3 + 1 = **4 Br** | [R28] |
| \|K₉(K₄I)(Ag₄I₄)₀.₅\|-LTA (as cited in [R29]) | Ag₄I₄ in half of sodalite cavities | K₄I³⁺ | 2 + 1 = **3 I** | C [R29] |
| \|(WCl)₀.₆(Tl₇Cl₆)₀.₅Tl₈.₅\|-LTA | — | half of α-cages: Tl₇Cl₆⁺ — central Tl⁺ octahedrally coordinated by 6 Cl⁻, each Cl also bonded to an 8-ring Tl⁺; other half: WCl⁵⁺ opposite 6-rings | 0.6 + 3 = **3.6 Cl** (**6 Cl in one α-cage**) | [R30] |

Take-aways for siting:

- In the α-cage, Cl⁻ is anchored to cations at 6-ring sites (on the 3-fold axes), to cations in the 8-ring windows (on the 4-fold axes), or to a cage-centre cation. It does not appear as isolated anions.
- The β-cage can hold **more than one** halide when Ag⁺ is present (M₄X₄ clusters). For Na/K/Li the sodalite-like [M₄Cl]³⁺ unit (1 Cl per β) remains the natural default, consistent with the sodalite reaction stoichiometry [R4].

---

## 5. Derived numbers (Grade D) — reproducible

Atomic masses: Na 22.990, Li 6.941, K 39.098, Al 26.982, Si 28.086, O 15.999, Cl 35.453, Cs 132.905.

| Quantity | Formula | Result |
|---|---|---|
| M(Na₁₂Al₁₂Si₁₂O₄₈) | — | 1704.6 g/mol (Si 19.77 wt% — matches 19.8 % in [R6]) |
| M(Li₁₂Al₁₂Si₁₂O₄₈) | — | 1512.1 g/mol |
| [R6] Cl per cell | q / (0.66 / 1704.6), q = 4.7, 5.0, 5.3 meq/g | 12.1, **12.9**, 13.7 |
| [R7] Cl per cell, per g zeolite | 3.5e-3 × 1512.1 | 5.3 |
| [R7] Cl per cell, per g product | n = 3.5e-3·M_LiA / (1 − 3.5e-3·M_LiCl) | 6.2 |
| LiCl–KCl eutectic mean M (58.5 : 41.5 mol) | 0.585·42.39 + 0.415·74.55 | 55.7 g/mol |
| wt% salt at 10 / 11 / 12 Cl per cell (NaA basis) | 55.7n / (1704.6 + 55.7n) | 24.6 / 26.5 / 28.2 % |
| Cl per cell at 25 wt% [R8] | w·M_NaA / ((1−w)·55.7) | 10.2 |
| α-cage volume, 11.4 Å free diameter [R31] | (4/3)π(5.7 Å)³ | ≈ 776 Å³ |
| Ion pairs per α-cage at bulk-melt packing | V_α / (M / ρ_melt / N_A); ρ ≈ 1.50 (LiCl), 1.55 (NaCl), 1.52 (KCl), 2.79 (CsCl) g/cm³ | LiCl 16.5, NaCl 12.4, KCl 9.5, CsCl 7.7 |

The melt densities are approximate handbook values near the melting point; they were **not** retrieved for this file. The packing estimate is an order-of-magnitude size argument only. The observed ~10 Cl/cell for LiCl–KCl is ~70–80 % of bulk-melt packing in the α-cage, and the ordering matches the nitrate trend (Li > Na ≫ K).

Geometry check for candidate α-cage Cl positions at radius r = 3.5 Å from the cage centre:

- 8 sites on the 3-fold axes: nearest Cl–Cl 4.04 Å.
- 6 sites on the 4-fold axes: nearest Cl–Cl 4.95 Å.
- **Between** a 3-fold and a 4-fold site: 3.22 Å — too short for simultaneous occupation.

---

## 6. Implications for the model (recommendations, not evidence)

1. **β-cage:** keep one Cl site at the cage centre as default for Li/Na/K (sodalite-like [M₄Cl]³⁺). Make the capacity a per-cation parameter so Ag-type M₄X₄ clusters (up to 4 X per β) can be represented [R27, R28].
2. **α-cage:** replace the single cage-centre site with a multi-site set. Candidates:
   - centre (1);
   - 3-fold-axis sites inward of 6-ring cations (8);
   - 4-fold-axis sites near 8-ring windows (6);
   - sites facing 4-rings (12).
   
   Enforce a Cl–Cl exclusion distance, or use off-lattice/grand-canonical sampling, so the capacity emerges rather than being hard-coded.
3. **Charge balance:** insert or remove Cl⁻ together with a cation; track occluded-site cation identity separately from framework sites [R14, R16].
4. **Calibration targets at ~773–823 K:**
   - total ≈ 11 ± 1 Cl per pseudo-cell for alkali-rich LiCl–KCl [R1];
   - total ≈ 10 ± 1 with Ca/Sr/U [R2];
   - upper bound ≈ 12–13 [R5, R6].
   
   Use the separation factors in §2 and the INL two-site fits [R14–R16] to constrain the cation partition.
5. **Cation-size penalty:** K⁺ (and larger) occupancy of 8-ring windows should suppress α-cage loading (KNO₃ ≈ 0 in K-A [R11]; blocking by a few K⁺ [R13]).

---

## 7. Corrections relative to the earlier chat summary

- **Pure LiCl:** the earlier summary said the framework "breaks down" in molten LiCl. Seo et al. 2006 [R9] actually report 10–11.5 salt molecules occluded in molten LiCl. The breakdown and transformation reports [R22–R24] are time- and temperature-dependent. Both are now listed (§3).
- **Gougar et al. two-site model [R14]:** the author list is not fully confirmed. A patent cites it as Simpson & Gougar; OSTI confirms journal, volume, issue and pages.
- **K₄I³⁺ geometry:** the earlier "near square planar" description came from a search-snippet of a citing work. Only the formula is retained here (Grade C).
- **Unchanged:** all other numbers (11 ± 1, 10 ± 1, 12, ≈ 13, 5–6, nitrate 11.7 / 10 / 10 / ≈ 0, 3.4–3.9) were re-checked against the sources and stand.

---

## 8. References

Status column: what was checked (abstract/record vs full text).

| Tag | Reference | Supports | Status |
|---|---|---|---|
| R1 | D. Lexa, I. Johnson, "Occlusion and ion exchange in the molten (lithium chloride–potassium chloride–alkali metal chloride) salt + zeolite 4A system with alkali metal chlorides of sodium, rubidium, and cesium," *Metall. Mater. Trans. B* (2001), doi:10.1007/s11663-001-0028-4 *(volume/pages to confirm)*. OSTI 943218: https://www.osti.gov/biblio/943218 | 11 ± 1 Cl/cell, 823 K; separation factors | Abstract |
| R2 | D. Lexa, "Occlusion and ion exchange in the molten (LiCl + KCl + alkaline-earth chloride) salt + zeolite 4A system … and … with UCl₃," *Metall. Mater. Trans. B* **34**, 201–208 (2003), doi:10.1007/s11663-003-0007-z. OSTI 949686: https://www.osti.gov/biblio/949686 | 10 ± 1 Cl/cell; separation factors; UCl₃ reaction | Abstract |
| R3 | D. Lexa, *J. Chem. Thermodyn.* (1999) *(details to confirm; cited in R1)* | Binary LiCl–KCl occlusion thermodynamics | Not checked |
| R4 | US Patent 10,600,527, "Ceramic waste form production from used nuclear fuel." https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10600527 | Sodalite stoichiometry; 3.4–3.9 Cl/cell process loading | Full text (excerpts) |
| R5 | ANL-West, "Ceramic waste form production and development at ANL-West." https://corpora.tika.apache.org/base/docs/govdocs1/230/230873.pdf | 12 Cl/cell capacity vs 3.8 process loading | Full text (excerpts) |
| R6 | M. A. Lewis, T. R. Johnson, US Statutory Invention Registration H1227 (1993), "Method for immobilizing mixed waste chloride salts containing radionuclides and other hazardous wastes." https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/H1227 | NaA + IFR salt Cl analysis; prior-art summary of Šušić and Liquornik | Full text |
| R7 | M. V. Šušić, N. A. Petranović, D. A. Mioč, *J. Inorg. Nucl. Chem.* **33**, 2667 (1971) | Li-A + LiCl–KCl: 3.5 meq Cl/g, LiCl only | Secondary via R6 and R21 |
| R8 | M. A. Lewis, C. Pereira, US Patent 5,613,240 (1997), "Method of preparing sodalite from chloride salt occluded zeolite." https://patents.justia.com/patent/5613240 | "up to about 25 % by weight" occluded salt | Full text (excerpts) |
| R9 | C. S. Seo, B. H. Park, S. B. Park, K.-J. Jung, S. W. Park, S. H. Kim, "Study on the characteristics of the ion exchange of zeolite 4A in a molten LiCl system," *J. Chem. Eng. Japan* **39**(1), 27–33 (2006), doi:10.1252/jcej.39.27 | 4A in molten LiCl occludes 10–11.5 salt molecules | Abstract |
| R10 | M. A. Lewis, D. F. Fischer, … Smith *(full author list to confirm)*, "Salt-occluded zeolites as an immobilization matrix for chloride waste salt," *J. Am. Ceram. Soc.* **76**(11), 2826–2832 (1993) | Foundational ANL paper (capacity, leach) | Record only |
| R11 | M. Liquornik, Y. Marcus, "Ion exchange in molten salts. II. The occlusion of lithium, sodium, potassium and silver nitrates in the respective forms of zeolite A," *Isr. J. Chem.* **6**, 115–121 (1968), doi:10.1002/ijch.196800017 | 11.7 / 10 / 10 / ≈ 0; Na₂₂ model | Abstract |
| R12 | M. Liquornik, Y. Marcus, "Ion exchange in molten salts. I. …sodium zeolite A in molten NaNO₃…," *J. Phys. Chem.* **72**, 2885–2889 (1968) | Nitrate occlusion in α-cage; exchange limited without occlusion | Abstract |
| R13 | M. Liquornik, Y. Marcus, "Ion exchange in molten salts. V. Potassium zeolite A as an ion exchanger in nitrate melts," *J. Phys. Chem.* **75**, 2523–2525 (1971) | K⁺ blocks channels; occlusion only after 7 K → Na | Abstract |
| R14 | M. L. D. Gougar, M. F. Simpson et al. *(author list to confirm)*, "Two-site equilibrium model for ion exchange between monovalent cations and zeolite-A in a molten salt," *Ind. Eng. Chem. Res.* **42**(18), 4208–4212 (2003). OSTI 961100: https://www.osti.gov/biblio/961100 | Cs prefers cage (occluded) sites; Li prefers framework | Abstract |
| R15 | T.-S. Yoo, S. M. Frank, M. F. Simpson, P. A. Hahn, T. J. Battisti, S. Phongikaroon, "Salt-zeolite ion-exchange equilibrium studies for a complete set of fission products in molten LiCl-KCl," *Nucl. Technol.* **171**(3), 306–315 (2010), doi:10.13182/NT10-A10866 | 500 °C / 24 h equilibrium; fitted total occlusion capacity | Abstract |
| R16 | S. Phongikaroon, M. F. Simpson, "Equilibrium model for ion exchange between multivalent cations and zeolite-A in a molten salt," *AIChE J.* **52**(5), 1736–1743 (2006); conf. version INL/CON-05-00460 | Sr prefers IX site; Cs prefers occlusion site | Abstract |
| R17 | "Mixing of zeolite powders and molten salt" (ANL), abstract via https://worldwidescience.org/topicpages/s/salt+occluded+zeolite.html *(full citation to confirm)* | > 775 K improves occlusion; < 80 % of capacity reduces free salt | Abstract |
| R18 | D. Lexa et al. *(co-authors to confirm)*, "Occlusion of the LiCl–KCl eutectic salt by zeolite 4A," APS Activity Report 2001. https://www-stg.aps.anl.gov/sites/www.aps.anl.gov/files/APS-sync/activity_reports/apsar2001/LEXAD1.PDF | Solid-state then molten occlusion at 628 K; lattice contraction | Report (excerpts) |
| R19 | D. Lexa et al. *(co-authors to confirm)*, "The interaction between UCl₃ and zeolite 4A by in situ temperature-resolved powder XRD," APS report 1999. https://www.aps.anl.gov/sites/www.aps.anl.gov/files/APS-Uploads/APS-Science/1999/lexad1.pdf | Salt diffusion into lattice from 450 K | Report (excerpts) |
| R20 | J. W. Richardson Jr., "Salt-occluded zeolite waste forms: crystal structures and transformability," *MRS Proc.* **465**, 395–400 (1997), doi:10.1557/PROC-465-395. OSTI 541788 | Neutron diffraction: cation clusters fill pores; 4A → sodalite, 5A resists | Abstract |
| R21 | M. A. Lewis, D. F. Fischer, J. J. Laidler, "Leach resistance properties and release processes for salt-occluded zeolite A," *MRS Proc.* **294** (1993), doi:10.1557/PROC-294-95 | Cites R7 | Record |
| R22 | J.-G. Kim, J.-H. Lee, J.-H. Yoo, J.-H. Kim, "Immobilisation of molten LiCl waste using zeolite A," OECD-NEA IEM 2004 poster PS1-03. https://test.oecd-nea.org/pt/docs/iem/lasvegas04/09_Poster_Session_I/PS1_03.pdf | 4A → Li-A in molten LiCl at 920 K | Poster (excerpts) |
| R23 | T. Yang, H. Park, I. Kim, "Thermal reaction of zeolite and Cs removal in molten salt," KNS meeting paper. https://kns.org/files/pre_paper/14/596양태훈.pdf | 4A → eucryptite in molten LiCl (620–680 °C) | Paper (excerpts) |
| R24 | US Patent 8,535,636. https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8535636 | "zeolite structure completely breaks down" in LiCl; cites R2, R14, R16 | Full text (excerpts) |
| R25 | "Crystal structure and morphology of fully hydrated zeolite Na-A" (*Fm-3c*, a = 24.550 Å, Na₉₄Al₉₄Si₉₈O₃₈₄). https://www.researchgate.net/publication/258724844 *(journal details to confirm)* | True-cell / pseudo-cell relation | Abstract |
| R26 | N. A. Petranović, M. V. Šušić, I. Krstanović, "Structural change in the zeolite inclusion complex through high temperature treatment. Denitrated zeolite," *J. Inorg. Nucl. Chem.* **38** (1976). https://cer.ihtm.bg.ac.rs/handle/123456789/4406 | NaA·NaNO₃: nitrate in α-cage | Abstract |
| R27 | S. H. Kim, T. N. Nguyen Thi, N. H. Heo, G. H. Kim, S. B. Hong, J. D. Head, K. Seff, "Single crystal structure of zeolite A (LTA) containing Ag₄Cl₄ nanoclusters and reduced 1,3,5-tripyrylium dimers…," *J. Phys. Chem. C* **112**(30), 11181–11193 (2008), doi:10.1021/jp801717e | Ag₄Cl₄ in β; K₃Cl²⁺ in α | Abstract |
| R28 | "Synthesis and crystal structure of Ag₄Br₄ nanoclusters in the sodalite cavities of fully K⁺-exchanged zeolite A (LTA)," *Bull. Korean Chem. Soc.* **26**(7), 1090–1096 (2005), doi:10.5012/bkcs.2005.26.7.1090 | Ag₄Br₄ in β; K₄Br³⁺ in α | Abstract |
| R29 | "Synthesis and crystal structure of zinc iodide in the sodalite cavities of zeolite A (LTA)," *Bull. Korean Chem. Soc.* **28**(4), 567 (2007), doi:10.5012/bkcs.2007.28.4.567 | Cites earlier \|K₉(K₄I)(Ag₄I₄)₀.₅\|-LTA | Full text (excerpts) |
| R30 | H. S. Lim, J. Y. Kim, N. H. Heo, K. Seff, "Crystal structure of zeolite LTA containing extraframework tungsten(VI) ions," *J. Phys. Chem. C* **122**(12), 6661 (2018) | Tl₇Cl₆⁺: 6 Cl in one α-cage | Abstract |
| R31 | "Zeolite–salt occlusion: a potential route for the immobilisation of iodine-129?" https://www.researchgate.net/publication/240638899 *(journal details to confirm)* | α-cage effective diameter 11.4 Å | Abstract |
| R32 | A. M. Harward, K. Carlson, T.-S. Yoo, G. Fredrickson, M. Patterson, M. Simpson (MS&T22 abstract); journal version "Synthesis and characterization of super occluded LiCl-KCl in zeolite-4A as a chloride salt waste form intermediate," *J. Nucl. Mater.* (2023/24) — https://www.sciencedirect.com/science/article/abs/pii/S0022311523005834 *(journal author list/volume to confirm)* | Loadings > 20 wt% exceed capacity; excess salt free, mainly NaCl | Abstract |

### Background (no occlusion numbers taken from these)

- R. M. Barrer, W. M. Meier, *J. Chem. Soc.* 299 (1958) — early molten-salt occlusion in zeolites (cited by R11, US 4,670,412).
- R. M. Barrer, W. M. Meier, *Trans. Faraday Soc.* **54**, 1074 (1958), doi:10.1039/TF9585401074 — structural and ion-sieve properties of zeolite A.
- R. M. Barrer, A. J. Walker, *Trans. Faraday Soc.* **60**, 171 (1964), doi:10.1039/TF9646000171 — salt imbibition from aqueous solution.

---

## 9. Full texts worth obtaining

These would replace Grade C entries and give site tables and composition curves:

- R1, R2, R3 — occlusion vs. composition curves.
- R7 — the only multi-halide, multi-cation-form melt dataset.
- R9 — to confirm the per-cell basis.
- R20 — neutron site tables.
- R14–R16 — two-site model parameters, including total occlusion capacity.
