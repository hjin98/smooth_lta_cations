# Lattice (volume) contraction of zeolite A on chloride-salt loading: evidence file

**Revision 3.0, 2026-10-17 (new file).** Companion to `lta_salt_occlusion_evidence.md` (references R1–R34) and `molten_salt_pair_distances.md` (P1–P26). New references are numbered **Q1, Q2, …** and are shared with `lta_coordination_numbers.md`. §9 lists what the independent reviews changed.

**Bottom line.** No source readable in this session gives a numeric lattice parameter for zeolite A loaded with molten chloride salt. The only direct evidence is qualitative: 4A contracts on LiCl–KCl occlusion, first gradually from about 450 K and then by growth of a second LTA phase with a "significantly smaller lattice parameter" from about 650 K [R18, R19]. **No matched dry/loaded pair at the same temperature exists**, so f_V for 4A + LiCl–KCl is **unavailable**. §5 gives a bounding estimate (grade E).

---

## 0. Conventions and grades

Same conventions as `lta_salt_occlusion_evidence.md` §0:

- **Pseudo-cell** (AlSiO₄)₁₂ = M₁₂Al₁₂Si₁₂O₄₈, *Pm-3m*, a ≈ 12.3 Å (1 α-cage + 1 β-cage).
- **True cell** of Na-A, *Fm-3c*, a ≈ 24.55 Å = 8 pseudo-cells. a_true ≈ 2 × a_pseudo.
- **Every a below names its cell.** Ratios f_a are only formed between values in the same cell setting, or after converting with a_true = 2 a_pseudo (exact for the metric; the Si/Al-ordering superstructure does not change it).
- **Contraction factors:** f_a = a_loaded / a_dry, f_V = f_a³.

| Grade | Meaning |
|---|---|
| **A** | Read from an abstract, database record (COD/IZA) or search-result snippet |
| **B** | Read from the full text (paper, report, ESI or accepted manuscript) |
| **C** | Secondary: reported by another source citing the original |
| **D** | Own calculation; formula given |
| **E** | Estimate/bound, not a measurement; reasoning given |

"Read" means read in this session (2026-10-17). Values relayed by a research sub-agent are marked †; for those the agent's quote and URL are in the session notes and the number was re-checked by the independent review (§9).

---

## 1. Reference cells: dehydrated ("dry") and hydrated zeolite A

### 1a. Na-A (zeolite 4A)

| # | State | T of measurement | a (Å) | Cell | Method | Grade | Ref |
|---|---|---|---|---|---|---|---|
| 1 | Hydrated, 54 % RH (NIST RM 8851, Na₉₇.₃Al₉₆.₂Si₉₅.₅O₃₈₄·209.4H₂O) | ambient (not stated) | **24.589 ± 0.002** (1 sd, 13 specimens) | Fm-3c | synchrotron XRD, Le Bail | B | [Q1, Table 8] |
| 2 | Dehydrated (200 °C, 24 h; sealed under He), same material | ambient (not stated) | **24.560 ± 0.006** (1 sd, 11 specimens) | Fm-3c | neutron powder, Rietveld | B | [Q1, Table 9] |
| 3 | Dehydrated | not stated | **24.555(2)** | Fm-3c | single crystal | C via [Q1, Table 17]; A via IZA-SC CIF (24.555) | [Q4] |
| 4 | Hydrated | room temperature | 24.550(3) | Fm-3c | single crystal (synchrotron) | A | [R25] |
| 5 | Hydrated | not stated | 24.61(1) (as 2.461(1) nm) | Fm-3c | — | C via [Q1, Table 17] | Gramlich & Meier 1971 |
| 6 | Dehydrated Na₁₂-A | not stated | **12.263** (no esd on record) | Pm-3m (by size; SG not on the record) | single crystal | A† (COD record) | [Q6] |

Notes:

- [Q1] prints the space group as "Fd-3c". That symbol does not exist for this structure; it is read as a typo for *Fm-3c*.
- RM 8851 contains a small sodalite impurity (a = 9.976 Å, *P-43n*): "not quantified by XRD but presumed small", ≤ ~2 % from neutron data to ~3 % from composition [Q1].
- Row 6 doubled is 24.526 Å, 0.14 % below row 2 and 0.12 % below row 3. The records give no temperature, so this spread is treated as a scatter between crystals/methods, not as a physical effect.
- **Dehydration barely changes a for Na-A:** f_a(hydrated/dry) = 24.589 / 24.560 = 1.0012 (D, rows 1–2, same material). Water loss during in-situ heating therefore cannot masquerade as the occlusion contraction in [R18].

### 1b. Other cation forms (dry references needed for the single-crystal analogues in §3)

| # | Form (dehydrated unless stated) | T | a (Å) | Cell | Grade | Ref |
|---|---|---|---|---|---|---|
| 7 | K₁₁.₁-A | not stated | **12.317** | Pm-3m | A† (COD 1530597) | [Q5] |
| 8 | K₁₁-A (COD formula field gives O42 vs. calculated O48) | not stated | 12.309 | Pm-3m | A† (COD 1536625, read via mirror) | [Q5] |
| 9 | K₉₁ per Fm-3c cell (Si,Al-ordered refinement) | not stated | 24.6 (truncated on record) | Fm-3c | A† | [Q5] |
| 10 | Ag₅.₆K₆.₄-A | 21 °C | 12.255(1) | Pm-3m | A | [Q9] |
| 11 | Tl₁₂-A (anhydrous) | 294 K | **12.153(4)** | Pm-3m | A | [Q7] |
| 12 | Ag (dehydrated Ag-A) | 300 K | **24.6636(2)** | Fm-3c | B | [Q2] |
| 13 | Ag (hydrated Ag-A) | 300 K | 24.5970(2) | Fm-3c | B | [Q2] |
| 14 | Ca₆-A | not stated | 12.278 | Pm-3m | A† | Firor & Seff 1978 (COD 4124433) |
| 15 | Sr-A (Al-deficient, "presence of Al complex") | not stated | 24.68 | Fm-3c | A† | Pluth & Smith 1982 (COD 4124454) |
| 16–18 | Ba₁Na₁₀-A / Ba₃.₅Na₅-A / partly dehydrated Ba₆-A | not stated | 12.262(3) / 12.267(7) / 12.189(2) | Pm-3m | A† | [Q26] (defined in the CN file) |
| 19 | Cs₇.₃Ag₄.₇-A | not stated | 12.282(1) | Pm-3m | A | Kim & Seff 1984 |
| 20 | Pure-silica LTA (ITQ-29) | 300 K | 11.853(2) | Pm-3m | B | [Q2] |

Which of the K-A COD records (rows 7–9) corresponds to which crystal or refinement, and at what temperature, could not be determined without the paper. Rows 7 and 8 are both K₁₁-type records. Row 7 is used below as the dry K-A reference; row 8 gives the sensitivity (0.07 %).

---

## 2. Salt-loaded zeolite A: what the sources say

### 2a. Molten LiCl–KCl and 4A (the system to be modelled)

| Source | Sample | What is reported | Grade |
|---|---|---|---|
| [R18] Lexa, Leibowitz & Kropf, APS 2001 | 4A (UOP, < 0.2 wt % H₂O) + 30 wt % LiCl–KCl eutectic (salt-mass basis: 0.30 × 1704.6 / (0.70 × 55.5) = **13.2 Cl/pseudo-cell** if all occluded, i.e. above the 10.5 ± 1 capacity; D, formula of `lta_salt_occlusion_evidence.md` §5); 300 → 800 K at 10 K/min; λ = 1.1273 Å; Rietveld (GSAS) | 300 to ~450 K: zeolite peaks shift to lower 2θ (thermal expansion). ~450 to 800 K: shift reversed, "implying a zeolite lattice contraction". ~625 K: occlusion exotherm (eutectic m.p. 628 K). From ~650 K a second set of peaks grows "at the high-2θ side of every zeolite peak": "zeolite with a significantly smaller lattice parameter", growing "at the expense of the original zeolite". The original phase's a "continues to decrease but never reaches the lattice parameter of the new zeolite phase". **No numeric a in the text; values are only in Figs. 3–4.** | B |
| [R19] Lexa, Leibowitz & Kropf, APS (1999 file) | 4A + 30.14 wt % UCl₃–LiCl–KCl; same set-up | Below 450 K expansion; above 450 K "the zeolite lattice contraction becomes apparent"; LiCl/KCl diffusion into the lattice from 450 K causes "a significant zeolite 4A lattice contraction". No 4A decomposition products (Si, SiO₂, Al₂O₃, NaCl). **Numbers only in Fig. 3.** | B |
| [R20] Richardson, MRS 465 (neutron Rietveld, Li/K/Cs/Sr/Ba/Y chlorides) | salt-occluded 4A and 5A | 4A transforms to sodalite; 5A "more readily" retains LTA. Abstract only; lattice parameters not reached. | A |
| [Q10] Pereira & Babcock, ANL 1996 | fully loaded zeolite A | "approximately 11 to 13" Cl per unit cell (loading only; no a) | B† |

Interpretation [R18, authors' speculation]: solid-state/vapour occlusion between ~450 and 625 K gives partially loaded particles; molten occlusion above 625 K gives particles with a highly occluded crust (new small-a phase) and a partially occluded core (original phase, continuously contracting). **The contraction is therefore not uniform: a two-phase (or graded) state forms, and the fully loaded phase is the one with the smaller a.**

### 2b. Single-crystal halide-containing LTA (analogues, not melt-equilibrated)

These come from solution or vapour routes (see `lta_salt_occlusion_evidence.md` §4).

| # | Structure | Halide / pseudo-cell (D) | T | a (Å) | Cell | Grade | Ref |
|---|---|---|---|---|---|---|---|
| 21 | K₉(K₄Br)·0.75Ag₄Br₄-LTA (from Ag₁₂-A + KBr/CH₃OH) | 4 Br | 294 K | **12.186(1)** | Pm-3m | A | [R28] |
| 22 | \|K₂.₃₅Ag₁.₁(Ag₄Cl₄)₀.₄₅(C₃H₃O₃)₁.₁(K₃Cl)₃.₄₅(K₃(OH)₂)₀.₅₅(H₂O)≥₃\|-LTA | 5.25 Cl | 294 K | **12.292(1)** | Pm-3m | A | [R27] |
| 23 | \|Ta²⁺₀.₁₂₃(TaCl₆Tl₃)²⁺₀.₂₄₆Tl⁺₉.₈₄\|[Si₁₃.₄Al₁₀.₆O₄₈]-LTA (from Tl₁₂-A + TaCl₅ vapour, 473 K) | 1.5 Cl | not stated | **12.099(1)** | Pm-3m | A | [Q8] |
| 24 | \|(WCl)₀.₆(Tl₇Cl₆)₀.₅Tl₈.₅\|-LTA (from Tl₁₂-A + WCl₆ vapour, 453 K, 48 h) | 3.6 Cl | not stated | **12.061(1)** | Pm-3m | A (verified by review 3) | [R30] |

Caveats that make these poor proxies for 4A + LiCl–KCl: row 22 is hydrated and contains organics; rows 21–22 contain Ag⁺ clusters; row 23's framework is **dealuminated** (Si/Al = 13.4/10.6), which by itself shortens a; rows 23–24 are Tl/Ta/W systems.

---

## 3. Contraction factors from available pairs (grade D unless stated)

f_a = a_loaded / a_dry, f_V = f_a³. Both members of each pair are in the pseudo-cell setting. Rows 21 and 11 were measured at 294 K; **rows 7, 8, 23 and 24 give no temperature** (ambient is assumed), so every f_a here carries "(T of one member not stated)". The hydration state of row 21 is not stated in the [R28] abstract.

| Pair (loaded / dry) | Rows | f_a | f_V | Usable as a 4A + LiCl–KCl proxy? |
|---|---|---|---|---|
| K₉(K₄Br)·0.75Ag₄Br₄-A / dry K₁₁.₁-A | 21 / 7 | **0.9894** | **0.968** | Weak: Br, Ag₄Br₄ in β, K form, methanol route |
| same / dry K₁₁-A (second record) | 21 / 8 | 0.9900 | 0.970 | sensitivity to the K-A reference |
| TaCl₆Tl₃-LTA / dry Tl₁₂-A | 23 / 11 | 0.9956 | 0.987 | No: framework dealuminated, Tl form |
| (WCl)(Tl₇Cl₆)-LTA / dry Tl₁₂-A | 24 / 11 | 0.9924 | 0.977 | No: W/Tl system, vapour route |
| Ag-A hydrated / dehydrated (300 K, same material) | 13 / 12 | 0.9973 | 0.992 | Shows guest filling can **contract** LTA (cf. Na-A, where water expands it by +0.12 %) |

**No pair matches the cation form, the salt and the temperature of 4A + LiCl–KCl.** The row 22 structure [R27] has no dry reference of the same mixed K/Ag/organic composition.

---

## 4. Temperature: separating thermal expansion from occlusion

| Quantity | Value | T range | Grade | Ref |
|---|---|---|---|---|
| α_V, dehydrated Na-A | **−6.34 × 10⁻⁶ K⁻¹** (negative thermal expansion) | 100–300 K | B† | [Q3] |
| α_V, dehydrated K-A | +2.03 × 10⁻⁶ K⁻¹ | 100–300 K | B† | [Q3] |
| α_V, dehydrated Li-A | +0.79 × 10⁻⁶ K⁻¹ | 100–300 K | B† | [Q3] |
| α_V, Rb₀.₇₉Na₀.₂₁-A / Cs₀.₅₈Na₀.₄₂-A | −0.78 / −1.76 × 10⁻⁶ K⁻¹ | 100–300 K | B† | [Q3] |
| α_a, dehydrated Ag-A | −3.3 × 10⁻⁶ K⁻¹ (D: (24.6636 − 24.6801)/24.6801/200 K) | 100–300 K | D from B | [Q2] |
| α_a, ITQ-29 | −8.0 × 10⁻⁶ K⁻¹ (D, same formula) | 100–300 K | D from B | [Q2] |
| Dry 4A + salt mixture, in situ | expands from 300 to ~450 K (sign only) | 300–450 K | B | [R18, R19] |

- **No α for dehydrated 4A above 300 K, and none for salt-loaded 4A, was found.**
- Size of the thermal term (E, extrapolated outside the measured 100–300 K range): if the Na-A value held from 300 to 800 K, ΔV/V = −6.34 × 10⁻⁶ × 500 = −0.3 %. Even the sign is unsupported above 300 K: [R18, R19] saw the 4A + salt mixture **expand** from 300 to ~450 K. It is a size estimate, not a bound. That is an order of magnitude smaller than the bound in §5, so thermal expansion does not limit how [R18] is interpreted, but **a 300 K dry cell must not be divided into an 800 K loaded cell when |1 − f_V| < 1 %**.
- [R18] measured the original and the new phase in the same pattern at the same temperature (650–800 K), **but the original phase is not dry**: above ~450 K it is itself partly loaded and contracting ("continues to decrease but never reaches the lattice parameter of the new zeolite phase"). Because a_orig(T) < a_dry(T), the ratio a_new/a_orig at one T is only an **upper bound on f_a** (a lower bound on the contraction): f_a ≤ a_new/a_orig. A loaded/dry ratio needs the dry 4A cell at the same T. It could be approximated from the original-phase peak at its ~450 K maximum plus a thermal correction (≈ 0.3 % in V, sign uncertain; E). The figure data were not available.

---

## 5. Recommended f_V for 4A + LiCl–KCl at ≈ 10–11 Cl per pseudo-cell

**Status: unavailable (no matched dry/loaded pair at the same temperature).**

### 5a. Bound from [R18] (grade E, derived D)

[R18] Fig. 3 is a "detail of a temperature-resolved XRD pattern" whose 2θ axis spans **17.3–17.8°** (λ = 1.1273 Å). It shows the original zeolite peak and the new high-2θ peak of the contracted phase. The window comes from the **axis labels only**; the text says nothing about it. Assuming every relevant peak position (the original peak at its ~450 K maximum and the new peak up to 800 K) lies inside the window, the ratio of any two d-spacings in it obeys

f_a ≥ sin(17.3°/2) / sin(17.8°/2) = **0.972**, i.e. **f_V ≥ 0.919** (D).

Because the ~450 K original peak is the closest available stand-in for the dry cell, this bounds a_new(≤ 800 K)/a_orig(~450 K), i.e. it mixes temperatures; the thermal term (≈ 0.3 % in V, sign uncertain, §4) is small against the 8 % bound.

Upper bound: f_a < 1, because the new phase is stated to have a smaller lattice parameter at the same temperature [R18].

### 5b. Working value for the model (grade E)

| Quantity | Value | Basis |
|---|---|---|
| f_a | **0.99 ± 0.01** | E: the K/Br analogue (rows 21/7 and 21/8) gives 0.989–0.990 (§3). The rejected Tl pairs (0.992–0.996) are context only, but fall inside the same ±0.01. All lie inside the [R18] bound |
| f_V | **0.97 ± 0.03** | E: f_a³; the ± covers 0.94–1.00 |
| Hard bounds | 0.919 ≤ f_V < 1 | §5a |
| Temperature | Applies at 650–800 K (in-situ, [R18]) and, by assumption, after quench to room temperature | not separated by data |
| Loading | Applies to the fully loaded ("crust") phase. [R18] used 30 wt % salt, above the ≈ 10.5 ± 1 Cl capacity, so the new phase is presumably at or near capacity (inference) | — |

**Use this only as a provisional input.** For a calibrated value, either (i) digitise [R18] Fig. 3/4 or obtain its Rietveld output (taking the dry reference from the ~450 K maximum, per the caveat in §4) (Lexa et al., *Thermochim. Acta*, "submitted" per [R18]; not found), or (ii) compute it: relax the dry and the salt-loaded Na-A pseudo-cell with the same functional (DFT or the fine-tuned MLFF) at 0 K and at 800 K (NPT), at 0, 4, 8 and 11 Cl/cell, and report f_V(n) directly. Route (ii) is within the user's toolset and would also give the loading dependence that no source provides.

### 5c. Dependence on loading, salt and cation

- **Loading:** no data. [R18] implies a continuous contraction of the partly loaded core and a lower limit reached by the fully loaded phase.
- **Salt:** UCl₃–LiCl–KCl contracts 4A from 450 K, like LiCl–KCl [R19]. No numbers for other salts.
- **Cation form:** 5A (Na/Ca) retains LTA more readily than 4A under salt loading [R20]; no lattice data.

---

## 6. Mechanism notes

1. **Not uniform.** Contraction first proceeds within the original phase (450–625 K, solid-state/vapour occlusion), then a distinct small-a LTA phase nucleates and grows during molten occlusion (≥ 625–650 K) [R18].
2. **No amorphisation or framework decomposition** is seen up to 800 K in either APS study: no Si, SiO₂, Al₂O₃ or NaCl peaks [R19]; the new phase is still "zeolite" [R18].
3. **Sodalite.** Salt-loaded 4A tends to transform to sodalite; 5A resists [R20]. In the ANL waste-form route this happens during hot pressing with glass, and "the products depended on the initial salt loading" (the rate showed "large variations"); heated salt-occluded zeolite converted to "nepheline, salt, and sodalite", and glass–zeolite mixtures "were almost always converted to sodalite" [Q12] (grade A).
4. **Hydration is not the cause.** Dehydrated vs hydrated Na-A differ by only 0.12 % in a [Q1], and the APS samples were < 0.2 wt % H₂O [R18, R19].

---

## 7. Data file

`lta_salt_loading_lattice.csv` lists every lattice parameter above, one row each (row numbers match the tables here), plus the qualitative in-situ entries (rows 25–26), the bound (row 27) and the working value (row 28). Columns: `row, material, cation_form, salt, loading, T_K, a_A, a_esd_A, cell, method, grade, source_tag, note`. Blank T_K means "not stated".

---

## 8. References

Status records what was read this session. † = read by a research sub-agent; quote and URL in the session notes.

| Tag | Reference | Supports | Status |
|---|---|---|---|
| Q1 | S. Turner, J. R. Sieber, T. W. Vetter, R. Zeisler, A. F. Marlow, M. G. Moreno-Ramirez, M. E. Davis, G. J. Kennedy, W. G. Borghard, S. Yang, A. Navrotsky, B. H. Toby, J. F. Kelly, R. A. Fletcher, E. S. Windsor, J. R. Verkouteren, S. D. Leigh, "Characterization of chemical properties, unit cell parameters and particle size distribution of three zeolite reference materials: RM 8850 – zeolite Y, RM 8851 – zeolite A and RM 8852 – ammonium ZSM-5 zeolite," *Microporous Mesoporous Mater.* **107**, 252–267 (2008), doi:10.1016/j.micromeso.2007.03.019. https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=831471 | Hydrated and dehydrated Na-A cells; literature comparison | Full text read |
| Q2 | T. Carey, A. Corma, F. Rey, C. C. Tang, J. A. Hriljac, P. A. Anderson, "The effect of extra framework species on the intrinsic negative thermal expansion property of zeolites with the LTA topology," *Chem. Commun.* **48**(47), 5829–5831 (2012), doi:10.1039/C2CC30582H (pages from the RSC article page). ESI: https://www.rsc.org/suppdata/cc/c2/c2cc30582h/c2cc30582h.pdf | Ag-A (dry, hydrated) and ITQ-29 cells, 100–300 K | ESI read; Crossref rate-limited, bibliography from OpenAlex |
| Q3 | T. Carey, C. C. Tang, J. A. Hriljac, P. A. Anderson, "Chemical control of thermal expansion in cation-exchanged zeolite A," *Chem. Mater.* **26**(4), 1561–1566 (2014), doi:10.1021/cm403312q (confirmed by the reviewer via the Birmingham repository page). Accepted manuscript: https://research.birmingham.ac.uk/files/23935278/Alkali_metal_exchanged_Zeolite_A_paper_Final_Revised.pdf | α_V of dehydrated Na/K/Li/Rb/Cs-A | Accepted manuscript read† |
| Q4 | J. J. Pluth, J. V. Smith, *J. Am. Chem. Soc.* **102**, 4704–4708 (1980) (to confirm) | Dehydrated Na-A, 24.555(2) Å | Via [Q1] Table 17 and IZA-SC CIF snippet; not opened |
| Q5 | J. J. Pluth, J. V. Smith, "Crystal structure of dehydrated potassium-exchanged zeolite A. Absence of supposed zero-coordinated potassium. Refinement of Si,Al-ordered superstructure," *J. Phys. Chem.* **83**(6), 741–749 (1979), doi:10.1021/j100469a021 | Dehydrated K-A cells | COD records only†; DOI via Crossref† |
| Q6 | R. Y. Yanagida, A. A. Amaro, K. Seff (author order to confirm: the COD record lists Yanagida, Seff, Amaro), "Redetermination of the crystal structure of dehydrated zeolite 4A," *J. Phys. Chem.* **77**, 805–809 (1973) (DOI to confirm) | Dehydrated Na₁₂-A pseudo-cell 12.263 Å | COD record only† |
| Q7 | N. H. Heo, H. C. Choi, S. W. Jung, M. Park, K. Seff, "Complete redox exchange of indium for Tl⁺ in zeolite A. Crystal structures of anhydrous Tl₁₂-A and In₁₀-A·In. Indium appears as In²⁺, In⁺, and In⁰. The clusters (In₅)⁸⁺ and (In₃)²⁺ are proposed," *J. Phys. Chem. B* **101**(28), 5531–5539 (1997), doi:10.1021/jp970206v | Anhydrous Tl₁₂-A, 12.153(4) Å, 294 K | Repository abstract (POSTECH OASIS) |
| Q8 | H. S. Lim, J. Y. Kim, N. H. Heo, K. Seff, "Preparation, crystal structure, and luminescence properties of zeolite LTA containing extraframework tantalum(V), tantalum(II), thallium(I), and chloride," *J. Phys. Chem. C* **120**(22), 12139–12148 (2016), doi:10.1021/acs.jpcc.6b03943 | TaCl₆Tl₃-LTA, 12.099(1) Å, dealuminated framework | Abstract (x-mol) |
| Q9 | Y. Kim, M. S. Jeong, K. Seff, "Crystallographic studies of dehydrated Ag⁺ and K⁺ exchanged zeolite A reacted with alkali metal vapor," *Bull. Korean Chem. Soc.* **14**(5), 603–610 (1993), doi:10.5012/bkcs.1993.14.5.603 (DOI to confirm) | Dehydrated Ag₅.₆K₆.₄-A, 12.255(1) Å | Abstract |
| Q10 | C. Pereira, B. D. Babcock, "Fission product removal from molten salt using zeolite," ANL (1996). https://digital.library.unt.edu/ark:/67531/metadc682020/m2/1/high_res_d/379094.pdf (symposium paper, Scottsdale, 27–30 Oct 1996; no report number printed) | "approximately 11 to 13" Cl per unit cell, fully loaded | Full text read† |
| Q11 | D. Lexa, L. Leibowitz, A. J. Kropf, "On the reactive occlusion of the (uranium trichloride + lithium chloride + potassium chloride) eutectic salt in zeolite 4A," *J. Nucl. Mater.* **279**, 57–64 (2000) (as cited in [R18]; record to confirm) | May contain numeric 4A lattice parameters vs T | **Not read** (ScienceDirect) |
| Q12 | M. A. Lewis, M. C. Hash, C. Pereira, J. P. Ackerman, "Effects of heating on salt-occluded zeolite," ANL (1996) (report number to confirm), INIS abstract via https://worldwidescience.org/topicpages/s/salt+occluded+zeolite.html | 4A → nepheline/salt/sodalite on heating; products depend on initial salt loading | Abstract |
| — | Firor & Seff, *J. Am. Chem. Soc.* **100**, 3091 (1978); Pluth & Smith, *J. Am. Chem. Soc.* **104**, 6977 (1982); Kim, Subramanian, Firor, Seff, *ACS Symp. Ser.* **135**, ch. 7 (1980), doi:10.1021/bk-1980-0135.ch007; Kim & Seff, *Bull. Korean Chem. Soc.* **5**(3), 117 (1984) | Rows 14–19 (dry Ca, Sr, Ba, Cs/Ag forms; context only) | COD records / abstracts† |

Existing tags used: R18, R19, R20, R25, R27, R28, R30, R32 (see `lta_salt_occlusion_evidence.md` §8); Q26 is defined in `lta_coordination_numbers.md` §7. Bibliographies confirmed by review 3: R27 *J. Phys. Chem. C* **112**(30), 11181–11193; R30 *J. Phys. Chem. C* **122**(12), 6661–6668.

### Pages that could not be read

| Page | Reason |
|---|---|
| Richardson MRS 465 full text (OSTI purl 541788; Cambridge) | robots-blocked / paywalled |
| ANL "Relating structural parameters to leachability in a glass-bonded ceramic waste form" (INIS asyms-e3n93; OSTI DE00008917) — XRD lattice parameters of 4A/5A through waste-form production | robots-blocked |
| [Q11] J. Nucl. Mater. 279, 57 (2000) | ScienceDirect |
| Lexa, *Thermochim. Acta* (cited as "submitted" in [R18]) | not found in 4+ searches |
| [R32] Harward et al., *J. Nucl. Mater.* 2023 (super-occluded LiCl–KCl in 4A; likely XRD) | ScienceDirect |
| Crossref record for [Q2] | HTTP 429 (rate-limited) |
| COD CIFs (koreascience PDFs) | robots-blocked; record pages and abstracts used instead |
| Hand-off leads not reported in this revision: Lexa, *Rev. Sci. Instrum.* 70, 2242 (1999); Lexa, *J. Chem. Thermodyn.* (1999); R10 (Lewis, Fischer & Smith 1993); patents US 10,600,527 and US 5,613,240; Battisti et al. 2002; molten-nitrate comparison (Petranović; Liquornik & Marcus) | not reached or not searched |

**Full texts most worth obtaining:** [Q11], [R32], Richardson [R20], the ANL leachability report, and the APS raw data behind [R18] Fig. 3–4.

---

## 9. Change log

- **3.0 (2026-10-17):** new file. Independent-review changes are appended below.
- **3.0, review 1 fixes:**
  - §4: [R18]'s original phase is partly loaded above ~450 K, so a_new/a_orig at one T is only an upper bound on f_a, not a loaded/dry ratio.
  - §5a: the 17.3–17.8° window comes from axis labels only; the bound assumes all relevant peaks lie inside it and compares a_new(≤ 800 K) with a_orig(~450 K).
  - §5b: the working range now rests on the K/Br pairs 21/7 and 21/8; the rejected Tl pairs are context only.
  - §3: T of rows 7, 8, 23, 24 and hydration of row 21 flagged as not stated.
  - Row 6 space group "by size"; 0.14 %/0.12 % spread; sodalite impurity wording; Q3 title confirmed; Q6 author order flagged; row 24 grade "unverified".
- **3.0, review 2 (blocking-only):** no blocking errors in this file.
- **3.0, review 3 (re-check of sources review 1 could not reach, plus a consistency audit):** no lattice value contradicted.
  - Row 24 (R30) verified: 12.061(1) Å, grade A. Row 8 is K₁₁-A, Pm-3m (not "K₁₂", SG listed). Row 4 esd (3) and room temperature added.
  - §4 thermal term relabelled as an extrapolated size estimate (sign unsupported above 300 K), not a bound.
  - §6 Q12 quote corrected ("the products depended on the initial salt loading").
  - Bibliographies completed: Q2 pages 5829–5831, Q7, Q8, R27, R30; Q4, Q9 DOI, Q12 marked "to confirm"; Q26 and R32 tagged; hand-off leads not reported listed in §8.
- **3.0, final blocking-only check:** Q12 reference row wording aligned with §6 (products, not rate, depend on loading).
