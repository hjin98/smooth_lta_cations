# Coordination numbers of Li⁺, Na⁺, K⁺, Rb⁺, Cs⁺, Sr²⁺, Ba²⁺, Cl⁻ and O²⁻ in molten chlorides (and inside LTA)

**Revision 3.0, 2026-10-17 (new file).** Companion to `molten_salt_pair_distances.md` (P1–P26) and `lta_salt_occlusion_evidence.md` (R1–R34). New references are **Q-series**, shared with `lta_salt_loading_lattice_contraction.md` (Q1–Q12 are used there; this file uses Q13–Q28). §9 lists the independent-review changes.

Machine-readable table: `lta_coordination_numbers.csv` (one row per value; columns `centre, neighbour, environment, T_K, composition, CN, CN_uncertainty, cutoff_A, method, grade, source_tag, note`).

---

## 0. Conventions and how to use the numbers

**Definition.** N(A–B) is the mean number of B within the first coordination shell of an A centre. Centre is always written first.

**CNs depend strongly on the definition.** For the same neutron data on molten NaCl, Biggin & Enderby (via [P1] ESI) obtained **3.9–5.3** from g(r) and **4.0–4.5** from G(r), depending on which of four integration schemes was used (symmetric r·g(r), symmetric r²·g(r), integration of g(r) or r²·g(r) to the first minimum). Across all studies in [P1] Table S3, molten NaCl N(Na–Cl) spans **3.7–5.8**. **Never average across definitions.**

**Number conservation (used for all D rows).** In a melt with mole fractions x, N(A–B)·x_A = N(B–A)·x_B for the same cutoff. Hence N(Cl–M) = N(M–Cl) for MCl, and N(Cl–M) = N(M–Cl)/2 for MCl₂.

**For the structure generator — which values share the `r_min1` boundary.**

- **Only the §1a values** are defined by a first-shell boundary taken from the cation–anion g(r). For Li/Na/K, [P2] integrates to r†, "the boundary of the first solvation shell … defined from the cation-anion g(r)"; `r_min1` in the pair-distance file was taken from r†, so the two match by construction. For Rb/Cs, [P1] integrates WBK g(r) "to the first minimum", but the pair file's `r_min1` for Rb and Cs (4.74, 5.01 Å) is an E ratio estimate, not that WBK minimum.
- **§1b values** are neutron compilations whose integration rule is **not stated** in any source read.
- **§1c values** are EXAFS/XAFS first-shell fits (no cutoff: Rb, Sr) or AIMD integrations of the RDF to an R_min that [Q14] marks in its figures without giving a number or calling it the first minimum.
- Use §1a values as targets for the first-shell count inside `r_min1`. Use §1b and §1c values **with their spread**, not as hard targets inside `r_min1`. Do not compare any of them with diffraction CNs obtained by peak-symmetric integration.

| Grade | Meaning |
|---|---|
| **A** | Abstract, database record or search-result snippet |
| **B** | Full text (paper, ESI, author version, report) |
| **C** | Secondary: quoted by another source; original not opened |
| **D** | Own calculation (formula given) |
| **M / S / E** | Measured in the melt / simulation / estimate (as in the pair-distance file) |

"Read" means read in this session (2026-10-17). Values relayed by a research sub-agent are marked †; their quote and URL are in the session notes, and the independent review re-checked them (§9). Values marked **(rev 2.1)** were read from [P22] in an earlier session and are carried forward without being re-read.

---

## 1. Recommended values (selection rule stated per block)

### 1a. Pure alkali chlorides — cation-centred N(M–Cl) = anion-centred N(Cl–M)

**Selection rule:** polarizable-ion simulations integrated to the first-shell boundary of g₊₋(r) (a single, reproducible definition; see §0 for how it relates to `r_min1`). Use [P2] PIM at 1100 K for Li/Na/K (its r†, the first-solvation-shell boundary from g₊₋(r), is the pair file's `r_min1`); for Rb and Cs, where [P2] has no data, use the [P1] WBK polarizable model at the melting point. Uncertainty: ±0.5 (the [P1] stated simulation uncertainty), plus the definition spread in the last column.

| Salt | Recommended N(M–Cl) = N(Cl–M) | Basis | T (K) | Spread across all sources (definition-dependent) |
|---|---|---|---|---|
| LiCl | **4.3 ± 0.5** | [P2] PIM 4.34 (r† = 3.6 Å) | 1100 | 3.5–5.5 (exp.); 3.9–5.2 (sim.) |
| NaCl | **5.2 ± 0.5** | [P2] PIM 5.16 (r† = 4.1 Å); [P10] PIM 5.06 | 1100 | 3.7–5.8 (exp.); 4.76–5.3 (sim.) |
| KCl | **5.2 ± 0.5** | [P2] PIM 5.18 (r† = 4.5 Å) | 1100 | 3.5–6.1 (exp.); 4.89–5.8 (sim.; 4.7–5.8 with the two Takagi 1979 rows that [P1] S4 labels MD but its reference list calls experimental) |
| RbCl | **5.4 ± 0.5** | [P1] WBK 5.4 at T_m, 5.3 at T_m + 100 K; XAFS 5.3 [P3] | 991–1091 | 4.2–6.9 (exp.); 5.3–5.8 (sim.) |
| CsCl | **5.5 ± 0.5** | [P1] WBK 5.5 at T_m, 5.3 at T_m + 100 K | 918–1018 | 3–5.9 (exp.); 5.3–5.8 (sim.) |

Temperature dependence: **0 to −0.2 per 100 K** in [P1] WBK and JC (e.g. WBK LiCl 5.2 → 5.0, NaCl 5.3 → 5.1; JC LiCl 3.9 → 3.9).

The simulation spreads include every chloride row of [P1] Tables S1, S2 and S4, [P2] and [P10]. Some S4 "MC" values (Baranyai, Ruff & McGreevy 1986) may be fits constrained by diffraction data rather than free simulations; this was not checked.

Model dependence matters most for LiCl: the non-polarizable JC model gives 3.9, WBK 5.2, PIM 4.34, RIM 3.90 [P1, P2].

**Like-ion pairs (alkali chlorides):** within a first shell of like ions (to the first minimum of g_MM or g_ClCl), simulations give **N(M–M) ≈ N(Cl–Cl) ≈ 14–16** (WBK) or 11–17 (JC) [P1 S1/S2]; neutron data give 7.4–16.3 [P1 S3]. These are second-neighbour shells; use them only as loose checks.

### 1b. Pure alkaline-earth chlorides — N(M–Cl); N(Cl–M) = N(M–Cl)/2

**Selection rule:** neutron-diffraction values as compiled by Enderby (1987) and quoted in [P23] (read in full this session), cross-checked against [P22] (rev 2.1) and [Q13]. Where they conflict, the conflict is reported, not averaged.

| Salt | T (K) | Recommended N(M–Cl) | N(Cl–M) (D) | Sources | Status |
|---|---|---|---|---|---|
| MgCl₂ (calibration) | 998 | **4.3** | 2.15 | P23 (C), P22 (rev 2.1); AIMD 4.4 at 1200 K [P6, B]; PIM 5.21 → 4.72 with rising T [P10]; BELION MD 5.2 ± 0.9 [P23] | consistent |
| CaCl₂ (calibration) | 1093 | **5.3–5.4** | 2.65–2.7 | 5.3 [P23, C]; 5.4 [P22 rev 2.1; P8 quoting Biggin, citation defective]; MD 6.0 at 1200 K [P8] | consistent within 0.1 |
| SrCl₂ | 1198 | **6.9** | 3.45 | 6.9 [P23, C]; 6.9 "at 925 °C" quoted from McGreevy & Mitchell by [Q13] (C, second independent quotation); 6.9 in [P22] Table 2 but **6.6 in [P22] text** (rev 2.1) | two quotations of one measurement |
| BaCl₂ | 1298 | **unresolved: 6.4 or 7.7** | 3.20 or 3.85 | 6.4 [P23, C, via Enderby 1987]; 7.7 [P22 Table 2, rev 2.1]; "eightfold" [Q27, weak C] | **conflict open** |

**BaCl₂ conflict — what was established this session:**

1. Both numbers trace to the same neutron study, Edwards, Howe, Enderby & Page 1978 [Q20]. Its abstract (via OpenAlex) gives **no CN** and no temperature; the full text is closed. The correct DOI is **10.1088/0022-3719/11/6/012** (the hand-off's /014 resolves to an unrelated MgO paper)†.
2. [P23] states that the M–Cl CN "steadily increases from MgCl₂ to BaCl₂" and in the same sentence lists 4.3, 5.3, 6.9, **6.4**. The list contradicts its own "steadily increases" statement. The 7.7 of [P22] is consistent with the monotonic size trend; the 6.4 is not.
3. A citing snippet says Ba is "eightfold-coordinated to chlorine in the melt" [Q27]; its source is unknown.
4. A difference of this size is typical of different integration definitions (§0), but no source read here says which definition each value uses.

**Working recommendation for the generator (E):** N(Ba–Cl) = **7 ± 1** (covers 6.4–7.7 and the qualitative 8), with N(Cl–Ba) = 3.5 ± 0.5. Resolve with the Edwards PhD thesis (Leicester 1977, open on figshare [Q22]; the PDF could not be fetched here) or [P22] Table 2 footnotes.

**Like pairs in MCl₂:** no N(M–M) or N(Cl–Cl) was found in any readable source (§5).

### 1c. Mixed and dilute melts

| System | Pair | Recommended | T (K) | Basis | Grade |
|---|---|---|---|---|---|
| LiCl–KCl eutectic (nominal 59.2 mol % LiCl; [Q14] box 36 Li / 25 K / 64 Cl + 1 Ln³⁺ = 59.0 %) | Li–Cl | **4.2** (4.1–4.3) | 773 | AIMD; "integration of RDFs" to an R_min marked in figures (number not given; "first minimum" not stated); three runs (Sm, Ce, Eu) [Q14]. Source says "Li⁺ CN" without naming the partner; Cl assumed | B, S |
| same | K–Cl | **7.0** (6.8–7.2) | 773 | same; "K⁺ CN", partner not stated (presumably Cl). **High vs. pure KCl (5.2–5.4, at 1043–1143 K, i.e. ~300 K hotter): treat as definition-dependent; the temperature difference also contributes** | B, S |
| same | Cl–Li | 2.34 (box) / 2.47 (ideal melt) | 773 | D: 4.17 × 36/64; ideal Ln-free melt 4.17 × 0.592 | D |
| same | Cl–K | 2.75 (box) / 2.87 (ideal melt) | 773 | D: 7.03 × 25/64; ideal 7.03 × 0.408 | D |
| same | Cl–(Li+K) | 5.09 (box) / 5.34 (ideal melt) | 773 | D: sum; inherits the partner and R_min caveats above | D |
| same | Li–Cl (independent) | ≈ 4 | n/g | X-ray + MD, quoted by a citing paper (cited DOI 10.1023/a:1023310006090, to confirm)† | C |
| Rb⁺ dilute in LiCl–KCl (15 % RbCl) | Rb–Cl | **5.8** | 873 | XAFS; pure RbCl 5.3 at 1050 K by the same method [P3] | B, M |
| Sr²⁺ dilute in NaCl (5 mol % SrCl₂), melt | Sr–Cl | **5.8 ± 1.1** | 1273 | EXAFS; N = 5.8 ± 1.1 [P5 main text and SI Table S1]; R = 3.03 Å [SI Table S1 only]; 5.9 ± 0.89 at 1073 K [SI Table S1 only] (see note below) | B, M |
| Sr²⁺ (SrO dissolved) in LiCl–KCl | Sr–Cl, Sr–O | **6 Cl + 1 O** | n/g | EXAFS [Q15] | A, M |
| Sr²⁺ in SrCl₂–KCl–MgCl₂ | Sr–Cl | 6–7 (dominant fold) | 873–1173 | DP-MD [P19] | A, S |
| Cs⁺, Ba²⁺, Na⁺ dilute in LiCl–KCl or NaCl | M–Cl | **unavailable** | — | no numeric value in any readable source | — |

Observations from the data (not recommendations):

- **Dilution raises Rb–Cl**, from 5.3 (pure, 1050 K) to 5.8 (in LiCl–KCl, 873 K), with the same XAFS method [P3]. Part of the rise may be the 177 K lower temperature.
- **EXAFS Sr–Cl (5.8–5.9 in molten NaCl) is lower than the neutron value for pure SrCl₂ (6.9).** The methods differ (EXAFS first-shell fit with third cumulant vs. neutron partial g(r)), so this is not evidence of a real dilution effect.
- **For the generator, a dilute Sr²⁺ in LiCl–KCl at 773–823 K:** use **6 ± 1** Cl (E): it brackets [P5] (5.8–5.9, melt), [Q15] (6 Cl, LiCl–KCl) and [P19] (6–7). The [Q13] genetic-algorithm structures for the ≤ 10 % Sr subset of Sr–Na–Cl give 6.2 → 5.8 → 5.6 at 25 / 600 / 1000 °C (Table 3), so a slightly higher CN at lower T is plausible.

**Note on [P5] SI Table S1.** The SI was read in this session from the Clemson-hosted PDF of [P5] (https://clemson.edu/cecas/departments/mse/research/emrl/images/pdfs/jphyschem_in-situ-determination-of-speciation-and-local-structure.pdf). In that fetch the SI followed the main text as pages S1–S3, and Table S1 is on p. S3. A later fetch of the same URL by the independent reviewer returned only the main text (ending at the references), and the ACS SI link failed, so the reviewer **could not re-verify** the SI values. Main-text values (7.5 ± 1.3 at 25 °C, 5.8 ± 1.1 at 1000 °C) were verified by both; 6.5 ± 1.9 at 600 °C is also quoted in [Q13] Table 3 (C). The 400 °C and 800 °C N values and all R values rest on the single in-session reading. Melting point of the 5 mol % mixture: 782 ± 2 °C (= 1055 K; [P5] Table 1, origin not stated), so only the 800 and 1000 °C points are melt.
- **Trivalent reference** (not a target species): Y³⁺ (15 mol % YCl₃) in LiCl–KCl at 973 K, 6.0 ± 0.3 (EXAFS; 6.1 and 5.9 at 30 and 50 mol %) [Q19]; dilute M³⁺ CN "overwhelmingly likely to be 6 in chloride melts" [Q28].

### 1d. Oxide O²⁻ in chloride melts

| Pair | Value | Basis | Grade |
|---|---|---|---|
| Sr–O (SrO in LiCl–KCl and LiCl–LiF) | **1** O per Sr, plus 6 Cl | EXAFS [Q15]. **SrOCl** observed (298–923 K); its greater stability than SrO and SrCl comes from the DFT/thermodynamic estimates | A |
| Li–O, O–Li (Li₂O in LiCl) | **unavailable** | [P13] and [P15] are the direct sources; abstracts and full texts unreadable | — |
| Ca–O, O–Ca (CaO in CaCl₂) | **unavailable** (number) | speciation only: below 1104 K, CaO + CaCl₂ → **Ca₄OCl₆**, which then dissociates; above 1104 K "CaO directly dissociates to produce O²⁻" [Q16]; Ca²⁺ forms "a positive coordination layer" around O²⁻ (MD, 1073–1473 K) [Q17] | A |
| Ba–O, Na–O, K–O, Mg–O in chloride melts | **unavailable** | none found | — |
| Li₂O in LiCl–Li | — | Li and Li₂O "behave independently over a broader range of melt compositions" [Q18] | A |

**Speciation guidance (E):** in the data read, O²⁻ in chloride melts appears as a cation-bound oxychloride unit (SrOCl in alkali chlorides [Q15]; Ca₄OCl₆ below 1104 K in CaCl₂ [Q16]) rather than as a free ion surrounded by Cl⁻. For the generator, O²⁻ should be placed with a first shell of cations only (use the O-bridged rows of the pair-distance file). No source gives N(O–M), so the size of that shell is **unavailable**.

---

## 2. First-hand data: pure alkali chlorides

### 2.1 Experiment ([P1] ESI Table S3; originals as listed there, not opened) — grade C, M

| Salt | T (K) | N(M–Cl) | Like-pair N | Method | Original |
|---|---|---|---|---|---|
| LiCl | 923 | 4 | 12 (mainly Cl–Cl) | X | Ohno et al. 1978 |
| LiCl | n/g | 4 | 12 | X | Levy et al. 1960 |
| LiCl | 893 | 4.1 | 9 | X | Zarzycki 1958 |
| LiCl | 958 | 5.5 | 14 (Li–Li), 14 (Cl–Cl) | N | Howe & McGreevy 1988 |
| LiCl | n/g | 3.5 | 8.3, 8.3 | N | Levy et al. 1960 |
| NaCl | 1083 | 3.7 | — | X | Ohno & Furukawa 1981 |
| NaCl | 1093 | 4.7 | 9 | X | Zarzycki 1958 |
| NaCl | 1148 | 5.8 | 13, 13.3 | N | Edwards et al. 1975 |
| NaCl | 1148 | 4.0 | — | N | same data, reinterpreted by Biggin & Enderby 1982 |
| KCl | n/g | 4.1 / 3.7 / 5.8 | — | X | Takagi 1979 / Levy 1960 / Lark & Miller 1937 |
| KCl | 1083 | 5.2 | 10 | X | Zarzycki 1958 |
| KCl | 1073 | 6.1 | 12/15.9, 12.3/16.2 | N | Derrien & Dupuy 1975 |
| KCl | n/g | 3.5 | 7.4, 7.4 | N | Levy et al. 1960 |
| RbCl | n/g | 4.2 | — | X | Levy & Danford 1964 |
| RbCl | 1023 | 6.9 | 13, 14 | N | Mitchell, Poncet & Stewart 1976 |
| CsCl | n/g | 4.6 | 7.1 | X | Levy et al. 1960 |
| CsCl | 968 | 5.8 | 15.4, 16.3 | N | Locke et al. 1985 |
| CsCl | 973 | 3 / 5.9 | — | N | Derrien & Dupuy 1975 |

Also: NaCl 4.83 (experiment, quoted from McGreevy 1987 by [P10]; C).

### 2.2 Simulation — grade B, S ([P1] S1/S2, [P2], [P10]); [P1] S4 entries grade C, S

| Salt | [P1] WBK, T_m / T_m+100 K | [P1] JC, T_m / T_m+100 K | [P2] PIM / RIM, 1100 K (r†, Å) | Other |
|---|---|---|---|---|
| LiCl | 5.2 / 5.0 | 3.9 / 3.9 | 4.34 / 3.90 (3.6) | MD/MC 4.0–4.9 (7 rows) [P1 S4] |
| NaCl | 5.3 / 5.1 | 4.9 / 4.8 | 5.16 / 4.82 (4.1) | PIM 5.06 at 1100 K [P10]; MD/MC 4.76–5.2 [P1 S4] |
| KCl | 5.4 / 5.2 | 5.4 / 5.2 | 5.18 / 4.89 (4.5) | MD/MC 5.13–5.8; Takagi 4.7, 5.1 [P1 S4] |
| RbCl | 5.4 / 5.3 | 5.5 / 5.4 | — | MD 5.41; MC 5.7, 5.8 [P1 S4] |
| CsCl | 5.5 / 5.3 | 5.8 / 5.6 | — | MD 5.7; MC 5.8 [P1 S4] |

T_m (K): LiCl 878.15, NaCl 1074.15, KCl 1043.15, RbCl 991.15, CsCl 918.15 [P1].

- [P2]'s 4.34 / 5.16 / 5.18 are the **PIM** set only; RIM gives 3.90 / 4.82 / 4.89. [P2] reports no AIMD CN and no anion-centred CN (correction to the hand-off).
- [P2] ESI also lists free-energy-well CNs (n_R/n_P) for its rate theory. These are not averages and must not be used as CNs.

---

## 3. First-hand data: alkaline-earth chlorides

| Salt | T (K) | Value | Method | Source | Grade |
|---|---|---|---|---|---|
| MgCl₂ | 998 | 4.3 | neutron (Biggin, Gay & Enderby 1984) | P23 quoting Enderby 1987 | C |
| MgCl₂ | 998 | 5.2 ± 0.9 (SD of distribution); 19.5 % of Mg four-coordinate | BELION MD fitted to the diffraction data | P23 | B, S |
| MgCl₂ | 1200 | 4.4 (cutoff not stated for the pure salt) | AIMD | P6 | B, S |
| MgCl₂ | not given | 5.21 → 4.72 with rising T | PIM, first minimum | P10 | B, S |
| CaCl₂ | 1093 | 5.3 | neutron (Biggin & Enderby 1981) | P23 | C |
| CaCl₂ | 1093 | 5.4 | as above | P22 (rev 2.1); P8 (citation defective) | C |
| CaCl₂ | 1200 | 6.0 (integration limit not stated) | classical BMH MD | P8 | B, S |
| SrCl₂ | 1198 | 6.9 | neutron (McGreevy & Mitchell 1982) | P23; Q13 | C |
| SrCl₂ | 298 | 6.61 | USPEX/DFT crystal ensemble, Boltzmann-weighted (room T only; not a liquid simulation). Q13 also prints 5.5 (500 °C) and 4.88 (1000 °C), but in a sentence that follows its Zr–F discussion, so the attribution to Sr–Cl is not secure; omitted | Q13 | B, S |
| BaCl₂ | 1298 | 6.4 | neutron (Edwards et al. 1978) | P23 | C |
| BaCl₂ | 1298 | 7.7 | as above | P22 (rev 2.1) | C |

Experimental uncertainty for diffraction CNs from partial g(r): about ±0.2 (Waseda & Toguri, quoted in [P23]).

---

## 4. Inside LTA (optional block): single-crystal coordination

These come from dehydrated crystals or from solution/vapour loading, not from melt equilibrium.

| Centre | Neighbours | CN | Distances | System | Grade | Ref |
|---|---|---|---|---|---|---|
| Sr²⁺ (6-ring) | framework O3 | **3** | Sr–O3 2.21(2), 2.31(1) Å; O3–Sr–O3 118.7°, 120.0° | dehydrated Sr₁.₆Tl₈.₈-A, Sr₅.₄₅Tl₁.₁-A, 294 K | A† | Q25 |
| Ba²⁺ | framework O | **3** (forms on dehydration; may disrupt the framework) | — | dehydrated Ba-A series | A† | Q26 |
| Ag⁺ (6-ring) | framework O | 3 | 2.52(3) Å | K₉(K₄Br)·0.75Ag₄Br₄-A | A | R28 |
| Cs⁺ (8-ring centre, D₄h) | framework O | not stated | about 0.3 Å longer than the ionic-radius sum | Cs₇.₃Ag₄.₇-A | A† | Kim & Seff 1984 (full reference in the lattice file §8) |
| Br⁻ (α-cage, opposite a 4-ring) | K⁺ | **4** (K₄Br³⁺: 2 eight-ring + 2 six-ring K⁺) | — | K₉(K₄Br)·0.75Ag₄Br₄-A | A | R28 |
| Br⁻ (β-cage Ag₄Br₄) | K⁺ through 6-ring | 1 | 3.00(4) Å | same | A | R28 |
| I⁻ (α-cage) | K⁺ | 4 (K₄I³⁺) | — | K₉(K₄I)(Ag₄I₄)₀.₅-A | A | R34 |
| Cl⁻ (α-cage) | K⁺ | **3** (K₃Cl²⁺, by formula) | — | R27 composition | A | R27 |
| Cl⁻ (β-cage Ag₄Cl₄) | K⁺ through 6-ring | 1 | 2.70(3) Å | same | A | R27 |
| Tl⁺ (α-cage centre) | Cl⁻ | 6 (octahedral) | — | (WCl)₀.₆(Tl₇Cl₆)₀.₅Tl₈.₅-A | A | R30 |
| Cl⁻ in Tl₇Cl₆⁺ | Tl⁺ | 2 (central Tl + one 8-ring Tl) | — | same | A | R30 |

**Gaps:** Na⁺–O CN at the 6-, 8- and 4-ring sites of dehydrated Na₁₂-A, and K⁺/Cs⁺–O CNs in fully exchanged dehydrated K-A and Cs-A, were not reached. Titles record a dispute over "near-zero-coordinate" Sr²⁺ (Firor & Seff 1978) versus its absence (Pluth & Smith 1982)† (references in the lattice file §8).

**Comparison with the melt:** confined α-cage halides in these crystals have **3–4** cation neighbours (K₃Cl²⁺, K₄X³⁺), fewer than the eutectic-melt value N(Cl–(Li+K)) = 5.1–5.3 (D, §1c; inherits the [Q14] definition caveats). These crystals were not melt-equilibrated, so this is an observation, not a calibration target.

---

## 5. Data availability

| Status | Pairs |
|---|---|
| Measured or simulated, recommended | M–Cl for Li, Na, K, Rb, Cs (pure); Mg, Ca, Sr (pure); Li–Cl, K–Cl in LiCl–KCl; Rb–Cl in LiCl–KCl; Sr–Cl dilute in NaCl and (with oxide) in LiCl–KCl |
| Derived (D) | every N(Cl–M) in §1 |
| **Conflict open** | Ba–Cl in pure BaCl₂ (6.4 vs 7.7) |
| **Unavailable** | Cs–Cl, Ba–Cl, Na–Cl dilute in LiCl–KCl or NaCl; Sr–Cl dilute in LiCl–KCl without oxide; all N(M–M) and N(Cl–Cl) in MCl₂ and in LiCl–KCl; every M–O and O–M except Sr–O = 1 |

---

## 6. Recommended follow-ups

1. **Li–O in LiCl:** obtain [P13] and [P15] (ScienceDirect) through the library.
2. **BaCl₂:** Edwards thesis [Q22] (figshare, 37.9 MB) or [P22] Table 2 footnotes.
3. **Eutectic partials:** Guo et al. [Q23] (machine-learned LiCl–KCl potential validated by HEXRD) should give all partial CNs; or compute them with the user's own MACE model on the 59.2 : 40.8 eutectic at 773–823 K, integrating to `r_min1`.
4. **Dilute Cs⁺, Ba²⁺ in LiCl–KCl:** [P17] (Sr²⁺/Cs⁺ DP-MD) full text; otherwise compute.
5. **Consistency:** when the generator counts first-shell neighbours, use `r_min1` from the pair-distance file. Only the §1a Li/Na/K values are defined with that boundary; for Rb/Cs, `r_min1` is an estimate (§0).

---

## 7. References (new Q-series; existing P/R tags as in the companion files)

Status records what was read this session; † = read by a research sub-agent.

| Tag | Reference | Used for | Status |
|---|---|---|---|
| Q13 | M. S. Christian, T. J. Lynch, J. Schorne-Pinto, A. M. Mofrad, N. R. Birkner, K. S. Brinkman, W. K. S. Chiu, T. M. Besmann, "Modeling metallic halide local structures in salt melts using a genetic algorithm," *J. Phys. Chem. C* **126**(22), 9239–9247 (2022), doi:10.1021/acs.jpcc.2c00747 (Crossref). Author PDF: https://www.clemson.edu/cecas/departments/mse/research/emrl/images/pdfs/acs.jpcc_modeling-metallic-halide-local-structures.pdf | SrCl₂ 6.9 quoted from McGreevy & Mitchell; USPEX Sr–Cl trends | Full text read |
| Q14 | J. Fuller, J. Moon, J. Zhang, D. Chidambaram, Q. An, "Coordination and thermophysical properties of select trivalent lanthanides in LiCl-KCl," *Phys. Chem. Chem. Phys.* **24**(21), 13102–13109 (2022), doi:10.1039/D1CP03056F. https://pubs.rsc.org/en/content/getauthorversionpdf/d1cp03056f | Li–Cl, K–Cl in the eutectic, 773 K | Author version read† |
| Q15 | D. Kang, C. Kwon, W. Yang, S. Yoon, Y. Lee, J. T. M. Amphlett, S.-E. Bae, S. Kim, S. Choi, "Spectroscopic and theoretical analyses of the reaction of SrO in molten chloride and fluoride salts," *J. Nucl. Mater.* **592**, 154962 (2024), doi:10.1016/j.jnucmat.2024.154962 | Sr: 1 O + 6 Cl; SrOCl | Abstract (Scopus record)† |
| Q16 | C. Lin, X. Yang, "Insights into the dissolution behavior of CaO and mass transfer of oxide ions in molten CaCl₂," *Metall. Mater. Trans. B* (2024), doi:10.1007/s11663-024-03123-0 (volume/pages to confirm) | CaO speciation; Ca₄OCl₆ | Abstract† |
| Q17 | Dan L., Liu N., Li Z., Li X., "Study on the diffusion of oxygen ions in calcium chloride molten salt," *Rare Met. Mater. Eng.* **50**(7), 2409–2414 (2021), doi:10.12442/j.issn.1002-185X.20200631 | Ca²⁺ shell around O²⁻ (MD) | Abstract† |
| Q18 | A. Merwin, "Material interactions with molten LiCl-Li₂O-Li," PhD dissertation, Univ. Nevada Reno (2016). https://scholarwolf.unr.edu/items/0bd22049-90ce-43d0-adad-0dba73bff2a6 | Li and Li₂O independent | Abstract† |
| Q19 | Y. Okamoto, M. Akabori, H. Motohashi, H. Shiwaku, T. Ogawa, "X-ray absorption study of molten yttrium trihalides," *J. Synchrotron Rad.* **8**, 1191–1199 (2001) (DOI to confirm). https://journals.iucr.org/s/issues/2001/06/00/fh2004 | Y–Cl in LiCl–KCl (trivalent reference) | Full text† |
| Q20 | F. G. Edwards, R. A. Howe, J. E. Enderby, D. I. Page, "The structure of molten barium chloride," *J. Phys. C* **11**(6), 1053–1057 (1978), **doi:10.1088/0022-3719/11/6/012** (corrected) | Original BaCl₂ partials | Abstract only (no CN)† |
| Q21 | J. E. Enderby, "The structure of molten salts," in *Molten Salt Chemistry* (Reidel/Springer, 1987), pp. 1–15, doi:10.1007/978-94-009-3863-2_1 | Compilation behind P23's CNs | Not opened; bibliography confirmed† |
| Q22 | F. G. Edwards, "The structure of molten salts," PhD thesis, Univ. Leicester (1977). https://figshare.le.ac.uk/articles/thesis/The_structure_of_molten_salts_/10185695 | BaCl₂ at 1025 °C (1298 K); likely the CN definition | Record read; **PDF not readable** (fetch error) |
| Q23 | J. Guo et al., "Composition-transferable machine learning potential for LiCl-KCl molten salts validated by high-energy x-ray diffraction," *Phys. Rev. B* **106**, 014209 (2022), doi:10.1103/PhysRevB.106.014209; ChemRxiv v1 (4 Feb 2022) | Eutectic partials (future) | ChemRxiv page read; no CNs on it |
| Q24 | T. Ohkubo, E. Tsuchida, M. Gobet, V. Sarou-Kanian, C. Bessada, Y. Iwadate, *J. Phys. Chem. B* **117**(18), 5668–5674 (2013), doi:10.1021/jp312486m | Shows that P14 is mis-cited: this is xLi₂O–(1−x)**B₂O₃**, not LiCl | Crossref record† |
| Q25 | J. W. Yang, J. Y. Park, U. S. Kim, Y. Kim, "Two crystal structures of dehydrated Sr²⁺ and Tl⁺ exchanged zeolite A, SrₓTl₁₂₋₂ₓ-A (x = 1.6 and 5.45)," *Bull. Korean Chem. Soc.* **10**(6), 582–585 (1989), doi:10.5012/bkcs.1989.10.6.582; data at 21(1) °C | Sr²⁺ three-coordinate at 6-rings | Abstract† |
| Q26 | Y. Kim, V. Subramanian, R. L. Firor, K. Seff, *ACS Symp. Ser.* **135**, ch. 7 (1980), doi:10.1021/bk-1980-0135.ch007 | Three-coordinate Ba²⁺ in dehydrated Ba-A | Abstract† |
| Q27 | Edgar et al., *IEEE Trans. Nucl. Sci.* (2010) (details to confirm), as shown in a citing snippet on the ResearchGate page of Q20 | "eightfold-coordinated" Ba in molten BaCl₂ | Snippet only† |
| Q28 | R. Corradini, P. A. Madden, M. Salanne, "Coordination numbers and physical properties in molten salts and their mixtures," *Faraday Discuss.* (2016), doi:10.1039/C5FD00223K (volume/pages to confirm). Accepted manuscript: https://pubs.rsc.org/en/content/getauthorversionpdf/c5fd00223k | Dilute M³⁺ CN = 6 (trivalent reference) | Accepted manuscript; quote verified by review 3 |

**Existing references, with status updated this session:**

| Tag | Update |
|---|---|
| P1 | ESI Tables S1–S4 re-read. **Correction:** the tables are for pure alkali halides, not LiCl–CsCl mixtures as the hand-off said |
| P2 | Author version re-read†; values are PIM only; RIM values added |
| P3 | Re-read† |
| P5 | **SI Table S1 read once** in this session (distances, σ², N at 25–1000 °C), from the Clemson-hosted PDF; the reviewer's re-fetch did not return the SI (see note in §1c) |
| P6, P8, P10 | Re-read (P8 citation of Biggin is defective) |
| P11 | *J. Mol. Liq.* **238**, 236–247 (2017) (pages now confirmed†); closed; no CNs read |
| P13 | Authors confirmed (Crossref)†: L. Lin, J. Song, C. Zhu, S. Guo, W. Zhou, Y. Wang, "Li₂O concentration influenced local structure and properties of molten LiCl salt by machine learning driven molecular dynamics simulation," *J. Mol. Liq.* **413**, 125950 (2024). **Not read** |
| P14 | **Mis-cited** — see Q24. Drop |
| P15 | Confirmed†: J. Xu, Y. Wang, B. Yao, Y. Jia, Y. Xiao, L. Zhang, B. Li, H. He, B. Yue, L. Yan, *Nucl. Eng. Des.* **438**, 114052 (2025), doi:10.1016/j.nucengdes.2025.114052. **Not read** |
| P17 | X. Zhang, J. Song, L. Lin, Y. Ding, X. Xu, W. Zhou, S. Guo, D. Wang, *J. Radioanal. Nucl. Chem.* (online 27 Sep 2024)†; volume/pages 334(1), 329–341 to confirm; abstract has no CN |
| P19 | J. Zhao, T. Feng, G. Lu, *J. Chem. Theory Comput.* **20**(17), 7611–7623 (2024), **doi:10.1021/acs.jctc.4c00824**†; abstract read |
| P21 | Data (Zenodo) not reachable from this sandbox; no CN computed |
| P22 | **Not re-read** (paywalled; HAL blocked). Values marked (rev 2.1) |
| P23 | Full text re-read this session (pyrometallurgy.co.za PDF) |

### Pages that could not be read

| Page | Reason |
|---|---|
| P13, P15 (ScienceDirect); P11 | robots-blocked / closed access |
| P22 (RSC; HAL copy) | paywall; HAL bot wall |
| Q20 full text (IOP); NASA ADS record | IOP blocks automated fetch; ADS robots-blocked |
| Q22 thesis PDF (figshare ndownloader) | fetch error |
| P17 full text (Springer) | paywall; OpenAlex rate-limited |
| P21 Zenodo data; Crossref/Zenodo from the shell | sandbox egress blocked |
| Crossref, OpenAlex, Europe PMC | intermittent HTTP 429 (rate limits) |
| Guo 2022 full text, Okamoto 2004 (molten SrCl₂ XAFS), Zhou 2024 (Cs in LiCl–KCl) main text, Zhou 2025 *Chin. J. Chem. Eng.* (Li₂O in LiCl–UCl₃) | not reached / not accessible |
| Hand-off leads not reported in this revision: Marcus 2016 review (used only as cited in [P1]); Ohno et al. 1994 book; [P22] Table 1 alkali CNs (P22 not re-read) | not reached or not searched |
| [P15] and [P17] bibliographic re-check; [Q27] original source | rate limits; two searches found no source for the "eightfold" claim |

---

## 8. Data file columns

`centre, neighbour, environment, T_K, composition, CN, CN_uncertainty, cutoff_A, method, grade, source_tag, note`. `CN` is text where the source gives a range or a qualitative value (e.g. "6-7 (dominant)", "unavailable"). Blank `T_K` means not stated. `grade` holds the source grade, then the kind of value after a comma (e.g. `B,S`, `C,M`; M = measured, S = simulated, E = estimate rows carry grade `E`; D rows have no kind; `-` = no value; a bare `C` = method unknown, used only for the [Q27] "eightfold" row). `cutoff_A` holds only an integration cutoff or boundary; for EXAFS/XAFS fits it says "none (first-shell fit)" and for crystal structures "n/a", with the fitted distances in `note`. The CSV holds every value in §§1–4 (215 rows): all chloride rows of [P1] Tables S1–S4, all D rows, the three E working values (Ba–Cl 7 ± 1, Cl–Ba 3.5 ± 0.5, dilute Sr–Cl 6 ± 1), and an "unavailable" row for every pair §5 lists as unavailable.

---

## 9. Change log

- **3.0 (2026-10-17):** new file. Independent-review changes are appended below.
- **3.0, review 1 fixes:**
  - Removed the [Q13] "binary Sr–Cl" 5.5 (773 K) and 4.88 (1273 K) values: in the source they follow the Zr–F sentence. Only 6.61 (298 K) is kept. The dilute-Sr reasoning now cites the [Q13] Table 3 subset trend 6.2 → 5.8 → 5.6.
  - §0 rewritten: only §1a shares the first-shell boundary; Rb/Cs `r_min1` is an estimate; §1b definitions unknown; §1c are EXAFS fits or AIMD to an unnumbered R_min. [Q14] cutoff wording corrected; partner of the "K⁺ CN" flagged.
  - Eutectic D rows now give both the [Q14] box stoichiometry (36/25/64: 2.34, 2.75, 5.09) and the ideal melt (2.47, 2.87, 5.34).
  - [P5] SI provenance stated (single in-session reading; not re-verifiable by the reviewer); 782 °C origin "not stated".
  - RbCl simulation spread 5.3–5.5; temperature dependence 0 to −0.2 per 100 K; KCl MC 5.8 noted; "full [P1] set" → selected S4 rows; "supported twice" → "two quotations of one measurement"; Q25 bibliography completed. CSV now 181 rows.
- **3.0, review 2 (blocking-only):** no numerical errors; R = 3.03 Å re-attributed to SI Table S1 only (Table S1 page S3 confirmed from the session's fetch record).
- **3.0, review 3 (full re-check of sources review 1 could not reach, plus a consistency audit):**
  - [P1] Table S4 read in full: RbCl simulation spread restored to **5.3–5.8** (review 1's narrowing to 5.3–5.5 rested on a partial reading); KCl 4.89–5.8 (Baranyai MC 5.8 verified); all 23 chloride S4 rows now in the CSV (LiCl 4.0–4.9).
  - Neutron like-pair range 7.4–16.3 (was 8.3–16.3). KCl eutectic comparison now notes the ~300 K temperature difference.
  - [Q15]: SrOCl stability attributed to DFT/thermodynamics, not Raman; study range 298–923 K. [Q18]: "(Raman)" removed. [Q19]: 15 mol % YCl₃ stated.
  - CSV: E working values and missing quoted values added (Li–Cl ≈ 4, Y–Cl, dilute M³⁺, Cl–Ca 2.70); "unavailable" rows for every unavailable pair; oxide rows re-tagged (Ca → Q16/Q17; Ba/Na/K/Mg → none found); M/S/E kind added to `grade`; fitted distances moved from `cutoff_A` to `note`. 215 rows.
  - New tag Q28 (Corradini, Madden & Salanne). Bibliographies completed for Q14, Q15, Q18, Q19, P17, P19; hand-off leads not reported listed in §7.
- **3.0, final blocking-only check:** Q15 rows re-graded `A,M` (EXAFS; were mislabelled S); bare-`C` grade for Q27 documented.
