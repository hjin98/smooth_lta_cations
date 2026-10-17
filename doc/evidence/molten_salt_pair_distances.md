# Pair distances in molten salts: reference table for the random molten-salt structure generator

**Revision 2.1, 2026-10-17.** This revision follows two independent reviews. §8 lists what changed and why.

Species covered:

- **Cations:** Li⁺, Na⁺, K⁺, Rb⁺, Cs⁺, Sr²⁺, Ba²⁺
- **Anions:** Cl⁻, O²⁻ (free oxide ion)

All 45 distinct pairs are covered. Like-sign pairs depend on the bridging counter-ion, so they are given per environment: cation pairs Cl-bridged and O-bridged, anion pairs around each cation. That gives 91 rows in total. Running the script overwrites the CSV.

| File | Contents |
|---|---|
| `molten_salt_pair_distances.csv` | Machine-readable table |
| `molten_salt_pair_table.py` | Generates the CSV; holds every input constant and formula, and prints the fitted constants when run |

## 0. Definitions and how the generator should use them

| Quantity | Definition | Generator use |
|---|---|---|
| **r_cut** | Closest approach: smallest r at which g(r) becomes non-zero (onset, g ≈ 0.01–0.1) | Hard-core rejection distance for random insertion. **Provisional estimate**, chosen so random sequential insertion reaches experimental density (§4.3) |
| **r_peak** | Position of the first maximum of g(r) | Target nearest-neighbour distance |
| **r_min1** | Position of the first minimum after the peak | First-shell boundary for coordination checks |

| Grade | Meaning |
|---|---|
| **M** | Measured in the molten salt (neutron/X-ray diffraction or XAFS), via the compilation or abstract cited |
| **S** | Molten-salt simulation (MD, polarizable-ion model (PIM), AIMD) read from the source |
| **E** | Estimate from the model in §4, not first-hand; each CSV row states its formula |

Grades refer to r_peak. Every r_cut is E. Every r_min1 is E except Li/Na/K–Cl, which are S.

**Reference temperature.** Values apply near the melting point, about 950–1300 K. In the polarizable simulations of [P1], temperature moves the peaks by:

| Pair type | Shift per 100 K |
|---|---|
| Cation–Cl⁻ | ≤ 0.02 Å |
| Like-ion pairs in the chlorides | up to +0.06 Å |

**No first-hand onset (r_cut) or measured first minimum could be read from any accessible source.** The raw partial g(r) of molten NaCl (neutron diffraction plus PIM) is openly published [P4]. Use it to calibrate the onset ratios (§6).

---

## 1. Recommended values (all pairs)

### 1a. Cation–anion

| Pair | Environment | r_cut (Å) | r_peak (Å) | r_min1 (Å) | Grade (peak) |
|---|---|---|---|---|---|
| Li–Cl | pure/any | 1.95 | 2.38 | 3.60 | M |
| Na–Cl | pure/any | 2.25 | 2.75 | 4.10 | M |
| K–Cl | pure/any | 2.53 | 3.08 | 4.50 | M |
| Rb–Cl | pure/any | 2.61 | 3.18 | 4.74 | M |
| Cs–Cl | pure/any | 2.78 | 3.39 | 5.01 | M |
| Sr–Cl | pure/any | 2.38 | 2.90 | 4.17 | M |
| Ba–Cl | pure/any | 2.54 | 3.10 | 4.46 | M |
| Li–O | pure/any | 1.62 | 1.98 | 2.85 | E |
| Na–O | pure/any | 1.95 | 2.38 | 3.42 | E |
| K–O | pure/any | 2.27 | 2.77 | 3.98 | E |
| Rb–O | pure/any | 2.38 | 2.90 | 4.17 | E |
| Cs–O | pure/any | 2.43 | 2.96 | 4.26 | E |
| Sr–O | pure/any | 2.09 | 2.55 | 3.67 | E |
| Ba–O | pure/any | 2.25 | 2.74 | 3.94 | E |

### 1b. Cation–cation sharing a Cl⁻ (chloride environment)

| Pair | Environment | r_cut (Å) | r_peak (Å) | r_min1 (Å) | Grade (peak) |
|---|---|---|---|---|---|
| Li–Li | Cl-bridged | 2.70 | 3.75 | 5.44 | M |
| Li–Na | Cl-bridged | 2.75 | 3.82 | 5.53 | E |
| Li–K | Cl-bridged | 2.94 | 4.09 | 5.92 | E |
| Li–Rb | Cl-bridged | 3.00 | 4.17 | 6.04 | E |
| Li–Cs | Cl-bridged | 3.13 | 4.35 | 6.30 | E |
| Li–Sr | Cl-bridged | 2.65 | 3.69 | 5.35 | E |
| Li–Ba | Cl-bridged | 2.76 | 3.84 | 5.57 | E |
| Na–Na | Cl-bridged | 2.83 | 3.93 | 5.70 | M |
| Na–K | Cl-bridged | 3.12 | 4.33 | 6.28 | E |
| Na–Rb | Cl-bridged | 3.18 | 4.41 | 6.40 | E |
| Na–Cs | Cl-bridged | 3.30 | 4.58 | 6.64 | E |
| Na–Sr | Cl-bridged | 2.83 | 3.93 | 5.69 | E |
| Na–Ba | Cl-bridged | 2.93 | 4.07 | 5.90 | E |
| K–K | Cl-bridged | 3.43 | 4.77 | 6.92 | M |
| K–Rb | Cl-bridged | 3.35 | 4.65 | 6.74 | E |
| K–Cs | Cl-bridged | 3.46 | 4.81 | 6.97 | E |
| K–Sr | Cl-bridged | 2.99 | 4.16 | 6.03 | E |
| K–Ba | Cl-bridged | 3.09 | 4.29 | 6.23 | E |
| Rb–Rb | Cl-bridged | 3.50 | 4.86 | 7.05 | M |
| Rb–Cs | Cl-bridged | 3.51 | 4.88 | 7.07 | E |
| Rb–Sr | Cl-bridged | 3.04 | 4.23 | 6.13 | E |
| Rb–Ba | Cl-bridged | 3.14 | 4.36 | 6.33 | E |
| Cs–Cs | Cl-bridged | 3.56 | 4.95 | 7.18 | M |
| Cs–Sr | Cl-bridged | 3.16 | 4.38 | 6.36 | E |
| Cs–Ba | Cl-bridged | 3.25 | 4.51 | 6.54 | E |
| Sr–Sr | Cl-bridged | 2.70 | 3.76 | 5.45 | E |
| Sr–Ba | Cl-bridged | 2.80 | 3.89 | 5.64 | E |
| Ba–Ba | Cl-bridged | 2.89 | 4.01 | 5.82 | E |

**How reliable the estimated rows are:**

- **Alkali–alkali pairs:** about −0.2 to +0.2 Å. This is the spread of the alkali bridge factor F (1.005–1.095, §4.4).
- **Alkali–Sr/Ba pairs:** about ±0.3 Å.
- **Sr–Sr, Ba–Ba, Sr–Ba: low confidence; use r_peak only as a lower estimate.**
  - These come from the CaCl₂ calibration (Ca–Cl–Ca ≈ 81°, face-sharing).
  - Crystal-like geometry gives much longer distances: crystalline SrCl₂ (fluorite) has Sr–Sr = a/√2 ≈ 4.9 Å, and BaCl₂ has Ba–Ba ≈ 5.2 Å.
  - The MgCl₂ calibration gives 4.6 / 4.9 Å.
  - **Plausible ranges: Sr–Sr 3.76–4.9 Å, Ba–Ba 4.0–5.2 Å.**
  - The r_cut values are permissive, so they do not block packing.
- **K–K and Cl–Cl in KCl:**
  - The neutron mean (4.77 / 4.76 Å) is used.
  - X-ray (4.30–4.40 Å) and all simulations in [P1] (4.20–4.61 Å) put them near 4.4–4.5 Å.
  - The 4.84 / 4.82 neutron row comes from the same source as the rejected CsCl rows (§3).
  - Treat these two values as 4.4–4.8 Å.

### 1c. Cation–cation sharing an O²⁻ (oxide environment)

| Pair | Environment | r_cut (Å) | r_peak (Å) | r_min1 (Å) | Grade (peak) |
|---|---|---|---|---|---|
| Li–Li | O-bridged | 1.65 | 2.29 | 3.32 | E |
| Li–Na | O-bridged | 1.82 | 2.53 | 3.67 | E |
| Li–K | O-bridged | 2.00 | 2.78 | 4.03 | E |
| Li–Rb | O-bridged | 2.06 | 2.87 | 4.16 | E |
| Li–Cs | O-bridged | 2.09 | 2.91 | 4.22 | E |
| Li–Sr | O-bridged | 2.11 | 2.93 | 4.25 | E |
| Li–Ba | O-bridged | 2.21 | 3.07 | 4.45 | E |
| Na–Na | O-bridged | 1.98 | 2.75 | 3.98 | E |
| Na–K | O-bridged | 2.15 | 2.98 | 4.32 | E |
| Na–Rb | O-bridged | 2.21 | 3.06 | 4.44 | E |
| Na–Cs | O-bridged | 2.23 | 3.10 | 4.50 | E |
| Na–Sr | O-bridged | 2.28 | 3.17 | 4.59 | E |
| Na–Ba | O-bridged | 2.37 | 3.30 | 4.78 | E |
| K–K | O-bridged | 2.30 | 3.20 | 4.64 | E |
| K–Rb | O-bridged | 2.36 | 3.27 | 4.75 | E |
| K–Cs | O-bridged | 2.38 | 3.31 | 4.80 | E |
| K–Sr | O-bridged | 2.46 | 3.42 | 4.96 | E |
| K–Ba | O-bridged | 2.55 | 3.54 | 5.13 | E |
| Rb–Rb | O-bridged | 2.41 | 3.35 | 4.86 | E |
| Rb–Cs | O-bridged | 2.44 | 3.38 | 4.91 | E |
| Rb–Sr | O-bridged | 2.53 | 3.51 | 5.09 | E |
| Rb–Ba | O-bridged | 2.61 | 3.62 | 5.25 | E |
| Cs–Cs | O-bridged | 2.46 | 3.42 | 4.96 | E |
| Cs–Sr | O-bridged | 2.55 | 3.55 | 5.15 | E |
| Cs–Ba | O-bridged | 2.64 | 3.66 | 5.31 | E |
| Sr–Sr | O-bridged | 2.60 | 3.61 | 5.23 | E |
| Sr–Ba | O-bridged | 2.69 | 3.74 | 5.43 | E |
| Ba–Ba | O-bridged | 2.79 | 3.87 | 5.62 | E |

These rows use the parent-crystal M–O–M geometry (§4.5). For example, Li–Li 2.29 Å matches a/2 = 2.31 Å in Li₂O.

**Rule for choosing between 1b and 1c:** use the O-bridged row when the two cations share an O²⁻ within the M–O r_min1. Otherwise use the Cl-bridged row. Applying the Cl-bridged cutoffs around O²⁻ would wrongly forbid real oxide coordination shells.

### 1d. Anion–anion (depends on which cation bridges them)

| Pair | Environment | r_cut (Å) | r_peak (Å) | r_min1 (Å) | Grade (peak) |
|---|---|---|---|---|---|
| Cl–Cl | around Li | 2.97 | 3.75 | 5.44 | M |
| O–O | around Li | 2.33 | 3.23 | 4.69 | E |
| Cl–O | around Li | 2.63 | 3.26 | 4.73 | E |
| Cl–Cl | around Na | 2.97 | 3.91 | 5.67 | M |
| O–O | around Na | 2.80 | 3.89 | 5.64 | E |
| Cl–O | around Na | 2.75 | 3.82 | 5.53 | E |
| Cl–Cl | around K | 3.43 | 4.76 | 6.90 | M |
| O–O | around K | 3.26 | 4.52 | 6.56 | E |
| Cl–O | around K | 3.13 | 4.35 | 6.30 | E |
| Cl–Cl | around Rb | 3.46 | 4.80 | 6.96 | M |
| O–O | around Rb | 3.41 | 4.74 | 6.87 | E |
| Cl–O | around Rb | 3.25 | 4.52 | 6.55 | E |
| Cl–Cl | around Cs | 3.49 | 4.85 | 7.03 | M |
| O–O | around Cs | 3.48 | 4.83 | 7.01 | E |
| Cl–O | around Cs | 3.40 | 4.72 | 6.85 | E |
| Cl–Cl | around Sr | 2.97 | 3.89 | 5.64 | E |
| O–O | around Sr | 2.60 | 3.61 | 5.23 | E |
| Cl–O | around Sr | 2.64 | 3.66 | 5.31 | E |
| Cl–Cl | around Ba | 2.99 | 4.16 | 6.03 | E |
| O–O | around Ba | 2.79 | 3.87 | 5.62 | E |
| Cl–O | around Ba | 2.83 | 3.93 | 5.69 | E |

For a mixed melt, interpolate the "around M" rows by the local cation mole fraction. Equivalently, evaluate §4.4 with composition-weighted distances.

---

## 2. First-hand data: cation–Cl⁻

### 2.1 Alkali chlorides

Compilation: Walz & van der Spoel [P1], ESI Table S3 (experiment; X = X-ray, N = neutron). Several rows in that table are themselves taken from the review by Marcus (2016), which is P1's ref. 10. The original papers are cited by number in [P1] and were not re-opened. Per the independent review, they include Edwards et al. 1975 (NaCl), Derrien & Dupuy 1975 (KCl, CsCl), Mitchell et al. 1976 (RbCl), Locke et al. 1985 (CsCl), McGreevy & Howe (LiCl), Levy et al. 1960 and Ohno et al. 1978 (X-ray).

| Salt | T (K) | r(M–Cl) Å | N(M–Cl) | Method |
|---|---|---|---|---|
| LiCl | 923 | 2.40 | 4 | X |
| LiCl | n/g | 2.47 | 4 | X |
| LiCl | 893 | 2.55 | 4.1 | X |
| LiCl | 958 | 2.30 | 5.5 | N |
| LiCl | n/g | 2.45 | 3.5 | N |
| NaCl | 1083 | 2.77 | 3.7 | X |
| NaCl | 1093 | 2.80/2.88 | 4.7 | X |
| NaCl | 1148 | 2.72 | 5.8 | N |
| NaCl | 1148 | 2.78 | 4.0 | N |
| KCl | n/g | 3.05 | 4.1 | X |
| KCl | n/g | 3.10 | 3.7 | X |
| KCl | 1083 | 3.10/3.20 | 5.2 | X |
| KCl | n/g | 3.14 | 5.8 | X |
| KCl | 1073 | 3.06 | 6.1 | N |
| KCl | n/g | 3.10 | 3.5 | N |
| RbCl | n/g | 3.30 | 4.2 | X |
| RbCl | 1023 | 3.18 | 6.9 | N |
| CsCl | n/g | 3.53 | 4.6 | X |
| CsCl | 973 | 3.38 | – | N |
| CsCl | 968 | 3.40 | 5.8 | N |
| CsCl | 973 | 3.38 | 3/5.9 | N |

**Selection rule for the recommended M–Cl peaks:** the arithmetic mean of the **neutron** rows.

| Pair | Recommended (Å) |
|---|---|
| Li–Cl | 2.38 |
| Na–Cl | 2.75 |
| K–Cl | 3.08 |
| Rb–Cl | 3.18 |
| Cs–Cl | 3.39 |

X-ray values run 0.02–0.15 Å longer for the heavy salts. A consensus review table [P22] gives Li–Cl 2.21, Na–Cl 2.70 and K–Cl 2.98 Å. It notes these are approximate and that LiCl coordination is debated. The Li–Cl peak therefore carries about ±0.1 Å uncertainty.

Polarizable-model (WBK) simulated peaks [P1, Table S1]:

| Salt | T (K) | r(M–Cl) Å |
|---|---|---|
| LiCl | 878 | 2.38 |
| NaCl | 1074 | 2.64 |
| KCl | 1043 | 2.96 |
| RbCl | 991 | 3.10 |
| CsCl | 918 | 3.28 |

For Na to Cs these are 0.1 Å shorter than the neutron values. For LiCl the simulations straddle the data: literature MD/MC gives 2.11–2.40 Å [P1, Table S4].

Other first-hand values:

- **Rb–Cl, XAFS** [P3]:
  - 3.20 Å, CN 5.3, pure molten RbCl at 1050 K;
  - 3.25 Å, CN 5.8, 15 % RbCl in the LiCl–KCl eutectic at 873 K.
- **First minimum (first-shell boundary r†)**, PIM at 1100 K [P2]:
  - LiCl 0.36 nm, NaCl 0.41 nm, KCl 0.45 nm. RIM gives the same r†. PIM average coordination is 4.34 / 5.16 / 5.18.
  - P2 defines r† as the first-shell boundary taken from the cation–anion g(r).
  - Independent simulation cross-check: BMH/Fumi–Tosi MD puts the Na–Cl first minimum at about 4.0 Å at 1066 K [P25].

### 2.2 Alkaline-earth chlorides

The compiled values come from X-ray and neutron diffraction [P22, Table 2; P23]:

| Salt | T (K) | r(M–Cl) Å | N(M–Cl) | Cl–Cl Å | M–M Å | Original paper |
|---|---|---|---|---|---|---|
| MgCl₂ | 998 | 2.42 | 4.3 | 3.56 | 3.81 | Biggin, Gay & Enderby, *J. Phys. C* 17, 977 (1984) |
| CaCl₂ | 1093 | 2.78 | 5.4 (5.3 in [P23]) | 3.73 | 3.6 | Biggin & Enderby, *J. Phys. C* 14, 3577 (1981) |
| **SrCl₂** | 1198 | **2.90** | 6.9 (6.6 in the text of [P22]) | – | – | McGreevy & Mitchell, *J. Phys. C* 15, 5537 (1982) |
| **BaCl₂** | 1298 | **3.10** | 7.7 (6.4 per Enderby 1987, via [P23]) | – | – | Edwards, Howe, Enderby & Page, *J. Phys. C* 11, 1053 (1978) |

Notes on this table:

- Cl–Cl and M–M for MgCl₂ and CaCl₂ come from [P23]. [P10] and [P8] give Mg–Mg 3.81 and Ca–Ca 3.73 Å respectively.
- The original *J. Phys. C* papers were not opened. The values are as compiled in [P22] and [P23].
- **Sr–Sr, Ba–Ba and Cl–Cl in SrCl₂/BaCl₂ were not found.** McGreevy & Mitchell (1982) determined the SrCl₂ partials, so the paper should contain them.

Further alkaline-earth data:

- **Sr²⁺ dilute in NaCl, EXAFS** [P5]: in 5 mol % SrCl₂–NaCl, N(Cl) = 7.5 ± 1.3 at 25 °C, falling to 5.8 ± 1.1 at 1000 °C. The distances are in that paper's SI Table S1, which was not accessible.
- **AIMD cutoffs** [P6], UClₙ–MgCl₂ at 1200 K:
  - Mg–Cl peak 2.40 Å with a CN cutoff of 3.4 Å. The cutoff is not stated to be the first minimum.
  - U(III)–Cl peak 2.74 Å with the cutoff at the RDF minimum, 4.0 Å.
- **Cross-check of the radius-sum model** (§4.1): it gives Sr–Cl 2.92 and Ba–Cl 3.09 Å, within 0.02 Å of the measured values.

## 3. First-hand data: like-ion pairs in pure alkali chloride melts

From [P1] Table S3 (experiment). Neutron entries give M–M and Cl–Cl separately.

| Salt | T (K) | r(M–M) Å | r(Cl–Cl) Å | Method |
|---|---|---|---|---|
| LiCl | 958 | 3.70 | 3.70 | N |
| LiCl | n/g | 3.80 | 3.80 | N |
| LiCl | 923 / 893 / n/g | 3.86 / 3.90 / 3.85 | – | X (mainly Cl–Cl) |
| NaCl | 1148 | 3.89 | 3.90 | N |
| NaCl | 1148 | 3.96 | 3.91 | N |
| NaCl | 1083 / 1093 | 4.15 / 4.20 | – | X |
| KCl | 1073 | 4.84 | 4.82 | N |
| KCl | n/g | 4.70 | 4.70 | N |
| KCl | 1083 / n/g | 4.30–4.40 | – | X (mainly Cl–Cl) |
| RbCl | 1023 | 4.86 | 4.80 | N |
| CsCl | 968 | 4.95 | 4.85 | N |
| CsCl | 973 (two rows) | 3.85 | 3.85 | N (rejected) |

The rejected CsCl rows are physically implausible, since they are shorter than √2 × r(Cs–Cl) = 4.79 Å. One of them comes from the same source (Derrien & Dupuy) as the KCl 4.84/4.82 row.

**Recommended values are the neutron means:**

| Pair | Recommended (Å) |
|---|---|
| Li–Li, Cl–Cl in LiCl | 3.75, 3.75 |
| Na–Na, Cl–Cl in NaCl | 3.93, 3.91 |
| K–K, Cl–Cl in KCl | 4.77, 4.76 |
| Rb–Rb, Cl–Cl in RbCl | 4.86, 4.80 |
| Cs–Cs, Cl–Cl in CsCl | 4.95, 4.85 |

Simulated ranges [P1, Table S4]:

| Salt | Like-pair range (Å) |
|---|---|
| LiCl | 3.65–4.10 |
| NaCl | 3.90–4.25 |
| KCl | 4.20–4.80 |
| RbCl | 4.50–4.80 |
| CsCl | 4.60–4.90 |

[P10] gives Na–Na 3.90 Å (experiment), 4.05 (PIM) and 4.10 (RIM).

**Unlike cation pairs in mixtures: no numerical first-hand value was readable.** The relevant studies are:

- MD of LiCl–RbCl, LiCl–CsCl, NaCl–RbCl and NaCl–CsCl at 1100 K [P11], paywalled;
- FPMD of LiCl–KCl, with figures only [P12];
- a neutron study of Li–Cl in molten (Li,K)Cl and (Li,Cs)Cl [P26], not readable.

---

## 4. Estimation model (grade-E values)

The constants below are printed when the script runs.

### 4.1 Radius-sum model (now a cross-check only)

r(M–Cl) ≈ r_Shannon,VI(M) + r_Shannon,VI(Cl) + δ.

| Cation | δ (Å) |
|---|---|
| Na | −0.08 |
| K | −0.11 |
| Rb | −0.15 |
| Cs | −0.09 |
| **Alkali mean** | **−0.1075** |
| Mg | −0.11 |
| Ca | −0.03 |
| **Divalent mean** | **−0.07** |

The divalent mean predicts Sr–Cl 2.92 and Ba–Cl 3.09 Å, against measured 2.90 and 3.10. The alkali mean offset is used only for **Cs–O**: 1.67 + 1.40 − 0.1075 = 2.96 Å.

Shannon radii [P7], CN VI unless stated:

| Ion | Radius (Å) |
|---|---|
| Li⁺ | 0.76 (CN IV: 0.59) |
| Na⁺ | 1.02 |
| K⁺ | 1.38 |
| Rb⁺ | 1.52 |
| Cs⁺ | 1.67 |
| Mg²⁺ | 0.72 |
| Ca²⁺ | 1.00 |
| Sr²⁺ | 1.18 |
| Ba²⁺ | 1.35 |
| Cl⁻ | 1.81 |
| O²⁻ | 1.40 (CN IV: 1.38; CN II: 1.35) |

### 4.2 First minimum

- **Rb–Cl and Cs–Cl.** The ratio r†/r_peak is formed inside the polarizable-model frame:
  - mean(3.6/2.38, 4.1/2.64, 4.5/2.96) = **1.529**, using P2's r† and P1's WBK peaks.
  - This ratio is applied to the WBK peaks of RbCl (3.10 Å) and CsCl (3.28 Å), giving r_min1 = **4.74 and 5.01 Å**.
  - The result continues smoothly from K–Cl (4.5 Å).
  - The measured peak is not used here because r† and the measured peak come from different sources.
- **Sr–Cl, Ba–Cl and all M–O pairs:** ratio **1.438**, the mean of the AIMD cutoff ratios Mg–Cl 3.4/2.40 and U(III)–Cl 4.0/2.74 [P6].
  - This is used for cations with charge ≥ 2 and, by assumption, for pairs with the doubly charged O²⁻.
  - The Mg cutoff is not stated to be the first minimum.
- **Like-sign pairs:** ratio 1.45. **This is an assumption; no source was found.**

### 4.3 Onset r_cut (all grade E)

- **Opposite-charge pairs:** r_cut = 0.82 × r_peak.
- **Same-charge pairs:** r_cut = 0.72 × r_peak.
- **Anion–anion pairs:** additionally r_cut ≥ 0.82 × (contact sum of Shannon radii), i.e. Cl–Cl ≥ 2.97, O–O ≥ 2.30, Cl–O ≥ 2.63 Å.
- **To be replaced by values read from the NaCl partials [P4]** (§6).

**Packing test (revision 2.1).** Random sequential insertion was run with these cutoffs at the experimental melt densities:

- about 300 ions per cubic periodic box (100 formula units for SrCl₂ and BaCl₂, 150 for the alkali chlorides), each placed cation then anion;
- up to 2 × 10⁵ attempts per ion;
- two random seeds.

Pure LiCl, NaCl, KCl, RbCl, CsCl, SrCl₂ and BaCl₂ all placed 100 % of their ions.

Revision 2 used 0.80 × peak for like pairs and 0.92 × contact for anions. With those, LiCl and SrCl₂ jammed at 83 % and KCl at 86 %.

LiCl, SrCl₂ and BaCl₂ still reach 97–100 % only near the floor. If a composition jams:
- insert at about 90 % density and then compress; or
- relax with Monte Carlo using a soft repulsion.

Do not raise r_cut without re-running a packing test.

### 4.4 Same-sign pairs bridged by a common counter-ion

r_XY = F × √(r_XA² + r_YA²). F = 1 corresponds to a 90° X–A–Y angle.

- **Alkali chlorides: F = 1.050.** Fitted to the measured like pairs of NaCl, KCl, RbCl and CsCl (range 1.005–1.095). LiCl is excluded because its Cl–Cl distance is anion-contact limited.
- **Divalent cations: calibrated on molten CaCl₂.** Its coordination (N ≈ 5.4) is closest to that of Sr/Ba.
  - F(cation–cation) = Ca–Ca / (√2 × Ca–Cl) = 3.6 / (√2 × 2.78) = **0.916**.
  - F(anion–anion around M²⁺) = 3.73 / (√2 × 2.78) = **0.949**.
  - For comparison, MgCl₂ (a network-forming melt, N ≈ 4.3) gives 1.113 and 1.040. These are 22 % (cation pairs) and 10 % (anion pairs) larger.
  - This calibration produces Sr–Sr 3.76, Ba–Ba 4.01, and Cl–Cl around Sr/Ba 3.89 / 4.16 Å.
  - The Cl–Cl values are ±0.3 Å.
  - Sr–Sr and Ba–Ba are **low confidence**. CaCl₂'s Ca–Cl–Ca ≈ 81° implies face-sharing, which may not transfer to CN 7–8 cations. The plausible ranges are given in §1b.
- **Mixed alkali–divalent cation pairs:** F = mean of the two class values.
- **Anion–anion floor:** these pairs are raised to at least the contact sum + 0.05 Å.

### 4.5 Free oxide ion O²⁻

**No accessible first-hand RDF of O²⁻ in a chloride melt was found** ([P13]–[P15] exist but were not readable).

r(M–O) is the parent-crystal nearest-neighbour distance minus 0.02 Å:

| Cation | Parent crystal | Crystal M–O (Å) | Source |
|---|---|---|---|
| Li | antifluorite M₂O, a√3/4 | 1.996 | [P9] |
| Na | antifluorite | 2.403 | [P9] |
| K | antifluorite | 2.792 | [P9] |
| Rb | antifluorite | 2.919 | [P9] |
| Sr | rock-salt MO | 2.57 | [P24] |
| Ba | rock-salt MO | 2.76 | [P24] |

These agree with the Shannon sums 2.58 and 2.75 Å. Cs–O uses the radius sum of §4.1.

**O–O** uses the O–M–O angle of the parent crystal. That gives 1.633 × r(M–O) for the antifluorite alkalis, i.e. the crystal O–O distance a/√2. It gives √2 × r(M–O) for rock-salt Sr/Ba.

**Cl–O** uses the class F of §4.4.

**Cation pairs sharing an O²⁻ (§1c)** use the parent-crystal M–O–M angle: r_XY = F_O × √(r_XO² + r_YO²).

- **Antifluorite alkalis:** F_O = (2/√3)/√2 = 0.816. The M–O–M angle is 70.5°, which gives M–M = a/2.
- **Rock-salt Sr/Ba:** F_O = 1, a 90° angle.
- **Mixed pairs:** the mean of the two class values.
- **Cs:** treated as antifluorite. Cs₂O is actually anti-CdCl₂, so this is an assumption.

Loose cross-checks from molten chlorates [P16], where the O belongs to ClO₃⁻ rather than being free:

- Na–O 2.412 Å at 560 K
- K–O 2.856 Å at 660 K

Speciation caveat: in chloride melts O²⁻ forms oxo-complexes with small or highly charged cations. Treat all M–O values as ±0.1 Å.

---

## 5. Summary of data availability

| Status | Pairs |
|---|---|
| Measured peak (M) | Li/Na/K/Rb/Cs/Sr/Ba–Cl; Li–Li, Na–Na, K–K, Rb–Rb, Cs–Cs; Cl–Cl in the five alkali chlorides |
| Simulated (S) | r_min1 of Li–Cl, Na–Cl, K–Cl |
| **Unavailable → estimated (E)** | every pair with O²⁻, including all O-bridged cation pairs; all unlike cation–cation pairs; Sr–Sr, Ba–Ba, Sr–Ba (low confidence); Cl–Cl around Sr/Ba; every r_cut; r_min1 except Li/Na/K–Cl |

## 6. Recommended follow-ups for the coder

1. **Calibrate r_cut and r_min1 from [P4].** Download `Fig6_rdf_pim+expt_forpaper_v4.agr` (https://researchdata.bath.ac.uk/1165/8/Fig6_rdf_pim%2Bexpt_forpaper_v4.agr, CC BY 4.0). Read the onset, peak and first minimum of g_NaNa, g_NaCl and g_ClCl. Replace the 0.82 / 0.72 / 1.45 ratios and the 0.82 floor in the script. Then re-run a packing test (§4.3). This sandbox could not reach that host.
2. **Divalent geometry.** McGreevy & Mitchell (1982, SrCl₂) and Edwards et al. (1978, BaCl₂) contain the partial g(r) needed to replace the CaCl₂-based Sr–Sr, Ba–Ba and Cl–Cl estimates. The PIM data for MgCl₂–BaCl₂ at 1273 K in [P21] (Zenodo, doi:10.5281/zenodo.20560046) would also serve.
3. **Mixed melts.** Compute anion–anion and cation–cation limits from the local cation mix (§4.4).
4. **Confinement.** Inside the LTA α-cage, use framework-oxygen values for cation–framework contacts, not this table. Kim et al. (2008), ref. R27 of the occlusion evidence file, report Ag–O(framework) = 2.493 Å and K⁺–Cl⁻ through a 6-ring = 2.70 Å.

---

## 7. References

| Tag | Reference | Used for | Access |
|---|---|---|---|
| P1 | M.-M. Walz, D. van der Spoel, "Molten alkali halides – temperature dependence of structure, dynamics and thermodynamics," *Phys. Chem. Chem. Phys.* **21**(34), 18516–18524 (2019), doi:10.1039/C9CP03603B. ESI: https://www.rsc.org/suppdata/c9/cp/c9cp03603b/c9cp03603b1.pdf | Alkali chloride peaks/CNs (Tables S1, S3, S4) | ESI read |
| P2 | S. Roy, F. Wu, H. Wang, A. S. Ivanov, S. Sharma, … V. S. Bryantsev (18 authors), "Structure and dynamics of the molten alkali-chloride salts from an X-ray, simulation, and rate theory perspective," *Phys. Chem. Chem. Phys.* **22**(40), 22900–22917 (2020), doi:10.1039/D0CP03672B | r† (LiCl, NaCl, KCl) | Author version read |
| P3 | Y. Okamoto, T. Yaita, H. Shiwaku, S. Suzuki, "Behavior of monovalent ion in LiCl-KCl eutectic melt," Photon Factory Activity Report 2007 #25 Part B, p. 14 (2008). https://pfwww.kek.jp/acr2007pdf/part_b/pf07b014.pdf | Rb–Cl XAFS | Read |
| P4 | A. Zeidler, P. S. Salmon, T. Usuki, S. Kohara, H. E. Fischer, M. Wilson, "Structure of molten NaCl and the decay of the pair-correlations," *J. Chem. Phys.* **157**, 094504 (2022), doi:10.1063/5.0107620; data doi:10.15125/BATH-01165 | Open NaCl partials for calibration | Record read; file not reachable |
| P5 | T. J. Lynch, N. R. Birkner, M. S. Christian, J. A. Wrubel, J. Schorne-Pinto, A. Van Veelen, J. R. Bargar, T. M. Besmann, K. S. Brinkman, W. K. S. Chiu, "In situ determination of speciation and local structure of NaCl–SrCl₂ and LiF–ZrF₄ molten salts," *J. Phys. Chem. B* **126**(7), 1539–1550 (2022), doi:10.1021/acs.jpcb.1c07552 | Sr–Cl EXAFS coordination | Main text read |
| P6 | B. Li, S. Dai, D.-e. Jiang, "First principles molecular dynamics simulations of UClₙ–MgCl₂ (n = 3, 4) molten salts," *Phys. Chem. Chem. Phys.* **24**(39), 24281–24289 (2022), doi:10.1039/D2CP02417A | Mg–Cl, U–Cl peaks and cutoffs | Author version read |
| P7 | R. D. Shannon, *Acta Cryst.* **A32**, 751–767 (1976); values via https://www.knowledgedoor.com/2/elements_handbook/shannon-prewitt_effective_ionic_radius.html (parts 1–4) | Ionic radii | Read |
| P8 | D. Chen, X. Wei, J. Ding, S. Liu, W. Wang, J. Lu, "Molecular dynamics simulation of molten CaCl₂ and NaCl–CaCl₂ eutectic molten salt," *Energy Proceedings* 21 (ICAE 2021) | Ca–Ca 3.73 (its citation of Biggin is garbled) | Read |
| P9 | Wikipedia, "Fluorite structure," antifluorite lattice constants at 300 K (Li₂O 4.61, Na₂O 5.55, K₂O 6.449, Rb₂O 6.74 Å). Unsourced there; primary source is Zintl, Harder & Dauth, *Z. Elektrochem.* 40, 588 (1934), not opened. Other compilations differ by ≤ 0.01 Å | M–O (alkali) | Read (secondary) |
| P10 | J. Lu, S. Yang, G. Pan, J. Ding, S. Liu, W. Wang, "Thermal and transport properties of molten chloride salts with polarization effect on microstructure," *Energies* **14**(3), 746 (2021), doi:10.3390/en14030746 | Na–Na, Mg–Mg | Read |
| P11 | J. Wang, J. Wu, G. Lu, J. Yu, "Molecular dynamics study of the transport properties and local structures of molten alkali metal chlorides. Part III," *J. Mol. Liq.* **238** (2017), doi:10.1016/j.molliq.2017.03.103 *(pages unconfirmed)* | Unlike cation pairs | Abstract only |
| P12 | A. Bengtson, H. O. Nam, S. Saha, R. Sakidja, D. Morgan, *Comput. Mater. Sci.* **83**, 362–370 (2014), doi:10.1016/j.commatsci.2013.10.043 | LiCl–KCl FPMD (figures only) | Draft read |
| P13 | L. Lin, J. Song, C. Zhu, S. Guo, W. Zhou, Y. Wang, *J. Mol. Liq.* **413**, 125950 (2024), doi:10.1016/j.molliq.2024.125950 | Oxide in LiCl | Not read |
| P14 | "First-principles molecular dynamics simulation and conductivity measurements of a molten xLi₂O–(1−x)LiCl system," *J. Phys. Chem. B* **117**(18), 5668 (2013) *(authors/DOI unconfirmed)* | Oxide in LiCl | Not read |
| P15 | "Molecular dynamics study of local structure and migration properties of LiCl–Li₂O–Li molten salts based on machine-learned deep potential," *Nucl. Eng. Des.* (2025), Elsevier PII S0029549325002298 *(details unconfirmed)* | Oxide in LiCl | Not read |
| P16 | Y. Iwadate, H. Ushida, K. Fukushima, M. Misawa, T. Fukunaga, *J. Phys. Soc. Jpn.* **70** Suppl. A, 353–355 (2001) | Na–O, K–O in chlorate melts | Read |
| P17 | X. Zhang, J. Song, L. Lin, Y. Ding, X. Xu, W. Zhou, S. Guo, D. Wang, "Research on the transport and thermodynamic properties of Sr²⁺ and Cs⁺ in LiCl–KCl using deep potential molecular dynamics," *J. Radioanal. Nucl. Chem.* (2024), doi:10.1007/s10967-024-09757-3 | Sr/Cs in eutectic | Not read (rate-limited) |
| P18 | "Molecular dynamics simulation of molten strontium chloride based on deep potential," *J. Mol. Liq.* (2021), PII S0167732221031056 *(details unconfirmed)* | SrCl₂ melt | Not read |
| P19 | J. Zhao, T. Feng, G. Lu, "Deep learning potential assisted prediction of local structure and thermophysical properties of the SrCl₂–KCl–MgCl₂ melt," *J. Chem. Theory Comput.* **20**(17), 7611–7623 (2024) | Sr in a mixed melt | Not read |
| P20 | Y. D. Fomin, "Structure and dynamics of molten SrCl₂," *Phys. Chem. Liq.* **58**(5), 693–700 (2020), doi:10.1080/00319104.2020.1757095 | SrCl₂ MD (RDFs in figures only) | Record read |
| P21 | Y. Chen, M. S. Emerson, H. H. Nguyen, R. Ogbodo, V. S. Bryantsev, J. F. Wishart, C. J. Margulis, *J. Phys. Chem. B* **130**(25), 6376–6382 (2026), doi:10.1021/acs.jpcb.6c01313; data doi:10.5281/zenodo.20560046 | PIM g(r) of Mg/Ca/Sr/BaCl₂ at 1273 K (SI) | Main text read |
| P22 | A.-L. Rollet, M. Salanne, "Studies of the local structures of molten metal halides," *Annu. Rep. Prog. Chem., Sect. C* **107**, 88–123 (2011), doi:10.1039/C1PC90003J | Sr–Cl, Ba–Cl, Ca–Cl, Mg–Cl (Table 2); alkali consensus values (Table 1) | Tables read |
| P23 | D. K. Belashchenko, O. I. Ostrovski, "Computer modelling of molten halides using diffraction data," *VII Int. Conf. on Molten Slags, Fluxes and Salts* (SAIMM, 2004), pp. 503–508. https://www.pyrometallurgy.co.za/MoltenSlags2004/503-Belashchenkko.pdf | Mg/Ca Cl–Cl and M–M; original-paper citations; Enderby 1987 CNs | Read |
| P24 | UCL crystal-data pages for SrO and BaO (nearest-neighbour distance 2.57 and 2.76 Å; unsourced there). https://www.homepages.ucl.ac.uk/~ucapahh/research/crystal/sro.htm, …/bao.htm | Sr–O, Ba–O | Read (secondary) |
| P25 | arXiv:cond-mat/0512282, "Physics of solid and liquid alkali halide surfaces near the melting point" | Na–Cl first minimum about 4.0 Å (MD) | Per independent review |
| P26 | Neutron study of the Li–Cl distance in molten (Li,K)Cl and (Li,Cs)Cl, *Mol. Phys.* **82**(5) (1994), doi:10.1080/00268979400100644 *(authors unconfirmed)* | Li–Cl in mixtures | Not read |

Original papers named in [P22]/[P23] but not opened:

- Biggin, Gay & Enderby, *J. Phys. C* 17, 977–985 (1984)
- Biggin & Enderby, *J. Phys. C* 14, 3577–3583 (1981)
- McGreevy & Mitchell, *J. Phys. C* 15, 5537–5550 (1982)
- Edwards, Howe, Enderby & Page, *J. Phys. C* 11, 1053–1057 (1978)
- J. E. Enderby, in *Molten Salt Chemistry* (Reidel, 1987), pp. 1–15

---

## 8. Changes in revision 2 (from the independent review)

1. **Sr–Cl and Ba–Cl are now measured.** Sr–Cl is 2.90 Å and Ba–Cl is 3.10 Å [P22, P23], replacing the model estimates 2.92 and 3.09. The model is kept as a cross-check, and it agrees within 0.02 Å.
2. **Divalent bridging recalibrated on CaCl₂.** The previous value came from MgCl₂ alone (F = 1.113), which gave Cl–Cl around Sr/Ba of about 4.3/4.5 Å and r_cut of 3.40/3.60 Å. Those were too long and would have rejected real Cl–Cl contacts. They are now 3.89/4.16 Å (r_cut 3.33 Å in revision 2; 2.97 Å in revision 2.1). Sr–Sr and Ba–Ba changed from 4.60/4.86 to 3.76/4.01 Å.
3. **O–O now uses the parent-crystal bridge angle.** For the alkali oxides this is tetrahedral, which reproduces the antifluorite O–O distance a/√2. Previously O–O was 0.37–0.54 Å too short.
4. **A single, stated selection rule for measured peaks:** the neutron mean.
   - Li–Cl changed from 2.40 to 2.38, Na–Cl from 2.76 to 2.75, K–Cl from 3.07 to 3.08, Rb–Cl from 3.19 to 3.18 Å.
   - Like pairs: Li 3.70 → 3.75, Na 3.90 → 3.93/3.91, K 4.50/4.45 → 4.77/4.76 Å.
   - The K values were a compromise graded M; they are now the actual neutron mean, with the X-ray disagreement flagged.
   - Alkali F changed from 1.031 to 1.050.
5. **Rb–Cl and Cs–Cl r_min1** are now computed consistently in the polarizable-model frame: 4.89 → 4.74 and 5.20 → 5.01 Å.
6. **Sr–O and Ba–O** now come from rock-salt crystal distances (2.51 → 2.55 and 2.68 → 2.74 Å). The Cs–O offset is stated.
7. **Smaller corrections:**
   - U(III) is triply charged, not doubly.
   - The Mg–Cl cutoff is not stated to be a minimum.
   - The like-pair temperature shift is +0.06 Å per 100 K, not 0.08.
   - "Simulations underestimate" does not hold for LiCl.
   - The crystal distances K 2.792 and Rb 2.919 Å are now given to the correct last digit.
   - Omitted P1 rows have been added.
   - The script comments now cite the correct sources.
   - The Ag–O value is cited to its true source.
8. **Bibliographic details completed** for P1, P2, P5, P6, P7, P10, P12, P13, P17, P19 and P21. P22–P26 are added.

**Revision 2.1 (blocking-error review):**

1. **r_cut loosened so random insertion can reach experimental density.** The like-pair onset changed from 0.80 to 0.72 × peak, and the anion contact floor from 0.92 to 0.82 × contact (Cl–Cl 3.33 → 2.97 Å). A packing test now places 100 % of ions for all seven pure melts (§4.3). The test was run by the reviewer and re-run here with two seeds.
2. **Added 28 O-bridged cation–cation rows** (§1c), with a rule for choosing between Cl-bridged and O-bridged rows. Without them the Cl-bridged cutoffs would forbid real oxide coordination shells.
3. **Sr–Sr, Ba–Ba and Sr–Ba flagged low confidence**, with plausible ranges of 3.76–4.9 and 4.0–5.2 Å. KCl like pairs are flagged as 4.4–4.8 Å.
4. The MgCl₂ comparison is now quantified (22 % / 10 %), and the history numbers in item 2 above are corrected.
5. The re-check confirmed that all M- and S-graded values match their sources. The CSV, tables and script are consistent.
