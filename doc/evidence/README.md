# LTA salt-occlusion evidence package (revision 3, 2026-10-17)

| File | Purpose |
|---|---|
| `lta_salt_occlusion_evidence.md` | Literature evidence for chloride-salt occlusion in zeolite A: calibration targets (≈10.5 ± 1 Cl per pseudo-cell, LiCl–KCl, 823 K), cation effects, framework stability, Cl siting, unit conversions, model recommendations, references R1–R34 (unchanged from rev 2.1) |
| `molten_salt_pair_distances.md` | Pair distances (r_cut, r_peak, r_min1) for Li⁺, Na⁺, K⁺, Rb⁺, Cs⁺, Sr²⁺, Ba²⁺, Cl⁻, O²⁻, with sources, grades (M/S/E), estimation model and references P1–P26 (unchanged from rev 2.1; see corrections below) |
| `molten_salt_pair_distances.csv` | Machine-readable pair table (91 rows) for the random molten-salt structure generator |
| `molten_salt_pair_table.py` | Regenerates the pair CSV (`python3 molten_salt_pair_table.py`; overwrites the CSV in the working directory) |
| **`lta_coordination_numbers.md`** | **New.** Coordination numbers of each ion in pure molten chlorides, LiCl–KCl and dilute melts, with dissolved oxide, and inside LTA; recommended values with selection rules; references Q13–Q28 |
| **`lta_coordination_numbers.csv`** | **New.** One row per value, 215 rows (columns `centre, neighbour, environment, T_K, composition, CN, CN_uncertainty, cutoff_A, method, grade, source_tag, note`); `grade` = source grade + M/S/E kind; includes the E working values and an "unavailable" row for each missing pair |
| **`lta_coordination_numbers_table.py`** | **New.** Regenerates the CN CSV, including the derived (D) anion-centred rows (`python3 lta_coordination_numbers_table.py OUTFILE`) |
| **`lta_salt_loading_lattice_contraction.md`** | **New.** Lattice parameters of hydrated, dehydrated and halide-loaded zeolite A; contraction factors f_a, f_V; bound and working value for 4A + LiCl–KCl; references Q1–Q12 |
| **`lta_salt_loading_lattice.csv`** | **New.** Lattice-parameter table (row numbers match the Markdown) |

Read §0 of each Markdown file first (conventions and how to use the numbers). Each file's change log lists what the three independent reviews corrected.

## Headline results of revision 3

- **Coordination numbers.** Recommended first-shell N(M–Cl) for pure alkali chlorides, integrated to the first-shell boundary of the simulated g(r) (for Li/Na/K this is the pair file's `r_min1`): LiCl 4.3, NaCl 5.2, KCl 5.2, RbCl 5.4, CsCl 5.5 (±0.5). Alkaline-earth chlorides (neutron, integration rule not stated): MgCl₂ 4.3, CaCl₂ 5.3–5.4, SrCl₂ 6.9; **BaCl₂ unresolved (6.4 vs 7.7; working 7 ± 1, E)**. LiCl–KCl eutectic at 773 K (AIMD, R_min not stated): Li–Cl ≈ 4.2, K–Cl ≈ 7.0 (definition-dependent). Dilute Sr²⁺: 6 ± 1 Cl (E). No numeric Li–O or other M–O CN was readable except Sr–O = 1 (SrO in LiCl–KCl).
- **Lattice contraction.** No numeric lattice parameter for salt-loaded 4A was readable. 4A contracts on LiCl–KCl occlusion, partly by growth of a second, smaller-a LTA phase [R18]. **f_V for 4A + LiCl–KCl is unavailable**; bound 0.919 ≤ f_V < 1 (from a figure's axis range, with stated assumptions) and provisional working value f_V = 0.97 ± 0.03 (grade E, from a K/Br-loaded LTA analogue).

## Corrections to rev 2.1 files found in this revision (not yet applied to those files)

| File / tag | Correction |
|---|---|
| `molten_salt_pair_distances.md` P5 | Lynch et al. SI Table S1 was read once in this session, appended to the Clemson-hosted PDF (Sr–Cl R = 2.92–3.03 Å, N at 25–1000 °C); the reviewer's re-fetch returned only the main text, so values not in the main text rest on that single reading. The rev 2.1 statement that it "was not accessible" is superseded |
| `molten_salt_pair_distances.md` P14 | **Mis-cited.** *J. Phys. Chem. B* 117(18), 5668 (2013) is Ohkubo et al. on molten xLi₂O–(1−x)**B₂O₃**, not LiCl. Drop P14 |
| P13, P15, P17, P19, P11 | Bibliographic details completed (see `lta_coordination_numbers.md` §7); P19 DOI is 10.1021/acs.jctc.4c00824 |
| Hand-off note for P1 | P1 ESI Tables S1–S4 are pure alkali halides, not LiCl–CsCl mixtures |
| Edwards et al. 1978 (BaCl₂) | Correct DOI 10.1088/0022-3719/11/6/012 |
| `lta_salt_occlusion_evidence.md` R30 | a = 12.061(1) Å now verified (abstract); bibliography *J. Phys. Chem. C* 122(12), 6661–6668 (2018) |
| `lta_salt_occlusion_evidence.md` R27 | Bibliography *J. Phys. Chem. C* 112(30), 11181–11193 confirmed |

## Full texts most worth obtaining (library access needed)

1. Lin et al. 2024 (*J. Mol. Liq.* 413, 125950) and Xu et al. 2025 (*Nucl. Eng. Des.* 438, 114052): Li–O coordination in LiCl–Li₂O.
2. Lexa, Leibowitz & Kropf, *J. Nucl. Mater.* 279, 57 (2000), and Harward et al., *J. Nucl. Mater.* (2023): numeric lattice parameters of salt-loaded 4A.
3. Edwards PhD thesis (Leicester 1977, figshare) or Rollet & Salanne 2011 Table 2: the BaCl₂ CN conflict.
4. Richardson, *MRS Proc.* 465 (neutron Rietveld of salt-occluded 4A/5A).
