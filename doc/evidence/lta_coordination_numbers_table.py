"""Build lta_coordination_numbers.csv from values read in session 2026-10-17.

Every non-D row is transcribed from a source (see lta_coordination_numbers.md).
D rows are computed here: N(Cl-M) = N(M-Cl) * x_M / x_Cl  (number conservation).
Usage: python3 lta_coordination_numbers_table.py OUTFILE

grade column: source grade A-D (E for estimates), then the kind of value after a comma:
M = measured, S = simulated/modelled, E = estimate (D rows carry no kind). "-" = no value.
A bare "C" (no kind) = method unknown (only the Q27 "eightfold" row).
cutoff_A holds an integration cutoff or boundary only; fitted bond distances are in note.
"""
import csv, sys

COLS = ["centre", "neighbour", "environment", "T_K", "composition", "CN",
        "CN_uncertainty", "cutoff_A", "method", "grade", "source_tag", "note"]
rows = []


SIM_KEYS = ("MD", "MC", "AIMD", "USPEX", "DFT", "PIM", "RIM", "BELION", "simulation")
MEAS_KEYS = ("X-ray", "neutron", "EXAFS", "XAFS", "experiment", "single-crystal", "diffraction", "Raman")


def kind(meth, grade):
    if grade in ("D", "-", ""):
        return ""
    if grade.startswith("E"):
        return "E"
    if any(k in meth for k in SIM_KEYS):
        return "S"
    if any(k in meth for k in MEAS_KEYS):
        return "M"
    return ""


def add(c, n, env, T, comp, cn, unc, cut, meth, grade, tag, note=""):
    if cn == "unavailable":
        grade = "-"
    k = kind(meth, grade)
    g = f"{grade},{k}" if k and grade != "E" else grade
    rows.append(dict(zip(COLS, [c, n, env, T, comp, cn, unc, cut, meth, g, tag, note])))


# ---- 1. Pure alkali chlorides: P1 Table S3 (experiment) ----
P1_EXP = [  # salt, T, N+-, N++, N--, method, original
    ("LiCl", "923", "4", "12", "", "X-ray", "Ohno et al. 1978; like-pair is mainly Cl-Cl"),
    ("LiCl", "", "4", "12", "", "X-ray", "Levy et al. 1960"),
    ("LiCl", "893", "4.1", "9", "", "X-ray", "Zarzycki 1958"),
    ("LiCl", "958", "5.5", "14", "14", "neutron", "Howe & McGreevy 1988 / McGreevy & Howe 1989"),
    ("LiCl", "", "3.5", "8.3", "8.3", "neutron", "Levy et al. 1960"),
    ("NaCl", "1083", "3.7", "", "", "X-ray", "Ohno & Furukawa 1981"),
    ("NaCl", "1093", "4.7", "9", "", "X-ray", "Zarzycki 1958"),
    ("NaCl", "1148", "5.8", "13", "13.3", "neutron", "Edwards et al. 1975"),
    ("NaCl", "1148", "4.0", "", "", "neutron", "Edwards et al. 1975 data as reinterpreted by Biggin & Enderby 1982"),
    ("KCl", "", "4.1", "", "", "X-ray", "Takagi et al. 1979"),
    ("KCl", "", "3.7", "", "", "X-ray", "Levy et al. 1960"),
    ("KCl", "1083", "5.2", "10", "", "X-ray", "Zarzycki 1958"),
    ("KCl", "", "5.8", "", "", "X-ray", "Lark & Miller 1937"),
    ("KCl", "1073", "6.1", "12/15.9", "12.3/16.2", "neutron", "Derrien & Dupuy 1975 (two like-pair values given)"),
    ("KCl", "", "3.5", "7.4", "7.4", "neutron", "Levy et al. 1960"),
    ("RbCl", "", "4.2", "", "", "X-ray", "Levy & Danford 1964"),
    ("RbCl", "1023", "6.9", "13", "14", "neutron", "Mitchell, Poncet & Stewart 1976"),
    ("CsCl", "", "4.6", "7.1", "", "X-ray", "Levy et al. 1960"),
    ("CsCl", "968", "5.8", "15.4", "16.3", "neutron", "Locke et al. 1985"),
    ("CsCl", "973", "3/5.9", "", "", "neutron", "Derrien & Dupuy 1975 (two values given)"),
]
for salt, T, npm, npp, nmm, meth, orig in P1_EXP:
    m = salt[:-2]
    note = f"original: {orig}; CN definition varies by study (P1 ESI)"
    add(m, "Cl", "pure melt", T, salt, npm, "", "", meth, "C", "P1 Table S3", note)
    if npp:
        lab = "Cl" if (meth == "X-ray") else m
        add(lab, lab, "pure melt", T, salt, npp, "", "", meth, "C", "P1 Table S3",
            note + ("; X-ray like-pair reflects mainly Cl-Cl" if meth == "X-ray" else ""))
    if nmm:
        add("Cl", "Cl", "pure melt", T, salt, nmm, "", "", meth, "C", "P1 Table S3", note)

# ---- P1 Table S1 (WBK polarisable) and S2 (JC non-polarisable): CN to first min, +/-0.5 ----
P1_SIM = {
    "WBK polarisable FF MD": [
        ("LiCl", "878.15", "5.2", "14.7", "14.8"), ("LiCl", "978.15", "5.0", "14.1", "14.1"),
        ("NaCl", "1074.15", "5.3", "15.2", "15.2"), ("NaCl", "1174.15", "5.1", "14.5", "14.7"),
        ("KCl", "1043.15", "5.4", "15.4", "15.7"), ("KCl", "1143.15", "5.2", "15.2", "15.1"),
        ("RbCl", "991.15", "5.4", "15.8", "15.9"), ("RbCl", "1091.15", "5.3", "14.9", "15.0"),
        ("CsCl", "918.15", "5.5", "16.0", "16.0"), ("CsCl", "1018.15", "5.3", "15.3", "15.6")],
    "JC non-polarisable FF MD": [
        ("LiCl", "878.15", "3.9", "11.4", "11.7"), ("LiCl", "978.15", "3.9", "11.4", "11.6"),
        ("NaCl", "1074.15", "4.9", "14.5", "14.3"), ("NaCl", "1174.15", "4.8", "14.1", "14.0"),
        ("KCl", "1043.15", "5.4", "16.0", "15.8"), ("KCl", "1143.15", "5.2", "15.3", "15.3"),
        ("RbCl", "991.15", "5.5", "16.4", "16.4"), ("RbCl", "1091.15", "5.4", "16.0", "15.8"),
        ("CsCl", "918.15", "5.8", "17.1", "17.2"), ("CsCl", "1018.15", "5.6", "16.6", "16.6")],
}
for meth, lst in P1_SIM.items():
    tab = "S1" if meth.startswith("WBK") else "S2"
    for salt, T, a, b, c in lst:
        m = salt[:-2]
        for cen, nb, v in ((m, "Cl", a), (m, m, b), ("Cl", "Cl", c)):
            add(cen, nb, "pure melt", T, salt, v, "0.5", "first min of g(r)", meth, "B",
                f"P1 Table {tab}", "S; source states CN uncertainty ~0.5")

# P1 Table S4 literature simulations: all chloride rows (8 read in session; the rest read by review 3)
for salt, T, v, meth, orig in [
        ("LiCl", "883", "4.9", "MC", "Lewis, Singer & Woodcock"), ("LiCl", "900", "4.05", "MC", "Wang et al. 2014 (P1 ref 27)"),
        ("LiCl", "1096", "4.05", "MD", "Lantelme & Turq 1982"), ("LiCl", "1050", "4.1", "MD", "Caccamo & Dixon 1980"),
        ("LiCl", "918", "4.2", "MD", "Caccamo & Dixon 1980"), ("LiCl", "923", "4.0", "MC", "Baranyai, Ruff & McGreevy 1986"),
        ("LiCl", "1073", "4.1", "MC", "Krogh-Moe"),
        ("NaCl", "1080", "4.76", "MD", "Wang et al. 2014 (P1 ref 27)"), ("NaCl", "1074", "5.2", "MC", "Lewis, Singer & Woodcock"),
        ("NaCl", "1123", "5.0", "MC", "Baranyai, Ruff & McGreevy 1986"),
        ("KCl", "1045", "5.3", "MC", "Woodcock & Singer"), ("KCl", "1306", "5.2", "MC", "Woodcock & Singer"),
        ("KCl", "1053", "5.18", "MD", "Wang et al. 2014"), ("KCl", "1096", "5.13", "MD", "Lantelme & Turq"),
        ("KCl", "1043", "5.2", "MD", "Caccamo & Dixon 1980"), ("KCl", "1073", "5.8", "MC", "Baranyai, Ruff & McGreevy 1986"),
        ("KCl", "1173", "4.7", "MD", "Takagi, Ohno & Furukawa 1979 (P1 ref list describes it as experimental)"),
        ("KCl", "1173", "5.1", "MD", "Takagi, Ohno & Furukawa 1979 (P1 ref list describes it as experimental)"),
        ("RbCl", "1000", "5.41", "MD", "Wang et al. 2014"), ("RbCl", "995", "5.8", "MC", "Lewis, Singer & Woodcock"),
        ("RbCl", "1023", "5.7", "MC", "Baranyai, Ruff & McGreevy 1986"),
        ("CsCl", "900", "5.7", "MD", "Wang et al. 2014"), ("CsCl", "973", "5.8", "MC", "Baranyai, Ruff & McGreevy 1986")]:
    note = f"original: {orig}"
    if "Baranyai" in orig:
        note += "; whether this MC was constrained by diffraction data (reverse MC) was not checked"
    add(salt[:-2], "Cl", "pure melt", T, salt, v, "", "", meth, "C", "P1 Table S4", note)

# ---- P2 Roy et al. 2020, n_ave to r_dagger, 1100 K ----
for salt, rd, pim, rim in [("LiCl", "3.6", "4.34", "3.90"), ("NaCl", "4.1", "5.16", "4.82"),
                           ("KCl", "4.5", "5.18", "4.89")]:
    add(salt[:2] if salt != "KCl" else "K", "Cl", "pure melt", "1100", salt, pim, "", rd,
        "PIM MD", "B", "P2 Table 2", "S; r_dagger = first-shell boundary from cation-anion g(r)")
    add(salt[:2] if salt != "KCl" else "K", "Cl", "pure melt", "1100", salt, rim, "", rd,
        "RIM MD", "B", "P2 Table 2", "S; same r_dagger as PIM")

# ---- P10 Lu et al. 2021 ----
add("Na", "Cl", "pure melt", "1100", "NaCl", "5.06", "", "first min of g(r) (value not given)", "PIM MD", "B", "P10", "S")
add("Na", "Cl", "pure melt", "", "NaCl", "4.83", "", "", "experiment", "C", "P10",
    "quoted from McGreevy, Solid State Phys. 40, 247 (1987); T not given by P10")
add("Mg", "Cl", "pure melt", "", "MgCl2", "5.21 to 4.72", "", "first min of g(r)", "PIM MD", "B", "P10",
    "S; decreases with T; end temperatures only in P10 Fig. 2")

# ---- 2. Alkaline-earth chlorides ----
for m, salt, T, v, orig in [("Mg", "MgCl2", "998", "4.3", "Biggin, Gay & Enderby 1984"),
                            ("Ca", "CaCl2", "1093", "5.3", "Biggin & Enderby 1981"),
                            ("Sr", "SrCl2", "1198", "6.9", "McGreevy & Mitchell 1982"),
                            ("Ba", "BaCl2", "1298", "6.4", "Edwards, Howe, Enderby & Page 1978")]:
    add(m, "Cl", "pure melt", T, salt, v, "", "", "neutron diffraction (isotopic substitution)", "C", "P23",
        f"value as compiled by Enderby 1987 [Q21] and quoted in P23; T from P23 for original {orig}")
add("Sr", "Cl", "pure melt", "1198", "SrCl2", "6.9", "", "", "thermal neutron scattering", "C", "Q13",
    "Q13 quotes McGreevy & Mitchell: CN 6.9 at 925 C (=1198 K); second independent quotation")
for m, salt, T, v, note in [
        ("Mg", "MgCl2", "998", "4.3", "carried from rev 2.1 (P22 Table 2 read in an earlier session); not re-read"),
        ("Ca", "CaCl2", "1093", "5.4", "carried from rev 2.1 (P22 Table 2); not re-read"),
        ("Sr", "SrCl2", "1198", "6.9", "carried from rev 2.1 (P22 Table 2; P22 text says 6.6); not re-read"),
        ("Ba", "BaCl2", "1298", "7.7", "carried from rev 2.1 (P22 Table 2); not re-read; CONFLICTS with 6.4 (P23)")]:
    add(m, "Cl", "pure melt", T, salt, v, "", "", "neutron diffraction (compilation)", "C", "P22", note)
add("Ba", "Cl", "pure melt", "", "BaCl2", "8 (\"eightfold\")", "", "", "unknown (qualitative)", "C", "Q27",
    "citing-paper snippet (Edgar et al. 2010) on ResearchGate; possibly paraphrased; source of the 8 unknown")
add("Mg", "Cl", "pure melt", "998", "MgCl2", "5.2", "0.9 (SD of per-ion distribution)", "", "MD fitted to diffraction (BELION)", "B", "P23", "S")
add("Mg", "Cl", "pure melt", "1200", "MgCl2", "4.4", "", "not stated for pure MgCl2", "AIMD", "B", "P6", "S")
add("Ca", "Cl", "pure melt", "1200", "CaCl2", "6.0", "", "not stated", "classical MD (BMH)", "B", "P8",
    "S; P8 also quotes 5.4 (exp, 1093 K) but its citation of Biggin is defective")
# Christian et al. USPEX/DFT (pure Sr-Cl binary)
add("Sr", "Cl", "pure Sr-Cl (crystal-structure ensemble)", "298", "Sr-Cl binary", "6.61", "", "", "USPEX/DFT Boltzmann-weighted crystal ensemble", "B", "Q13",
    "not an MD liquid; model average. Q13 also gives 5.5 (500 C) and 4.88 (1000 C) in a sentence whose Sr/Zr attribution is ambiguous; omitted")

# D rows: anion-centred CN in MCl2: N(Cl-M) = N(M-Cl)/2
for m, salt, T, v, src in [("Mg", "MgCl2", "998", 4.3, "P23"), ("Ca", "CaCl2", "1093", 5.3, "P23"),
                           ("Ca", "CaCl2", "1093", 5.4, "P22"),
                           ("Sr", "SrCl2", "1198", 6.9, "P23/Q13"), ("Ba", "BaCl2", "1298", 6.4, "P23"),
                           ("Ba", "BaCl2", "1298", 7.7, "P22")]:
    add("Cl", m, "pure melt", T, salt, f"{v/2:.2f}", "", "", "derived", "D", src,
        f"N(Cl-M) = N(M-Cl) x_M/x_Cl = N(M-Cl)/2 with N(M-Cl) = {v}")
for m in ["Li", "Na", "K", "Rb", "Cs"]:
    add("Cl", m, "pure melt", "", f"{m}Cl", "= N(M-Cl)", "", "", "derived", "D", "-",
        "1:1 salt: N(Cl-M) = N(M-Cl) for the same cutoff")

# ---- 3. Mixed / dilute melts ----
add("Rb", "Cl", "pure melt", "1050", "RbCl", "5.3", "not given", "none (first-shell fit)", "XAFS (Rb K-edge)", "B", "P3", "fitted <r> = 3.20 A")
add("Rb", "Cl", "dilute in LiCl-KCl", "873", "15% RbCl in LiCl-KCl eutectic", "5.8", "not given",
    "none (first-shell fit)", "XAFS (Rb K-edge)", "B", "P3", "fitted <r> = 3.25 A; % basis not stated (probably mol%); 177 K below the pure-RbCl point")
for T, v, u, R in [("298", "7.5", "1.3", "3.01"), ("673", "6.5", "1.4", "2.94"), ("873", "6.5", "1.9", "2.92"),
                   ("1073", "5.9", "0.89", "2.93"), ("1273", "5.8", "1.1", "3.03")]:
    state = "solid" if int(T) < 1055 else "melt"
    add("Sr", "Cl", f"dilute in NaCl ({state})", T, "5 mol% SrCl2 in NaCl", v, u,
        "none (first-shell fit)", "EXAFS (Sr K-edge)", "B", "P5 SI Table S1",
        f"fitted R = {R} A. SI Table S1 (p. S3) read once in session from the Clemson-hosted PDF; reviewer re-fetch did not return the SI. 298 and 1273 K N also in main text; 873 K N also in Q13 Table 3. Mix m.p. 782+/-2 C (=1055 K, P5 Table 1, origin not stated)")
for T, v in [("298", "6.2"), ("873", "5.8"), ("1273", "5.6")]:
    add("Sr", "Cl", "dilute in NaCl (model)", T, "<=10% Sr in Sr-Na-Cl", v, "", "", "USPEX/DFT Boltzmann crystal ensemble", "B", "Q13",
        "model average compared with P5")
add("Sr", "Cl", "dilute in LiCl-KCl, with oxide", "", "SrO in LiCl-KCl (and LiCl-LiF)", "6", "not given", "",
    "EXAFS", "A", "Q15", "study also used Raman and DFT; Sr environment = 1 O + 6 Cl; abstract gives 298-923 K for the study, EXAFS T not stated")
add("Sr", "Cl", "SrCl2-KCl-MgCl2 melt", "873-1173", "SrCl2-KCl-MgCl2", "6-7 (dominant)", "", "", "DP-MD", "A", "P19",
    "abstract: 6- and 7-fold dominant; no mean CN")
for ln, li, k in [("Sm", "4.3", "7.1"), ("Ce", "4.1", "6.8"), ("Eu", "4.1", "7.2")]:
    comp = f"LiCl-KCl eutectic 59.2 mol% LiCl + 1 {ln}3+ (36 Li/25 K/64 Cl box)"
    add("Li", "Cl", "LiCl-KCl eutectic", "773", comp, li, "", "RDF integrated to an unnumbered R_min (arrows in figures)", "AIMD (VASP PBE-D3)", "B", "Q14",
        "S; source says 'Li+ CN' without naming the partner; Cl assumed")
    add("K", "Cl", "LiCl-KCl eutectic", "773", comp, k, "", "RDF integrated to an unnumbered R_min (arrows in figures)", "AIMD (VASP PBE-D3)", "B", "Q14",
        "S; source says 'K+ CN' without naming the partner; Cl assumed; high vs pure KCl ~5.2-5.4 (which is at 1043-1143 K, ~300 K hotter); definition-dependent")
# D: eutectic anion-centred from Q14 means (Li 4.17, K 7.03), using Q14's box counts 36 Li / 25 K / 64 Cl
li_m = (4.3 + 4.1 + 4.1) / 3
k_m = (7.1 + 6.8 + 7.2) / 3
nli, nk, ncl = 36, 25, 64
comp = "Q14 box 36 Li/25 K/64 Cl (+1 Ln)"
cav = "inherits Q14 caveats: R_min not numbered; neighbour of 'K+ CN' assumed Cl"
add("Cl", "Li", "LiCl-KCl eutectic", "773", comp, f"{li_m*nli/ncl:.2f}", "", "", "derived", "D", "Q14",
    f"N(Cl-Li) = N(Li-Cl)*N_Li/N_Cl, N(Li-Cl)={li_m:.2f} (mean of 3 runs); {cav}")
add("Cl", "K", "LiCl-KCl eutectic", "773", comp, f"{k_m*nk/ncl:.2f}", "", "", "derived", "D", "Q14",
    f"N(Cl-K) = N(K-Cl)*N_K/N_Cl, N(K-Cl)={k_m:.2f}; {cav}")
add("Cl", "Li+K", "LiCl-KCl eutectic", "773", comp, f"{li_m*nli/ncl + k_m*nk/ncl:.2f}", "", "", "derived", "D", "Q14",
    "sum of the two rows above (Ln neighbours excluded); ideal Ln-free 59.2:40.8 melt would give 2.47 + 2.87 = 5.34")

# ---- 4. Oxide ----
add("Sr", "O", "dilute in LiCl-KCl (and LiCl-LiF)", "", "SrO in LiCl-KCl", "1", "not given", "",
    "EXAFS", "A", "Q15", "study also used Raman and DFT; Sr-centred; N(O-Sr) not derivable; SrOCl observed (298-923 K); 'more stable' from DFT/thermodynamic estimates")
for c, n, tag, note in [
        ("Li", "O", "P13/P15", "direct sources not readable (ScienceDirect)"),
        ("O", "Li", "P13/P15", "direct sources not readable (ScienceDirect)"),
        ("Ca", "O", "Q16/Q17", "speciation only: Ca4OCl6 below 1104 K; Ca2+ layer around O2-; no number"),
        ("O", "Ca", "Q16/Q17", "speciation only; no number"),
        ("Ba", "O", "none found", ""), ("O", "Ba", "none found", ""),
        ("Na", "O", "none found", ""), ("O", "Na", "none found", ""),
        ("K", "O", "none found", ""), ("O", "K", "none found", ""),
        ("Mg", "O", "none found", ""), ("O", "Sr", "Q15", "Sr-centred value only (Sr-O = 1)")]:
    add(c, n, "chloride melt with dissolved oxide", "", "", "unavailable", "", "", "", "-", tag, note)

# ---- Other values quoted in the Markdown ----
add("Li", "Cl", "LiCl-KCl eutectic", "", "LiCl-KCl eutectic", "~4", "", "", "X-ray diffraction + MD", "C", "-",
    "quoted by a citing paper; cited DOI 10.1023/a:1023310006090 (to confirm)")
add("Y", "Cl", "15 mol% YCl3 in LiCl-KCl (trivalent reference)", "973", "15 mol% YCl3 in LiCl-KCl", "6.0", "0.3", "none (first-shell fit)",
    "EXAFS (Y K-edge)", "B", "Q19", "fitted R = 2.68 A; 30 and 50 mol%: 6.1 +/- 0.4, 5.9 +/- 0.4; not a target species")
add("M3+", "Cl", "dilute M3+ in alkali chlorides (trivalent reference)", "", "", "6 (\"overwhelmingly likely\")", "", "", "review of simulations", "B", "Q28",
    "qualitative; not a target species")
for c, n, env in [("Cs", "Cl", "dilute in LiCl-KCl or NaCl"), ("Ba", "Cl", "dilute in LiCl-KCl or NaCl"),
                  ("Na", "Cl", "dilute in LiCl-KCl"), ("Sr", "Cl", "dilute in LiCl-KCl without oxide"),
                  ("M", "M", "pure MgCl2/CaCl2/SrCl2/BaCl2"), ("Cl", "Cl", "pure MgCl2/CaCl2/SrCl2/BaCl2"),
                  ("Li/K", "Li/K", "LiCl-KCl eutectic"), ("Cl", "Cl", "LiCl-KCl eutectic")]:
    add(c, n, env, "", "", "unavailable", "", "", "", "-", "none found", "no numeric value in any readable source")

# ---- E: working recommendations for the generator (lta_coordination_numbers.md section 1) ----
add("Ba", "Cl", "pure melt (working value)", "1298", "BaCl2", "7", "1", "", "estimate", "E", "P22/P23/Q27",
    "E: brackets the unresolved 6.4 (P23) vs 7.7 (P22) and the qualitative 8 (Q27); not an average")
add("Cl", "Ba", "pure melt (working value)", "1298", "BaCl2", "3.5", "0.5", "", "estimate", "E", "-", "N(Ba-Cl)/2")
add("Sr", "Cl", "dilute in LiCl-KCl, 773-823 K (working value)", "773-823", "dilute Sr2+ in LiCl-KCl", "6", "1", "", "estimate", "E", "P5/Q15/P19/Q13",
    "E: brackets P5 (5.8-5.9 in molten NaCl), Q15 (6 Cl with 1 O), P19 (6-7); Q13 trend favours slightly higher CN at lower T")

# ---- 5. Inside LTA (single crystal, dehydrated or solution-loaded) ----
add("Sr", "O(framework)", "dehydrated Sr/Tl-A, 6-ring", "294", "Sr1.6Tl8.8-A and Sr5.45Tl1.1-A", "3", "", "n/a (crystal structure)",
    "single-crystal XRD", "A", "Q25", "all Sr at 6-rings, 3-coordinate to O3; Sr-O3 2.21(2) and 2.31(1) A")
add("Ba", "O(framework)", "dehydrated Ba-A", "", "Ba-A series", "3", "", "", "single-crystal XRD", "A", "Q26",
    "abstract: 3-coordinate Ba2+ forms on dehydration; may disrupt framework")
add("Ag", "O(framework)", "K9(K4Br)0.75Ag4Br4-A, 6-ring", "294", "", "3", "", "n/a (crystal structure)", "single-crystal XRD", "A", "R28", "Ag-O 2.52(3) A; each Ag also bonded within Ag4Br4")
add("Br", "K", "K4Br3+ unit in alpha-cage", "294", "K9(K4Br)0.75Ag4Br4-A", "4", "", "", "single-crystal XRD", "A", "R28",
    "by unit formula; Br opposite a 4-ring, held by 2 8-ring + 2 6-ring K")
add("Br", "K", "Ag4Br4 cluster in beta-cage", "294", "K9(K4Br)0.75Ag4Br4-A", "1", "", "n/a (crystal structure)", "single-crystal XRD", "A", "R28",
    "Br-K 3.00(4) A; each cluster Br coordinates one K+ through a 6-ring")
add("Cl", "K", "K3Cl2+ unit in alpha-cage", "294", "R27 composition", "3", "", "", "single-crystal XRD", "A", "R27", "by unit formula")
add("Cl", "K", "Ag4Cl4 cluster in beta-cage", "294", "R27 composition", "1", "", "n/a (crystal structure)", "single-crystal XRD", "A", "R27", "Cl-K 2.70(3) A, through a 6-ring")
add("I", "K", "K4I3+ unit in alpha-cage", "", "K9(K4I)(Ag4I4)0.5-A", "4", "", "", "single-crystal XRD", "A", "R34", "by unit formula")
add("Tl", "Cl", "Tl7Cl6+ in alpha-cage", "", "(WCl)0.6(Tl7Cl6)0.5Tl8.5-A", "6", "", "", "single-crystal XRD", "A", "R30", "central Tl at the alpha-cage centre, octahedrally coordinated by 6 Cl")
add("Cl", "Tl", "Tl7Cl6+ in alpha-cage", "", "(WCl)0.6(Tl7Cl6)0.5Tl8.5-A", "2", "", "", "single-crystal XRD", "A", "R30", "1 central Tl + 1 8-ring Tl, by the described topology")

with open(sys.argv[1], "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "rows")
