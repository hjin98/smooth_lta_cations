"""Build the pair-distance table for the random molten-salt structure generator.

Revision 2 (2026-10-17), after independent review. Every number traces to
molten_salt_pair_distances.md (reference tags [Pn]).

Grades:
  M  = measured in a molten salt (neutron/X-ray diffraction or XAFS)
  S  = simulation of the molten salt (MD / PIM / AIMD), read from the source
  E  = estimate from the model in this script (formula stated in 'basis')

Distances in angstrom. r_cut = g(r) onset (closest approach),
r_peak = first maximum, r_min1 = first minimum after the peak.
Reference temperature: near the melting point (~950-1300 K).
"""
import csv
import math

# ---------------------------------------------------------------------------
# 1. Inputs taken from sources
# ---------------------------------------------------------------------------
# Selection rule for measured peaks: arithmetic mean of the NEUTRON values in
# the P1 compilation (Table S3); X-ray / XAFS values are listed in the md only
# as cross-checks. Divalent salts: compiled XRD/ND values [P22].
R_MCL = {
    # ion: (r_peak, grade, basis)
    "Li": (2.38, "M", "mean of neutron 2.30 (958 K) and 2.45 [P1]; X-ray 2.40-2.55"),
    "Na": (2.75, "M", "mean of neutron 2.72 and 2.78 (1148 K) [P1]; X-ray 2.77-2.80"),
    "K":  (3.08, "M", "mean of neutron 3.06 (1073 K) and 3.10 [P1]; X-ray 3.05-3.14"),
    "Rb": (3.18, "M", "neutron 3.18 (1023 K) [P1]; Rb K-edge XAFS 3.20 (1050 K) [P3]"),
    "Cs": (3.39, "M", "mean of neutron 3.38, 3.40, 3.38 (968-973 K) [P1]"),
    "Sr": (2.90, "M", "molten SrCl2 1198 K, N=6.9 (McGreevy & Mitchell 1982 via [P22])"),
    "Ba": (3.10, "M", "molten BaCl2 1298 K, N=7.7 (Edwards et al. 1978 via [P22])"),
}

# First minimum of M-Cl where reported: PIM first-shell boundary r†, 1100 K [P2]
R_MIN1_MCL_REPORTED = {"Li": 3.6, "Na": 4.1, "K": 4.5}
# Polarizable-model (WBK) first peaks, lower temperature of each pair [P1, Table S1]
R_PEAK_POLARIZABLE = {"Li": 2.38, "Na": 2.64, "K": 2.96, "Rb": 3.10, "Cs": 3.28}

# Divalent min/peak ratio from AIMD coordination cutoffs [P6]:
# U(III)-Cl 4.0/2.74 (stated to be the RDF minimum); Mg-Cl 3.4/2.40 (cutoff
# used for CN counting, not explicitly stated to be the minimum)
RATIO_MIN_PEAK_DIVALENT = (3.4 / 2.40 + 4.0 / 2.74) / 2.0

# Shannon effective ionic radii, CN VI (angstrom) [P7]
SHANNON_VI = {"Li": 0.76, "Na": 1.02, "K": 1.38, "Rb": 1.52, "Cs": 1.67,
              "Mg": 0.72, "Ca": 1.00, "Sr": 1.18, "Ba": 1.35,
              "Cl": 1.81, "O": 1.40}

# Measured molten MCl2 data used to calibrate divalent geometry [P22, P23]
DIVALENT_MELTS = {
    #       M-Cl,  Cl-Cl, M-M
    "Mg": (2.42, 3.56, 3.81),   # 998 K, Biggin, Gay & Enderby 1984
    "Ca": (2.78, 3.73, 3.60),   # 1093 K, Biggin & Enderby 1981
}

# Like-pair first peaks in pure molten alkali chlorides: mean of neutron values
# [P1 Table S3]; CsCl rows giving 3.85/3.85 (973 K) rejected as anomalous.
LIKE_MEASURED = {
    # salt: (r_MM, r_ClCl)
    "Li": (3.75, 3.75),   # neutron 3.70 (958 K) and 3.80
    "Na": (3.93, 3.91),   # neutron 3.89/3.90 and 3.96/3.91 (1148 K)
    "K":  (4.77, 4.76),   # neutron 4.84/4.82 (1073 K) and 4.70/4.70; X-ray 4.30-4.40 disagrees
    "Rb": (4.86, 4.80),   # neutron (1023 K)
    "Cs": (4.95, 4.85),   # neutron (968 K)
}

# Crystal nearest-neighbour M-O distances
ANTIFLUORITE_A = {"Li": 4.61, "Na": 5.55, "K": 6.449, "Rb": 6.74}   # M2O, a [P9]
ROCKSALT_NN = {"Sr": 2.57, "Ba": 2.76}                              # MO, a/2 [P24]

CATIONS = ["Li", "Na", "K", "Rb", "Cs", "Sr", "Ba"]
CHARGE = {"Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Sr": 2, "Ba": 2,
          "Cl": -1, "O": -2}

# ---------------------------------------------------------------------------
# 2. Model constants
# ---------------------------------------------------------------------------
ALK_OFF = {m: R_MCL[m][0] - (SHANNON_VI[m] + SHANNON_VI["Cl"]) for m in ["Na", "K", "Rb", "Cs"]}
ALK_OFF_MEAN = sum(ALK_OFF.values()) / len(ALK_OFF)
DIV_OFF = {m: v[0] - (SHANNON_VI[m] + SHANNON_VI["Cl"]) for m, v in DIVALENT_MELTS.items()}
DIV_OFF_MEAN = sum(DIV_OFF.values()) / len(DIV_OFF)

# M-Cl r_min1/r_peak in the polarizable-model frame (r† and WBK peaks)
RATIO_MIN_PEAK_MONO = sum(R_MIN1_MCL_REPORTED[m] / R_PEAK_POLARIZABLE[m]
                          for m in R_MIN1_MCL_REPORTED) / 3.0

# Onset ratios r_cut / r_peak (no measured onset readable; md §4.3).
# Like-pair ratio and anion contact floor lowered in rev 2.1 so that a random
# insertion generator can reach experimental density (packing test, md §4.3).
ONSET_UNLIKE = 0.82
ONSET_LIKE = 0.72
FLOOR_FACTOR = 0.82          # anion-anion r_cut >= FLOOR_FACTOR * contact sum
RATIO_MIN_PEAK_LIKE = 1.45   # assumption, md §4.2

# Bridge factor F: r_XY = F*sqrt(r_XA^2 + r_YA^2) for same-sign X, Y sharing A.
# Alkali: fitted to measured like pairs of NaCl, KCl, RbCl, CsCl (Li excluded:
# Cl-Cl in LiCl is anion-contact limited).
def fit_F():
    vals = []
    for m in ["Na", "K", "Rb", "Cs"]:
        r = R_MCL[m][0]
        mm, xx = LIKE_MEASURED[m]
        vals += [mm / (math.sqrt(2) * r), xx / (math.sqrt(2) * r)]
    return sum(vals) / len(vals), min(vals), max(vals)

F_ALK, F_LO, F_HI = fit_F()
# Divalent: calibrated on molten CaCl2 (closest to Sr/Ba in coordination);
# MgCl2 (network melt, CN ~4) shown for comparison only.
F_DIV_CATION = DIVALENT_MELTS["Ca"][2] / (math.sqrt(2) * DIVALENT_MELTS["Ca"][0])
F_DIV_ANION = DIVALENT_MELTS["Ca"][1] / (math.sqrt(2) * DIVALENT_MELTS["Ca"][0])
F_MG_CATION = DIVALENT_MELTS["Mg"][2] / (math.sqrt(2) * DIVALENT_MELTS["Mg"][0])
F_MG_ANION = DIVALENT_MELTS["Mg"][1] / (math.sqrt(2) * DIVALENT_MELTS["Mg"][0])

# Oxide geometry: O-M-O bridge angle of the parent crystal
TETRA = 2 * math.sin(math.radians(109.4712) / 2)   # 1.633: antifluorite, O-O = a/sqrt(2)
OCTA = math.sqrt(2)                                 # rock salt, O-O = a/sqrt(2)
# O-bridged cation-cation: M-O-M angle of the parent crystal
ANTIFL_MM = 2 / math.sqrt(3)   # 1.155: antifluorite M-M = a/2, M-O-M = 70.5 deg
ROCK_MM = math.sqrt(2)         # rock salt M-M = a/sqrt(2), M-O-M = 90 deg

FLOOR = {("Cl", "Cl"): 2 * SHANNON_VI["Cl"], ("O", "O"): 2 * SHANNON_VI["O"],
         ("Cl", "O"): SHANNON_VI["Cl"] + SHANNON_VI["O"]}

# ---------------------------------------------------------------------------
# 3. Cation-anion distances
# ---------------------------------------------------------------------------
def mcl_min1(m, rpk):
    if m in R_MIN1_MCL_REPORTED:
        return R_MIN1_MCL_REPORTED[m], "S", "PIM first-shell boundary r†, 1100 K [P2]"
    if CHARGE[m] == 1:
        r = RATIO_MIN_PEAK_MONO * R_PEAK_POLARIZABLE[m]
        return round(r, 2), "E", (f"{RATIO_MIN_PEAK_MONO:.3f} x polarizable-model peak "
                                  f"{R_PEAK_POLARIZABLE[m]} [P1,P2]")
    r = RATIO_MIN_PEAK_DIVALENT * rpk
    return round(r, 2), "E", f"{RATIO_MIN_PEAK_DIVALENT:.3f} x r_peak (AIMD cutoff ratio [P6])"

def mo_peak(m):
    if m in ANTIFLUORITE_A:
        d = ANTIFLUORITE_A[m] * math.sqrt(3) / 4
        return round(d - 0.02, 2), "E", (f"M2O crystal M-O {d:.3f} (a={ANTIFLUORITE_A[m]} [P9]) "
                                         f"- 0.02 melt contraction")
    if m in ROCKSALT_NN:
        d = ROCKSALT_NN[m]
        return round(d - 0.02, 2), "E", f"MO rock-salt M-O {d:.2f} [P24] - 0.02 melt contraction"
    s = SHANNON_VI[m] + SHANNON_VI["O"]
    return round(s + ALK_OFF_MEAN, 2), "E", (f"Shannon VI sum {s:.2f} + mean alkali M-Cl "
                                             f"offset {ALK_OFF_MEAN:+.3f}")

def f_cation_pair(x, y):
    fx = F_DIV_CATION if CHARGE[x] == 2 else F_ALK
    fy = F_DIV_CATION if CHARGE[y] == 2 else F_ALK
    return (fx + fy) / 2

# ---------------------------------------------------------------------------
# 4. Assemble rows
# ---------------------------------------------------------------------------
def row(i, j, rcut, rpk, rmin, grade, basis, env=""):
    return {"ion_i": i, "ion_j": j, "environment": env,
            "r_cut_A": f"{rcut:.2f}", "r_peak_A": f"{rpk:.2f}",
            "r_min1_A": f"{rmin:.2f}", "grade_peak": grade, "basis": basis}

def build():
    rows, mcl, mo = [], {}, {}
    for m in CATIONS:
        rpk, g, b = R_MCL[m]
        rmin, gm, bm = mcl_min1(m, rpk)
        mcl[m] = rpk
        rows.append(row(m, "Cl", ONSET_UNLIKE * rpk, rpk, rmin, g,
                        f"peak: {b}; min ({gm}): {bm}; cut (E): {ONSET_UNLIKE} x peak"))
    for m in CATIONS:
        rpk, g, b = mo_peak(m)
        mo[m] = rpk
        rows.append(row(m, "O", ONSET_UNLIKE * rpk, rpk, RATIO_MIN_PEAK_DIVALENT * rpk, g,
                        f"peak: {b}; min (E): {RATIO_MIN_PEAK_DIVALENT:.3f} x peak; "
                        f"cut (E): {ONSET_UNLIKE} x peak"))
    for a in range(len(CATIONS)):
        for b in range(a, len(CATIONS)):
            x, y = CATIONS[a], CATIONS[b]
            if x == y and x in LIKE_MEASURED:
                rpk, g = LIKE_MEASURED[x][0], "M"
                basis = "mean of neutron values [P1]"
                if x == "K":
                    basis += "; X-ray gives 4.30-4.40, treat as +/-0.25"
            else:
                f = f_cation_pair(x, y)
                rpk, g = f * math.sqrt(mcl[x] ** 2 + mcl[y] ** 2), "E"
                basis = f"{f:.3f} x sqrt(r_{x}Cl^2 + r_{y}Cl^2)"
            if CHARGE[x] == 2 and CHARGE[y] == 2:
                basis += "; LOW CONFIDENCE: crystal-like geometry would give up to ~1.2 A longer"
            rows.append(row(x, y, ONSET_LIKE * rpk, rpk, RATIO_MIN_PEAK_LIKE * rpk, g,
                            f"peak: {basis}; min (E): {RATIO_MIN_PEAK_LIKE} x peak; "
                            f"cut (E): {ONSET_LIKE} x peak", env="Cl-bridged"))
    # cation-cation sharing an O2- (oxide-bridged)
    for a in range(len(CATIONS)):
        for b in range(a, len(CATIONS)):
            x, y = CATIONS[a], CATIONS[b]
            fx = (ROCK_MM if CHARGE[x] == 2 else ANTIFL_MM) / math.sqrt(2)
            fy = (ROCK_MM if CHARGE[y] == 2 else ANTIFL_MM) / math.sqrt(2)
            f = (fx + fy) / 2
            rpk = f * math.sqrt(mo[x] ** 2 + mo[y] ** 2)
            basis = (f"{f:.3f} x sqrt(r_{x}O^2 + r_{y}O^2) (parent-crystal M-O-M angle: "
                     f"antifluorite 70.5 deg / rock-salt 90 deg)")
            rows.append(row(x, y, ONSET_LIKE * rpk, rpk, RATIO_MIN_PEAK_LIKE * rpk, "E",
                            f"peak: {basis}; min (E): {RATIO_MIN_PEAK_LIKE} x peak; "
                            f"cut (E): {ONSET_LIKE} x peak", env="O-bridged"))
    for m in CATIONS:
        div = CHARGE[m] == 2
        for (p, q) in [("Cl", "Cl"), ("O", "O"), ("Cl", "O")]:
            if p == q == "Cl" and m in LIKE_MEASURED:
                rpk, g, basis = LIKE_MEASURED[m][1], "M", "mean of neutron values [P1]"
            elif p == q == "O":
                ang = OCTA if div else TETRA
                rpk, g = ang * mo[m], "E"
                basis = (f"{ang:.3f} x r_{m}O ({'rock-salt 90 deg' if div else 'antifluorite 109.5 deg'} bridge)")
            else:
                f = F_DIV_ANION if div else F_ALK
                rp = mcl[m] if p == "Cl" else mo[m]
                rq = mcl[m] if q == "Cl" else mo[m]
                rpk, g = f * math.sqrt(rp ** 2 + rq ** 2), "E"
                basis = f"{f:.3f} x sqrt(r_{m}{p}^2 + r_{m}{q}^2)"
            floor = FLOOR[(p, q)]
            if rpk < floor + 0.05:
                rpk = floor + 0.05
                basis += f"; raised to contact floor {floor:.2f}+0.05"
            rcut = max(ONSET_LIKE * rpk, FLOOR_FACTOR * floor)
            rows.append(row(p, q, rcut, rpk, RATIO_MIN_PEAK_LIKE * rpk, g,
                            f"peak: {basis}; min (E): {RATIO_MIN_PEAK_LIKE} x peak; "
                            f"cut (E): max({ONSET_LIKE} x peak, {FLOOR_FACTOR} x contact)",
                            env=f"around {m}"))
    return rows

def main():
    print("Alkali M-Cl offsets:", {k: round(v, 3) for k, v in ALK_OFF.items()},
          "mean", round(ALK_OFF_MEAN, 4))
    print("Divalent offsets:", {k: round(v, 3) for k, v in DIV_OFF.items()},
          "mean", round(DIV_OFF_MEAN, 3))
    for m in ["Sr", "Ba"]:
        est = SHANNON_VI[m] + SHANNON_VI["Cl"] + DIV_OFF_MEAN
        print(f"  cross-check {m}-Cl model {est:.2f} vs measured {R_MCL[m][0]}")
    print("min/peak mono:", round(RATIO_MIN_PEAK_MONO, 4), "divalent:", round(RATIO_MIN_PEAK_DIVALENT, 4))
    print("F_alk:", round(F_ALK, 4), "range", round(F_LO, 4), round(F_HI, 4))
    print("F_div (CaCl2) cation/anion:", round(F_DIV_CATION, 4), round(F_DIV_ANION, 4),
          "| MgCl2 cation/anion:", round(F_MG_CATION, 4), round(F_MG_ANION, 4))
    rows = build()
    with open("molten_salt_pair_distances.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        assert float(r["r_cut_A"]) < float(r["r_peak_A"]) < float(r["r_min1_A"]), r
    print(len(rows), "rows written")
    for r in rows:
        print(f"{r['ion_i']:>3}-{r['ion_j']:<3} {r['environment']:<12} "
              f"{r['r_cut_A']:>5} {r['r_peak_A']:>5} {r['r_min1_A']:>5} {r['grade_peak']}")

if __name__ == "__main__":
    main()
