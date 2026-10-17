---
title: "Current Group-I/II Dry-Volume Reference Model (v7.15)"
date: "2026-10-05"
toc: true
toc-depth: 3
geometry: margin=0.78in
fontsize: 9.5pt
header-includes:
  - |
    \usepackage{booktabs}
    \usepackage{longtable}
    \usepackage{microtype}
  - |
    \setlength{\emergencystretch}{3em}
---

# 1. Status and authority

This document is the **accepted-current scientific/numerical authority** for the dry framework-compensating composition-volume model retained by v7.15. It owns:

- the supported framework normalization for the empirical volume model;
- charge-equivalent interpolation;
- the `/24` denominator;
- species-specific empirical/effective endpoint values and classifications;
- placeholder behavior;
- source-to-model transformations and their limitations;
- `volume_factor` semantics.

The chemistry, supported species, exact framework-compensating population typing, and formal-charge table are owned by `v7_13_group_i_ii_cation_extension.md`. The sparse chloride extension and the explicit exclusion of occlusion-counterion/Cl populations from this automatic volume model are owned by `v7_15_chloride_occlusion_extension.md`. Module specifications and code concretize these authorities.

# 2. Background and terminology

For the supported Al24Si24O96 LTA cell, let:

- $N_i=N_i^{\rm fw}$ be the exact integer number of **framework-compensating** cations of species $i$;
- $z_i$ be its formal positive charge (+1 for Group I, +2 for Group II);
- $V_i^{(24)}$ be the dry-LTA **reference corner volume** assigned to species $i$, normalized to one 24-Al reference cell;
- $w_i=z_iN_i/24$ be the species' charge-equivalent fraction;
- $f$ be `volume_factor`.

An endpoint may be **DIRECT**, **EXTRAPOLATED**, **PLACEHOLDER**, or **QUESTIONABLE**. These classifications describe the evidence/source-to-model status of the endpoint; they are not uncertainty distributions.

# 3. Supported framework scope

The model is defined only for the current production normalization:

```text
Al24Si24O96
Si:Al = 1:1
sum_i z_i N_i^fw = 24
```

The empirical $V_i^{(24)}$ values are 24-Al reference-cell volumes. Non-24-Al LTA cells, other Si:Al ratios, and alternative cell normalizations are outside v7.15 scope.

A future extension must explicitly define both the endpoint normalization and the mapping from normalized composition volume to generated-cell volume. Support must not be inferred merely by replacing the fixed denominator 24 with an arbitrary $N_{\rm Al}$.

# 4. Governing equations

For every supported neutral composition,

$$
\sum_i z_iN_i=24,
\qquad
w_i=\frac{z_iN_i}{24},
\qquad
\sum_i w_i=1,
$$

and

$$
V_{\rm ref}
=\sum_i w_iV_i^{(24)}
=\frac{1}{24}\sum_i z_iN_iV_i^{(24)}.
$$

The final target volume is

$$
V_{\rm target}=fV_{\rm ref}.
$$

Isotropic lattice scaling is

$$
s=\left(\frac{V_{\rm target}}{V_{\rm old}}\right)^{1/3},
\qquad
A_{\rm new}=sA_{\rm old}.
$$

Fractional coordinates are unchanged by this pure lattice scaling.

For Li/Na/K-only cells, all $z_i=1$, so this expression is exactly the historical count-over-24 interpolation.

# 5. Source-to-model transformations

## 5.1 Cell normalization

For a conventional cubic Fm-3c LTA cell with lattice parameter $a$, the package's 24-Al reference volume is

$$
V^{(24)}_{\rm Fm\bar 3c}(a)=\frac{a^3}{4}.
$$

For the 12-Al Pm-3m pseudocell convention used for the Ba source,

$$
V^{(24)}_{\rm Pm\bar 3m}(a)=2a^3.
$$

These two transformations prevent 12 Å / 24 Å LTA convention confusion.

## 5.2 Affine extraction of a target corner from a dry Na/M mixture

When a source provides a dry Na/M mixed composition rather than a pure endmember, the effective target corner is derived in **charge-equivalent space**. Let $q_M$ and $q_{\rm Na}$ be the target and Na charge equivalents in the source normalization, with $q_{\rm tot}=q_M+q_{\rm Na}$. Then

$$
V_M^{(24)}
=\frac{q_{\rm tot}V_{\rm mix}^{(24)}-q_{\rm Na}V_{\rm Na}^{(24)}}{q_M}.
$$

For monovalent alkalis, charge equivalents equal ion counts. For Mg and Ba, the Group-II population contributes twice its ion count.

## 5.3 Placeholder rule

If the package has no bound defensible species-specific dry crystallographic endpoint for a supported species, that species alone is assigned the standard Na-A endpoint:

$$
V_i^{(24)}=V_{\rm Na}^{(24)}.
$$

This substitution acts only on that species' own $w_i$ contribution. A placeholder species never forces an entire mixed composition to use the Na volume.

# 6. Bound endpoint table

| Species | Source basis | Lattice parameter used | $V_i^{(24)}$ (A^3) | Classification |
|---|---|---:|---:|---|
| H+ | no species-specific endpoint currently bound; Na-A substitution | 24.55930 A via Na | 3703.292036 | PLACEHOLDER |
| Li+ | dry Li92Na4-A, 300 K | 23.88056 A | 3391.674204 | EXTRAPOLATED |
| Na+ | dry Na96-A, 300 K | 24.55930 A | 3703.292036 | DIRECT |
| K+ | dry K93Na3-A, 300 K | 24.55825 A | 3702.801747 | EXTRAPOLATED |
| Rb+ | dry Rb76Na20-A, 300 K | 24.55332 A | 3699.875813 | EXTRAPOLATED |
| Cs+ | dry Cs56Na40-A, 300 K | 24.49934 A | 3656.907061 | EXTRAPOLATED |
| Fr+ | no species-specific endpoint currently bound; Na-A substitution | 24.55930 A via Na | 3703.292036 | PLACEHOLDER |
| Be2+ | no species-specific endpoint currently bound; Na-A substitution | 24.55930 A via Na | 3703.292036 | PLACEHOLDER |
| Mg2+ | dehydrated approximately Mg2Na8-A; residual water reported | 24.4244(14) A | 3521.221887 | EXTRAPOLATED |
| Ca2+ | fully exchanged/dehydrated Ca48-A | 24.47 A | 3663.042156 | DIRECT |
| Sr2+ | two dehydrated Fm-3c refinements; nonideal Sr/extra-framework-Al chemistry | 24.674(1), 24.687(1) A -> package $a_{\rm eff}=24.68$ A | 3758.161808 | QUESTIONABLE |
| Ba2+ | fully dehydrated Ba3.5Na5-A, 12-Al pseudocell | 12.267(7) A | 3683.692333 | EXTRAPOLATED |
| Ra2+ | no species-specific endpoint currently bound; Na-A substitution | 24.55930 A via Na | 3703.292036 | PLACEHOLDER |

# 7. Reproducible endpoint derivations

## 7.1 Na

Carey's dry Na96-A conventional cell gives

$$
V_{\rm Na}^{(24)}=\frac{24.55930^3}{4}=3703.292036385714\;\mathrm{A^3}.
$$

The package stores this value at binary floating-point precision.

## 7.2 Li, K, Rb, and Cs

The Carey measurements are conventional dry cells with 96 charge equivalents. Their 24-Al mixed-cell volumes are $a^3/4$, followed by removal of the Na contribution:

$$
V_M^{(24)}=\frac{96V_{\rm mix}^{(24)}-N_{\rm Na}V_{\rm Na}^{(24)}}{N_M},
$$

because these species are monovalent.

The bound source compositions are:

```text
Li92Na4-A   a = 23.88056 A
K93Na3-A    a = 24.55825 A
Rb76Na20-A  a = 24.55332 A
Cs56Na40-A  a = 24.49934 A
```

This yields the Li/K/Rb/Cs endpoint values tabulated above.

## 7.3 Mg

For approximately Mg2Na8-A per 12 framework Al, the charge-equivalent composition is

```text
Mg: 2 ions x +2 = 4 equivalents
Na: 8 ions x +1 = 8 equivalents
total             = 12 equivalents
```

The reported conventional Fm-3c lattice parameter is `24.4244(14) A`, so

$$
V_{\rm mix}^{(24)}=\frac{24.4244^3}{4},
$$

and

$$
V_{\rm Mg}^{(24)}
=\frac{12V_{\rm mix}^{(24)}-8V_{\rm Na}^{(24)}}{4}
=3521.221887\;\mathrm{A^3}.
$$

The source reports a small residual water content despite dehydration. The endpoint is therefore an extrapolated dry-volume model corner, not an unqualified perfectly anhydrous pure-Mg measurement.

## 7.4 Ca

For fully exchanged/dehydrated Ca-A with conventional `a=24.47 A`,

$$
V_{\rm Ca}^{(24)}=\frac{24.47^3}{4}=3663.042156\;\mathrm{A^3}.
$$

No Na-removal extrapolation is applied.

## 7.5 Sr

Pluth and Smith report two relevant dehydrated Fm-3c lattice parameters:

```text
24.674(1) A
24.687(1) A
```

The package preserves its established v7.14 Sr corner by an explicit deterministic source reduction:

$$
a_{\rm midpoint}=\frac{24.674+24.687}{2}=24.6805\;\mathrm A,
$$

$$
a_{\rm eff}=\operatorname{round}(24.6805,2)=24.68\;\mathrm A,
$$

$$
V_{\rm Sr}^{(24)}=\frac{24.68^3}{4}=3758.161808\;\mathrm{A^3}.
$$

`24.68 A` is therefore a **package representative**, not a single directly reported source measurement. The source chemistry is also nonideal for a pure Sr6-A interpretation: the refinements contain approximately 5.7 Sr per 12-Al pseudocell and an extra-framework Al complex. The endpoint remains explicitly QUESTIONABLE.

## 7.6 Ba

For fully dehydrated Ba3.5Na5-A in the 12-Al Pm-3m pseudocell,

$$
V_{\rm mix}^{(24)}=2(12.267)^3.
$$

Its charge-equivalent composition is

```text
Ba: 3.5 ions x +2 = 7 equivalents
Na: 5 ions  x +1 = 5 equivalents
total              = 12 equivalents
```

so

$$
V_{\rm Ba}^{(24)}
=\frac{12V_{\rm mix}^{(24)}-5V_{\rm Na}^{(24)}}{7}
=3683.692333\;\mathrm{A^3}.
$$

The source reports instability of fully Ba-exchanged LTA toward dehydration. This endpoint is therefore an affine model corner inferred from a dry partially exchanged composition, not a claim that a stable dry Ba6-A endmember was measured.

# 8. Placeholder species

For H, Be, Fr, and Ra, the package currently binds no species-specific empirical dry-LTA volume endpoint. Each therefore uses

$$
V_i^{(24)}=V_{\rm Na}^{(24)}=3703.292036385714\;\mathrm{A^3}.
$$

This is a package-provenance statement, not an exhaustive assertion that suitable literature data do not exist. A future evidence review may replace a placeholder only by explicitly binding a defensible source and source-to-model transformation.

# 9. Scientific limitations

The dry model introduced in v7.14 and retained unchanged for framework-compensating counts in v7.15 is an affine empirical volume model for generating starting structures. It does not claim:

- excess/mixing volumes are zero in real materials;
- all endpoint measurements share identical temperature or hydration history;
- extrapolated corners correspond to stable isolable pure endmembers;
- the Sr QUESTIONABLE corner is an ideal pure-Sr measurement;
- the placeholder endpoints are physically equal to Na-A;
- random or maximum-symmetry cation positions are experimentally relaxed equilibrium positions;
- the model applies outside Al24Si24O96.

These limitations constrain physical interpretation but do not alter the defined numerical model.

# 10. Primary empirical sources

1. T. Carey, C. C. Tang, J. A. Hriljac, and P. A. Anderson, "Chemical Control of Thermal Expansion in Cation-Exchanged Zeolite A," *Chemistry of Materials* **26**, 1561-1566 (2014). DOI: `10.1021/cm403312q`.
2. T. Carey, *Chemical control of thermal expansion in zeolites with the LTA topology*, PhD thesis, University of Birmingham (2013), Chapter 3, especially Tables 3.1 and 3.5-3.9. Stable repository record: https://etheses.bham.ac.uk/id/eprint/4258/.
3. J. M. Adams and L. V. C. Rees, "The structure of a dehydrated partially magnesium exchanged zeolite A (~Mg2Na8-A) by neutron profile refinement," *Journal of Solid State Chemistry* **62**, 184-190 (1986). DOI: `10.1016/0022-4596(86)90230-6`.
4. F. Porcher, M. Souhassou, H. Graafsma, A. Puig-Molina, Y. Dusausoy, and C. Lecomte, "Refinement of framework disorder in dehydrated CaA zeolite from single-crystal synchrotron data," *Acta Crystallographica B* **56**, 766-772 (2000). DOI: `10.1107/S0108768100006893`.
5. J. J. Pluth and J. V. Smith, "Crystal structure of dehydrated strontium-exchanged zeolite A. Absence of near-zero-coordinate strontium(2+) ion. Presence of aluminum complex," *Journal of the American Chemical Society* **104**, 6977-6982 (1982). DOI: `10.1021/ja00389a015`.
6. Y. Kim, V. Subramanian, R. L. Firor, and K. Seff, "Four Crystal Structures of Ba_xNa_(12-2x)-A (1 <= x <= 6) Relating to the Instability of Barium-Exchanged Zeolite A Toward Dehydration," *ACS Symposium Series* **135**, 137-153 (1980). DOI: `10.1021/bk-1980-0135.ch007`.

# 11. Historical note

v7.14 replaced the temporary v7.13 input-template-volume fallback with the current evidence-bound endpoint table and charge-equivalent interpolation. The filename is retained for release lineage and compatibility; the semantics above are the current v7.15 dry-volume authority.

## v7.15 occlusion compatibility boundary

The accepted v7.14 dry-volume model is unchanged. In v7.15 its population
argument is explicitly the **framework-compensating** role only:

```
V_ref = (1/24) sum_i z_i N_i^fw V_i^(24)
```

Neither `occlusion_counterion` cations nor Cl cage guests participate in the
automatic interpolation. This preserves the qualified dry reference as an
initialization policy; it is not evidence that salt loading causes zero
physical volume change.

## v7.16 occlusion boundary (added on v7.16 D1 ratification, 2026-10-17)

V_ref above remains the dry reference of the framework-compensating
population. When the v7.16 cage-salt model is active,
`v7_16_molten_salt_occlusion_extension.md` §9 changes the *default* volume
factor:

```
f_V(n_bar) = 1 - (1 - f_ref) * min(n_bar / n_ref, 1)
f_ref = 0.97,   n_ref = 10.5 Cl per pseudo-cell
```

An explicit `--volume-factor` replaces it, so V_target = volume_factor · V_ref
as here. The dry path and the v7.15 legacy model keep the default of 1.0.
