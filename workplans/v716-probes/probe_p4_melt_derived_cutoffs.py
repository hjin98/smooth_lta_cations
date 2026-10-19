"""Probe 4 (not production): pair minima derived from evidence-file §5 melt densities.
Unlike pair: d_MX = (v/2)^(1/3), v = M/(rho N_A) (rock-salt-like spacing).
Like pair: sqrt(2) * mean of the two d_MX (rock-salt like-ion geometry). All scaled by s."""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/user/smooth_lta_cations")
import lta_cation_population_grouped as P
from scipy.optimize import minimize
NA = 0.602214076  # mol^-1 * 1e24 / cm^3->A^3 factor
melt = dict(Li=(42.39, 1.50), Na=(58.44, 1.55), K=(74.55, 1.52), Cs=(168.36, 2.79))
dMX = {m: ((M/(rho*NA))/2) ** (1/3) for m, (M, rho) in melt.items()}
print({k: round(v, 3) for k, v in dMX.items()}, {k: round(np.sqrt(2)*v, 3) for k, v in dMX.items()})
R = dict(Li=0.76, Na=1.02, K=1.38, Cl=1.81, O=1.40, Si=0.26, Al=0.39)
def base(a, b, dClCl):
    if a == "Cl" and b == "Cl": return dClCl
    if "Cl" in (a, b):
        m = b if a == "Cl" else a
        return dMX[m] if m in dMX else None
    if a in dMX and b in dMX: return np.sqrt(2)*(dMX[a]+dMX[b])/2
    return None  # framework: Shannon f=0.85 fallback (no evidence-file data)
def run(npair, s, s_ring=None, trials=12, seed=0, rconf=6.0, xli=0.585):
    s_ring = s if s_ring is None else s_ring
    vmix = sum(x*(melt[m][0]/(melt[m][1]*NA)) for m, x in (("Li", xli), ("K", 1-xli)))
    dClCl = np.sqrt(2)*(vmix/2)**(1/3)
    r = P.generate_lta_cation_structure_from_specs("tests/real_lta_primitive.cif",
            n4="Na:2", n6="Na:16", n8="Na:6", seed=seed)
    st = r.structure; L = st.lattice; M = L.matrix; Minv = np.linalg.inv(M)
    fx_el = [x.specie.symbol for x in st]; fx = np.array([x.coords for x in st])
    centers = [L.get_cartesian_coords(c) for c in ((.25,.25,.25), (.75,.75,.75))]
    nli = round(npair*xli); el = []; cc = []
    for c in centers:
        for e in ["Cl"]*npair + ["Li"]*nli + ["K"]*(npair-nli): el.append(e); cc.append(c)
    cc = np.array(cc)
    def dmin(a, b, ring=False):
        d = base(a, b, dClCl)
        if d is None: return 0.85*(R[a]+R[b])
        return (s_ring if ring else s)*d
    Dmm = np.array([[dmin(a, b) for b in el] for a in el])
    Dmf = np.array([[dmin(a, b, ring=True) for b in fx_el] for a in el])
    def mic(d):
        fr = d @ Minv; fr -= np.round(fr); return fr @ M
    def pen(x):
        X = x.reshape(-1, 3); g = np.zeros_like(X); e = 0.0
        d = mic(X[:, None]-X[None]); n = np.linalg.norm(d, axis=-1)+1e-12
        v = np.maximum(Dmm-n, 0); np.fill_diagonal(v, 0); e += 0.5*(v**2).sum(); g -= 2*((v/n)[..., None]*d).sum(1)
        d2 = mic(X[:, None]-fx[None]); n2 = np.linalg.norm(d2, axis=-1)+1e-12
        v2 = np.maximum(Dmf-n2, 0); e += (v2**2).sum(); g -= 2*((v2/n2)[..., None]*d2).sum(1)
        dc = X-cc; nc = np.linalg.norm(dc, axis=1)+1e-12; vc = np.maximum(nc-rconf, 0)
        e += (vc**2).sum(); g += 2*(vc/nc)[:, None]*dc
        return e, g.ravel()
    rng = np.random.default_rng(seed); ok = 0; best = 1e9
    for t in range(trials):
        u = rng.normal(size=cc.shape); u /= np.linalg.norm(u, axis=1)[:, None]
        x0 = cc + u*(rconf*rng.random(len(cc))**(1/3))[:, None]
        res = minimize(pen, x0.ravel(), jac=True, method="L-BFGS-B", options=dict(maxiter=5000))
        best = min(best, res.fun); ok += res.fun < 1e-6
    return ok, trials, best, dClCl
for s in (0.70, 0.75, 0.80, 0.85, 0.90, 1.0):
    for n in (8, 10, 11, 12):
        ok, t, b, dcc = run(n, s)
        print(f"s={s:.2f} n/alpha={n}: {ok}/{t} best={b:.2e} (ClCl_base={dcc:.3f})", flush=True)
print("-- ring cations relaxed: like-ion vs ring cations at s_ring=0.7")
for s in (0.85, 0.90, 1.0):
    for n in (10, 12):
        ok, t, b, _ = run(n, s, s_ring=0.7)
        print(f"s={s:.2f} s_ring=0.70 n/alpha={n}: {ok}/{t} best={b:.2e}", flush=True)
