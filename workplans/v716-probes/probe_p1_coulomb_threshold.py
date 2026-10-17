"""Feasibility probe (not production): pack n MCl per alpha cage under
hard-core f*(ri+rj) and same-sign d0*sqrt(zi zj), Packmol-style penalty minimization."""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/user/smooth_lta_cations")
import lta_cation_population_grouped as P
from scipy.optimize import minimize
R = dict(Li=0.76, Na=1.02, K=1.38, Cl=1.81, O=1.40, Si=0.26, Al=0.39)  # Shannon CN6 (Si/Al CN4)
Z = dict(Li=1, Na=1, K=1, Cl=-1, O=0, Si=0, Al=0)  # framework: radius rule only
def run(npair, f, d0, frac_li=0.585, trials=12, seed=0, rconf=6.0):
    r = P.generate_lta_cation_structure_from_specs("tests/real_lta_primitive.cif",
            n4="Na:2", n6="Na:16", n8="Na:6", seed=seed)
    s = r.structure; L = s.lattice; M = L.matrix
    fixed = [(x.specie.symbol, x.coords) for x in s]
    centers = [L.get_cartesian_coords(c) for c in ((.25,.25,.25), (.75,.75,.75))]
    nli = round(npair*frac_li); mob = []
    for c in centers:
        mob += [("Cl", c)]*npair + [("Li", c)]*nli + [("K", c)]*(npair-nli)
    el = [m[0] for m in mob]; cc = np.array([m[1] for m in mob])
    Minv = np.linalg.inv(M)
    def mic(d):
        fr = d @ Minv; fr -= np.round(fr); return fr @ M
    fx_el = [a for a, _ in fixed]; fx = np.array([b for _, b in fixed])
    def dmin(a, b):
        d = f*(R[a]+R[b])
        if Z[a]*Z[b] > 0: d = max(d, d0*np.sqrt(Z[a]*Z[b]))
        return d
    Dmm = np.array([[dmin(a, b) for b in el] for a in el])
    Dmf = np.array([[dmin(a, b) for b in fx_el] for a in el])
    iu = np.triu_indices(len(el), 1)
    def pen(x):
        X = x.reshape(-1, 3); g = np.zeros_like(X); e = 0.0
        d = mic(X[:, None, :]-X[None, :, :]); n = np.linalg.norm(d, axis=-1)+1e-12
        v = np.maximum(Dmm-n, 0); np.fill_diagonal(v, 0); e += 0.5*(v**2).sum()
        g -= 2*((v/n)[..., None]*d).sum(1)
        d2 = mic(X[:, None, :]-fx[None, :, :]); n2 = np.linalg.norm(d2, axis=-1)+1e-12
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
    return ok, trials, best
if __name__ == "__main__":
    for f in (0.8, 0.9):
        for d0 in (3.0, 3.5, 4.0):
            for n in (6, 8, 10, 11, 12):
                ok, t, b = run(n, f, d0)
                print(f"f={f} d0={d0} n/alpha={n}: feasible {ok}/{t} best_pen={b:.2e}", flush=True)
