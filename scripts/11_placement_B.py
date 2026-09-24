"""Placement B (diamond shifted by one lattice step): the full analysis at the same level
as placement A, and the comparison of the two (Sec. XI, Fig. 15).

    python scripts/11_placement_B.py [--recompute]      (~5 min with cached data)

Boundary defects exist for placement A with odd n and for placement B with even n.
"""
import argparse, os, math, json
import numpy as np
from _common import *
from aztecdefects.lattice import aztec, tile_center, depth
from aztecdefects import thermo
from aztecdefects.boundary import k1_outer, k2_split_statistics, boundary_tiles
from aztecdefects.pfaffian import tile_probabilities, local_statistics
from aztecdefects.kasteleyn import k1_weights_exact

ap = argparse.ArgumentParser(); ap.add_argument('--recompute', action='store_true'); a = ap.parse_args()
A, B = load_polys()
print('outer-layer diagonal tiles  n: A  B')
for n in range(6, 12):
    print(f'   {n}: {len(boundary_tiles(n, 0))}  {len(boundary_tiles(n, 1))}')
print('half-step  Z_B(2m+1) = 2^(2m+1) Z_B(2m):', all(B[2*m+1] == [2**(2*m+1)*x for x in B[2*m]] for m in range(1, 7) if 2*m+1 in B))

# ---- n = 10 local statistics (placement B)
mpath = os.path.join(DATA, 'B10_k1_marginals.npz')
if a.recompute or not os.path.exists(mpath):
    ls = local_statistics(aztec(10, 1, 0)); c = [int(x) for x in ls['c']]; T = len(ls['tiles'])
    rho1 = np.array([float(int(ls['ce'][1, e])) / c[1] for e in range(T)])
    np.savez(mpath, rho1=rho1, typ=np.array([t[2] == 'B' for t in ls['tiles']]), cen=np.array([tile_center(t) for t in ls['tiles']]))
m10 = np.load(mpath)
comp = json.load(open(os.path.join(DATA, 'placement_comparison_local.json')))
for key in ['A9', 'B10']:
    o = comp[key]; r1, r2 = o['rows'][0], o['rows'][1]
    print(f"{key}: c1/c0={o['c1c0']:.2f} outer(k=1) B {r1['outerB']:.3f} W {r1['outerW']:.3f}; top site B {o['top_B']:.3f} W {o['top_W']:.3f}; "
          f"<d_BW> {r1['d_BW']:.2f} (indep {r1['d_BW_ind']:.2f}); k=2 <d_BB> {r2['d_BB']:.2f} (indep {r2['d_BB_ind']:.2f})")

# ---- exact outer-layer asymptotics, both placements
ns = [20, 40, 80, 120, 200, 400, 800, 1600]; out = {}
for n in ns:
    for ox, nn in [(1, n), (0, n + 1)]:
        r = k1_outer(nn, ox); out[(ox, nn)] = (math.exp(r['ln_c1_over_c0'] + 3 * math.log(nn) - 2 * nn * math.log(2)), r['var_over_n'], nn * r['corr'])
x = np.array([1 / n for n in ns if n >= 200]); Af = np.vstack([np.ones(len(x)), x, x ** 2]).T
for ox, off, lab in [(1, 0, 'B'), (0, 1, 'A')]:
    arr = np.array([out[(ox, n + off)] for n in ns if n >= 200])
    ex = [np.linalg.lstsq(Af, arr[:, i], rcond=None)[0][0] for i in range(3)]
    print(f'{lab}: extrapolated n^3 c1/(c0 4^n) = {ex[0]:.7f} (2/pi {2/math.pi:.7f}), Var/n = {ex[1]:.7f}, n corr = {ex[2]:.4f}')
k2 = {n: k2_split_statistics(n, ox=1) for n in (40, 80, 120, 160, 200)}
for n, r in k2.items():
    print(f"B k=2 n={n}: n^6 c2/(c0 16^n) = {math.exp(r['ln_c2_over_c0']-4*n*math.log(2)+6*math.log(n)):.5f}, n corr(BB) = {r['n_corr_BB']:.3f}")

# ---- exact k=1 (all layers), u = 1, thresholds
ex1 = {}
for n in (20, 40, 60, 80):
    p = os.path.join(DATA, f'k1_exact_B_n{n}.npz')
    if a.recompute or not os.path.exists(p):
        L, cB, cW = k1_weights_exact(n, ox=1, verbose=False); np.savez(p, L=L, cB=cB, cW=cW)
    d = np.load(p); L = d['L']; P = np.exp(L - L.max()); P /= P.sum(); cB = d['cB']
    dB = n - (np.abs(cB[:, 0] - 1) + np.abs(cB[:, 1])); ex1[n] = (L.max() + math.log(np.exp(L - L.max()).sum()), P.sum(1)[dB == 0].sum())
    print(f'B n={n}: ln(c1/c0) = {ex1[n][0]:.4f}, outer-layer probability {ex1[n][1]:.4f}')
fu = {}
for n in (14, 20, 30, 40, 50, 60):
    _, p = tile_probabilities(aztec(n, 1, 0), 1.0); fu[n] = p.sum() / (n * (n + 1))
print('B u=1 diagonal fraction:', {n: round(v, 5) for n, v in fu.items()})
fa = {n: 2 * sum(k * x for k, x in enumerate(A[n])) / sum(A[n]) / (n * (n + 1)) for n in A if n % 2 == 1 and n >= 5}
fb = {n: 2 * sum(k * x for k, x in enumerate(B[n])) / sum(B[n]) / (n * (n + 1)) for n in B if n % 2 == 0 and n >= 4}

# ---- Fig. 15
fig, ax = plt.subplots(2, 3, figsize=(10.5, 6.4)); a_ = ax[0, 0]
R = aztec(10, 1, 0); cen = m10['cen']; typ = m10['typ']; rho1 = m10['rho1']
for (i, j) in R:
    a_.add_patch(plt.Rectangle((i, j), 1, 1, fc='#eeeeee' if (i + j) % 2 == 0 else 'white', ec='#ddd', lw=.3))
for isB, col in [(True, BCOL), (False, WCOL)]:
    k = typ == isB; a_.scatter(cen[k, 0], cen[k, 1], s=900 * rho1[k], c=col, alpha=.8)
a_.set_aspect('equal'); a_.axis('off'); a_.set_title('(a) placement B, $n=10$, $k=1$\n(area $\\propto$ occupation; black B, red W)', fontsize=8)
a_ = ax[0, 1]
for ox, off, m, col, lab in [(0, 1, 'o', BL, 'A, odd $n$'), (1, 0, 's', OR, 'B, even $n$')]:
    a_.plot([1 / (n + off) for n in ns], [out[(ox, n + off)][0] for n in ns], m + '-', color=col, ms=4, label=lab)
a_.axhline(2 / math.pi, color='k', ls='--', lw=.8); a_.set_xlabel('$1/n$'); a_.set_title(r'(b) $n^3(c_1/c_0)/4^n\to2/\pi$'); a_.legend(frameon=False, fontsize=8)
a_ = ax[0, 2]
for ox, off, m, col, lab in [(0, 1, 'o', BL, 'A'), (1, 0, 's', OR, 'B')]:
    a_.plot([1 / (n + off) for n in ns], [8 * out[(ox, n + off)][1] for n in ns], m + '-', color=col, ms=4, label=f'{lab}: $8\\,$Var$/n$')
    a_.plot([1 / (n + off) for n in ns], [out[(ox, n + off)][2] / 4 for n in ns], m + ':', color=col, ms=4, label=f'{lab}: $n\\,$corr$/4$')
a_.axhline(1, color='k', ls='--', lw=.8); a_.set_xlabel('$1/n$'); a_.set_title('(c) edge spread and correlation'); a_.legend(frameon=False, fontsize=7)
a_ = ax[1, 0]
a_.plot([1 / n for n in fa], list(fa.values()), 'o', color=BL, ms=4, label='A odd (exact)'); a_.plot([1 / n for n in fb], list(fb.values()), 's', color=OR, ms=4, label='B even (exact)')
a_.plot([1 / n for n in fu], list(fu.values()), '^', color=AQ, ms=4, label='B even (Pfaffian)'); a_.axhline(0.162636, color='k', ls='--', lw=.8)
a_.set_xlabel('$1/n$'); a_.set_title('(d) diagonal fraction at $u=1$'); a_.legend(frameon=False, fontsize=7)
a_ = ax[1, 1]
for P_, par, col, lab in [(A, 1, BL, 'A'), (B, 0, OR, 'B')]:
    nn = sorted(n for n in P_ if n % 2 == par and n >= 7)
    for k, m in [(1, 'o'), (3, '^')]:
        a_.plot(nn, [n * math.log(2) * thermo.thresholds(P_[n])[1][k - 1] for n in nn], m + '-', color=col, ms=4, label=f'{lab}, $T_{k}$')
a_.set_xlabel('$n$'); a_.set_ylabel(r'$T_k\,n\ln2/\varepsilon_d$'); a_.set_title('(e) threshold temperatures'); a_.legend(frameon=False, fontsize=7)
a_ = ax[1, 2]; tt = np.linspace(0.05, 1.2, 500)
for P_, n, col in [(A, 13, BL), (B, 14, OR), (A, 15, VI), (B, 16, YE)]:
    if n in P_:
        a_.step(tt, [thermo.k_star(P_[n], t) for t in tt], where='post', color=col, label=f"{'A' if n % 2 else 'B'}, $n={n}$")
a_.set_ylim(0, 9); a_.set_xlabel('$k_BT/\\varepsilon_d$'); a_.set_ylabel('$k^*$'); a_.set_title('(f) staircase $k^*(T)$'); a_.legend(frameon=False, fontsize=7)
for x_ in ax.ravel()[1:]:
    clean(x_)
fig.tight_layout(); save(fig, 'fig15_placementB')
