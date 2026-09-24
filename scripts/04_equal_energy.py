"""Equal creation energies, u = 1 (Sec. VI, Fig. 11).

    python scripts/04_equal_energy.py [--nmax 61]     (~3 minutes for nmax 61)
"""
import argparse, os, math
import numpy as np
from _common import *
import matplotlib.colors as mc
from aztecdefects.lattice import aztec, tile_center
from aztecdefects.pfaffian import tile_probabilities

ap = argparse.ArgumentParser(); ap.add_argument('--nmax', type=int, default=61); a = ap.parse_args()
polys, _ = load_polys(); b = np.load(os.path.join(DATA, 'bulk.npz')); finf = float(b['f1'])
exact = {}
for n in sorted(polys):
    c = polys[n]; Z = sum(c); k1 = sum(k * x for k, x in enumerate(c)) / Z; k2 = sum(k * k * x for k, x in enumerate(c)) / Z
    exact[n] = 2 * k1 / (n * (n + 1)); print(f'n={n:2d}  <k>={k1:.4f}  Var k={k2-k1*k1:.4f}  f={exact[n]:.5f}')
num = {}; maps = {}
for n in [x for x in (15, 21, 31, 41, 51, 61, 71, 81) if x <= a.nmax]:
    tiles, p = tile_probabilities(aztec(n), 1.0); _, p2 = tile_probabilities(aztec(n), 1.0, gauge_seed=1)
    num[n] = p.sum() / (n * (n + 1)); maps[n] = (tiles, p)
    print(f'n={n}: f={num[n]:.5f}  gauge difference {np.abs(p-p2).max():.1e}', flush=True)
ns = np.array(sorted(num)); f = np.array([num[n] for n in ns])
co = np.linalg.lstsq(np.vstack([np.ones(len(ns)), 1 / ns, 1 / ns ** 2]).T, f, rcond=None)[0]
print(f'fit f(n) = {co[0]:.5f} + {co[1]:.3f}/n + {co[2]:.2f}/n^2   (bulk {finf:.6f})')
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={'width_ratios': [1.1, 1, 1.1]}); a_ = ax[0]
ne = np.array(sorted(n for n in exact if n >= 3)); fe = np.array([exact[n] for n in ne])
a_.plot(1 / ne[ne % 2 == 1], fe[ne % 2 == 1], 'o', color=BL, ms=4, label='exact, odd $n$')
a_.plot(1 / ne[ne % 2 == 0], fe[ne % 2 == 0], 's', color=OR, ms=4, label='exact, even $n$')
a_.plot(1 / ns, f, '^', color=AQ, ms=4, label='Pfaffian'); a_.axhline(finf, color='k', ls='--', lw=.8)
a_.set_xlabel('$1/n$'); a_.set_ylabel('diagonal fraction $f(n)$'); a_.legend(frameon=False, fontsize=6.5); clean(a_); a_.set_title('(a) $u=1$')
nn = max(maps); tiles, p = maps[nn]; cen = np.array([tile_center(t) for t in tiles])
sc = ax[1].scatter(cen[:, 0], cen[:, 1], c=p, s=1.2, cmap=mc.LinearSegmentedColormap.from_list('b', ['#e3edfb', '#2a78d6', '#0d366b']), vmin=.1, vmax=.3)
ax[1].set_aspect('equal'); ax[1].axis('off'); ax[1].set_title(f'(b) occupancy, $n={nn}$'); fig.colorbar(sc, ax=ax[1], fraction=.045)
ax[2].semilogx(b['us'], b['fu'], 'o-', color=BL, ms=3); ax[2].set_xlabel('$u$'); ax[2].set_ylabel(r'$f_\infty(u)$'); ax[2].set_title('(c) bulk'); clean(ax[2])
fig.tight_layout(); save(fig, 'fig11_equal_energy')
