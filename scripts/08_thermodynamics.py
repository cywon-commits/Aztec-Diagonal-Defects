"""Thermodynamics and free-energy minimisation (Sec. VII, Figs. 13-14).

    python scripts/08_thermodynamics.py [--float-n 21 41]    (needs 01 and 03; ~5 minutes)
Temperatures in units of eps_d/k_B (energy of one diagonal tile above a domino).
"""
import argparse, os, math
import numpy as np
from scipy.optimize import brentq
from _common import *
from aztecdefects import thermo
from aztecdefects.lattice import aztec
from aztecdefects.pfaffian import tile_probabilities

ap = argparse.ArgumentParser(); ap.add_argument('--float-n', type=int, nargs='*', default=[21, 41]); a = ap.parse_args()
polys, _ = load_polys(); b = np.load(os.path.join(DATA, 'bulk.npz'))
odd = sorted(n for n in polys if n % 2 == 1 and n >= 5)
print('n   T_1    T_2    T_3    T_4   (units eps_d)   concave')
for n in odd:
    d, T, conc = thermo.thresholds(polys[n]); print(n, ' '.join(f'{t:6.3f}' for t in T[:4]), conc)
Tstar = {n: brentq(lambda t: thermo.k_mean(polys[n], t) - 1, 0.01, 5) for n in odd if n >= 7}
print('T* (<k>=1):', {n: round(t, 4) for n, t in Tstar.items()})
flt = {}
for n in a.float_n:
    rows = []
    for t in np.logspace(-1.2, 1, 24):
        u = math.exp(-2 / t); _, p1 = tile_probabilities(aztec(n), u); _, p2 = tile_probabilities(aztec(n), u, gauge_seed=3)
        k1, k2 = p1.sum() / 2, p2.sum() / 2
        if abs(k1 - k2) < 1e-6 * max(1, k1):
            rows.append((t, k1))
    flt[n] = np.array(rows); print(f'float n={n}: {len(rows)} reliable temperatures', flush=True)

# Fig. 13
c9 = polys[9]; S9 = thermo.entropy(c9); k9 = np.arange(len(S9))
fig, ax = plt.subplots(2, 2, figsize=(10, 7.6)); a_ = ax[0, 0]; a_.plot(k9, S9 - S9[0], 'k-', lw=2, label=r'$\ln(c_k/c_0)$')
for t, cc in zip([0.3, 0.5, 1, 2], [BL, OR, AQ, VI]):
    F = 2 * k9 - t * (S9 - S9[0]); a_.plot(k9, -F / t, '--', color=cc, label=f'$-\\Delta F/k_BT$, $T={t}$'); a_.plot(k9[F.argmin()], -F.min() / t, 'o', color=cc)
a_.set_ylim(-5, 20); a_.set_xlabel('$k$'); a_.legend(frameon=False, fontsize=7.5); a_.set_title('(a) landscape, $n=9$')
a_ = ax[0, 1]; taus = np.logspace(-1.7, 1, 200)
for n, cc in [(9, BL), (11, '#88aee8')]:
    a_.loglog(taus, [thermo.k_mean(polys[n], t) for t in taus], color=cc, label=f'$n={n}$ exact')
for (n, r), cc in zip(flt.items(), [OR, AQ]):
    a_.loglog(r[:, 0], r[:, 1], 'o-', color=cc, ms=3, label=f'$n={n}$')
a_.set_ylim(1e-3, 300); a_.set_xlabel('$k_BT/\\varepsilon_d$'); a_.set_ylabel(r'$\langle k\rangle$'); a_.legend(frameon=False, fontsize=7.5); a_.set_title('(b) mean pair number')
a_ = ax[1, 0]; a_.semilogx(b['tau'], b['c'], color=BL); a2 = a_.twinx(); a2.semilogx(b['tau'], b['s'], color=OR); a2.axhline(0.2915609, color=OR, ls=':')
a_.set_xlabel('$k_BT/\\varepsilon_d$'); a_.set_ylabel('$c/k_B$ per cell', color=BL); a2.set_ylabel('$s/k_B$ per cell', color=OR); a_.set_title('(c) bulk')
a_ = ax[1, 1]; nn = np.arange(5, 202, 2)
a_.loglog(nn, [thermo.onset_temperature_asymptotic(n) for n in nn], color=BL, label='$2\\varepsilon_d/\\ln(c_1/c_0)$')
a_.loglog(nn, 1 / (nn * math.log(2)), '--', color='#888', label='$\\varepsilon_d/(n\\ln2)$')
a_.loglog(list(Tstar), list(Tstar.values()), 'o', color=OR, label=r'exact $\langle k\rangle=1$'); a_.legend(frameon=False, fontsize=8); a_.set_title('(d) onset temperature')
for x in [ax[0, 0], ax[0, 1], ax[1, 1]]:
    clean(x)
fig.tight_layout(); save(fig, 'fig13_thermo')

# Fig. 14
cols = {7: '#9ec5f4', 9: '#5598e7', 11: '#2a78d6', 13: '#1c5cab', 15: '#104281', 17: '#0d366b'}
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={'width_ratios': [1, 1.15, 1]})
for n in [n for n in cols if n in polys]:
    d, T, _ = thermo.thresholds(polys[n]); ax[0].plot(np.arange(1, 9), d[:8], 'o-', color=cols[n], ms=3, label=f'$n={n}$')
ax[0].set_xlabel('$k$'); ax[0].set_ylabel(r'$\ln(c_k/c_{k-1})$'); ax[0].legend(frameon=False, fontsize=6.5, ncol=2); ax[0].set_title('(a) entropy gain')
tt = np.linspace(0.05, 1.2, 600)
for n in [n for n in (9, 13, 17) if n in polys]:
    ax[1].step(tt, [thermo.k_star(polys[n], t) for t in tt], where='post', color=cols[n], label=f'$k^*$, $n={n}$')
    ax[1].plot(tt, [thermo.k_mean(polys[n], t) for t in tt], ':', color=cols[n])
ax[1].set_ylim(0, 10.5); ax[1].set_xlabel('$k_BT/\\varepsilon_d$'); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title('(b) $k^*(T)$ and $\\langle k\\rangle$')
for k, m, cc in [(1, 'o', '#2a78d6'), (2, 's', '#86b6ef'), (3, '^', '#eb6834'), (4, 'v', '#f4a582')]:
    ax[2].plot([n for n in odd if n >= 7], [n * math.log(2) * thermo.thresholds(polys[n])[1][k - 1] for n in odd if n >= 7], m + '-', color=cc, ms=4, label=f'$k={k}$')
ax[2].set_xscale('log'); ax[2].set_ylim(0.8, 7); ax[2].set_xlabel('$n$'); ax[2].set_ylabel(r'$T_k n\ln2/\varepsilon_d$'); ax[2].legend(frameon=False, fontsize=6.5); ax[2].set_title('(c) thresholds')
for x in ax:
    clean(x)
fig.tight_layout(); save(fig, 'fig14_staircase')
