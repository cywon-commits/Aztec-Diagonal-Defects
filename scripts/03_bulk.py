"""Bulk properties of the periodic model (Sec. V, Fig. 6) and bulk thermodynamics data.

    python scripts/03_bulk.py      (~5 minutes)
"""
import os, math
import numpy as np
from _common import *
from aztecdefects import bulk

pw = bulk.bulk_pair_weight()
for d, (r, v) in pw.items():
    sl = np.polyfit(np.log(r[4:]), np.log(v[4:]), 1)[0]; print(f'direction {d}: fitted exponent {sl:.3f}')
ts = np.logspace(-2, 0.3, 13); gap = np.array([bulk.spectral_min(t) for t in ts])
print('min|P| / t^4 :', (gap / ts ** 4).round(6))

mats = bulk.momentum_matrices(500)
print('phi(0) =', bulk.thermo(0.0, mats=mats)[0], ' G/pi =', 0.915965594177219 / math.pi)
f1 = bulk.diagonal_fraction(1.0, mats=mats); print('bulk diagonal fraction at u=1:', f1)
us = np.logspace(-3, 2, 16); fu = np.array([bulk.diagonal_fraction(u, mats=mats) for u in us])
taus = np.logspace(np.log10(0.25), 1, 40); th = np.array([bulk.thermo(math.exp(-2 / t), mats=mats) for t in taus])
phi, rho = th[:, 0], th[:, 1]; s = phi + (2 / taus) * rho; cv = np.gradient(2 * rho, taus)
print(f'specific-heat maximum c = {cv.max():.4f} k_B/cell at T = {taus[cv.argmax()]:.3f} eps_d')
np.savez(os.path.join(DATA, 'bulk.npz'), ts=ts, gap=gap, us=us, fu=fu, f1=f1, tau=taus, phi=phi, rho=rho, s=s, c=cv,
         **{f'r{d[0]}{d[1]}': pw[d][0] for d in pw}, **{f'v{d[0]}{d[1]}': pw[d][1] for d in pw})

fig, ax = plt.subplots(1, 2, figsize=(7, 2.9)); a = ax[0]
for d, c_, lab in [((1, 0), BL, 'axis direction'), ((1, 1), OR, 'diagonal direction')]:
    a.loglog(*pw[d], 'o', color=c_, ms=4, label=lab)
r = np.logspace(0.3, 2, 20); a.loglog(r, 0.9 * pw[(1, 0)][1][-1] * (pw[(1, 0)][0][-1] / r) ** 2, 'k--', lw=1, label=r'$\propto r^{-2}$')
a.set_xlabel('B–W separation $r$'); a.set_ylabel(r'pair weight $\mathcal{W}_{\rm BW}(r)$'); a.legend(frameon=False, fontsize=8); a.set_title('(a) bulk B–W interaction')
a = ax[1]; a.loglog(ts, gap, 'o', color=BL, ms=4, label='min over torus'); a.loglog(ts, ts ** 4, 'k--', lw=1, label='$t^4$')
a.set_xlabel('diagonal-tile weight $t$ ($u=t^2$)'); a.set_ylabel(r'$\min|P|$'); a.legend(frameon=False, fontsize=8); a.set_title('(b) spectral gap')
for x in ax:
    clean(x)
fig.tight_layout(); save(fig, 'fig6_bulk')
