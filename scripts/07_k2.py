"""Two defect pairs from the determinant formula, Theorem 2 (Sec. X, Fig. 12).

    python scripts/07_k2.py [--n 41 81 121 161 201]     (n = 201: ~1 minute)
"""
import argparse, math
import numpy as np
from _common import *
from aztecdefects.boundary import k2_split_statistics, k1_outer_total

ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, nargs='+', default=[41, 81, 121, 161, 201]); a = ap.parse_args()
rows = []
for n in a.n:
    r = k2_split_statistics(n); k1 = k1_outer_total(n)
    C2 = math.exp(r['ln_c2_over_c0'] - 4 * n * math.log(2) + 6 * math.log(n)); kap = math.exp(r['ln_c2_over_c0'] - 2 * k1['ln_c1_over_c0'])
    rows.append((n, C2, kap, r['n_corr_BB'], -r['n_corr_BW']))
    print(f"n={n}: n^6 c2/(c0 16^n)={C2:.5f}  c2 c0/c1^2={kap:.5f}  Var/n={r['var_over_n']:.5f}  n corr(BB)={r['n_corr_BB']:.3f}  n|corr(BW)|={-r['n_corr_BW']:.3f}", flush=True)
rows = np.array(rows); x = 1 / rows[:, 0]
if len(rows) >= 4:
    A = np.vstack([x ** p for p in range(4)]).T
    print(f'extrapolated C2 = {np.linalg.lstsq(A, rows[:, 1], rcond=None)[0][0]:.5f}  ((2/pi)^2 = {(2/math.pi)**2:.5f});  '
          f'c2 c0/c1^2 -> {np.linalg.lstsq(A, rows[:, 2], rcond=None)[0][0]:.4f}')
fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.6))
ax[0].plot(x, rows[:, 1], 'o-', color=BL, label=r'$n^6 c_2/(c_0 16^n)$'); ax[0].plot(x, rows[:, 2] * (2 / math.pi) ** 2, 's--', color=OR, label=r'$(2/\pi)^2 c_2c_0/c_1^2$')
ax[0].axhline((2 / math.pi) ** 2, color='k', ls=':'); ax[0].legend(frameon=False, fontsize=7); ax[0].set_title('(a) weight of two pairs')
ax[1].plot(x, rows[:, 4], 's-', color=OR, label=r'$n|$corr(B,W)$|$'); ax[1].plot(x, rows[:, 3], '^-', color=BL, label='$n$ corr(B,B)')
ax[1].axhline(4, color='#bbb', ls=':'); ax[1].axhline(2, color='#bbb', ls=':'); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title('(b) residual correlations')
for x_ in ax:
    x_.set_xlabel('$1/n$'); clean(x_)
fig.tight_layout(); save(fig, 'fig12_k2')
