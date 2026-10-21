"""Two defect pairs from the determinant formula, Theorem 2 (Sec. XII, Fig. 16; file fig12).

    python scripts/07_k2.py [--n 41 81 121 161 201]     (n = 201: ~1 minute)
"""
import argparse, math
import numpy as np
from _common import *
from aztecdefects.boundary import k2_split_statistics, k1_outer_total, k2_all_sectors

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
    lnw = np.log(rows[:, 1]) - 6 * np.log(rows[:, 0])     # ln[(c2/c0)/16^n]
    print(f'fitted exponent of n with growth rate 16 fixed: {np.polyfit(np.log(rows[:, 0]), lnw, 1)[0]:.2f}')
same = {}
for n in [21, 41, 61, 81, 121]:
    lnc2, ps = k2_all_sectors(n); same[n] = n * n * ps
    print(f'n={n}: n^2 P(both B tiles on one edge) = {same[n]:.3f}   all-sector n^6 c2/(c0 16^n) = {math.exp(lnc2 - 4 * n * math.log(2) + 6 * math.log(n)):.5f}', flush=True)
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={'width_ratios': [1, 1.1, 1.1]})
a_ = ax[0]; m = 10
a_.plot([0, m, 0, -m, 0], [m, 0, -m, 0, m], 'k', lw=1); th = np.linspace(0, 2 * np.pi, 200)
a_.plot(m / np.sqrt(2) * np.cos(th), m / np.sqrt(2) * np.sin(th), ':', color='#888', lw=.8)
pts = {'NE': (m / 2, m / 2, 'B', '+2', BCOL), 'SW': (-m / 2, -m / 2, 'B', '+2', BCOL), 'SE': (m / 2, -m / 2, 'W', '\u22122', WCOL), 'NW': (-m / 2, m / 2, 'W', '\u22122', WCOL)}
for k_, (x_, y_, t_, q_, c_) in pts.items():
    a_.plot(x_, y_, 'o', ms=11, color=c_); a_.text(x_ * 1.62, y_ * 1.62, f'{t_} ({q_})', ha='center', va='center', fontsize=7.5, color=c_)
for u_, v_ in [('NE', 'SE'), ('SE', 'SW'), ('SW', 'NW'), ('NW', 'NE')]:
    a_.annotate('', xy=pts[v_][:2], xytext=pts[u_][:2], arrowprops=dict(arrowstyle='<->', color=AQ, lw=1, shrinkA=7, shrinkB=7))
for u_, v_ in [('NE', 'SW'), ('SE', 'NW')]:
    a_.plot([pts[u_][0], pts[v_][0]], [pts[u_][1], pts[v_][1]], '--', color='#aaa', lw=.8)
a_.set_aspect('equal'); a_.axis('off'); a_.set_title('(a) dominant $k=2$ configuration', fontsize=9); a_.set_xlim(-m * 1.25, m * 1.25); a_.set_ylim(-m * 1.25, m * 1.25)
b_ = ax[1]
b_.plot(x, rows[:, 1], 'o-', color=BL, ms=4, label=r'$n^6 (c_2/c_0)/16^n$'); b_.plot(x, rows[:, 2] * (2 / math.pi) ** 2, 's--', color=OR, ms=4, label=r'$(2/\pi)^2\,c_2c_0/c_1^2$')
b_.axhline((2 / math.pi) ** 2, color='k', ls=':', lw=.9); b_.text(0.019, (2 / math.pi) ** 2 - 0.012, r'$(2/\pi)^2$', fontsize=8)
b_.legend(frameon=False, fontsize=7, loc='lower left'); b_.set_xlim(0, 0.026); b_.set_title('(b) weight of two pairs', fontsize=9)
c_ = ax[2]; nn = sorted(same)
c_.plot([1 / n for n in nn], [same[n] for n in nn], 'o-', color=AQ, ms=4, label=r'$n^2 P(\mathrm{B,B\ same\ edge})$')
c_.plot(x, rows[:, 4], 's-', color=OR, ms=4, label=r'$n\,|\mathrm{corr}(\mathrm{B},\mathrm{W})|$, adjacent')
c_.plot(x, rows[:, 3], '^-', color=BL, ms=4, label=r'$n\,\mathrm{corr}(\mathrm{B},\mathrm{B})$, antipodal')
for yv in [2, 3, 4]:
    c_.axhline(yv, color='#bbb', lw=.6, ls=':')
c_.legend(frameon=False, fontsize=6.5); c_.set_xlim(0, 0.05); c_.set_ylim(1.5, 4.3); c_.set_title('(c) residual interactions', fontsize=9)
for x_ in (b_, c_):
    x_.set_xlabel('$1/n$'); clean(x_)
fig.tight_layout(); save(fig, 'fig12_k2')
