"""Large-n asymptotics of one pair from the exact boundary formula (Sec. XI, Fig. 15; file fig10).

    python scripts/06_boundary_asymptotics.py [--nmax 1601]     (~1 minute)
"""
import argparse, os, math
import numpy as np
from _common import *
from aztecdefects.boundary import k1_outer_total

ap = argparse.ArgumentParser(); ap.add_argument('--nmax', type=int, default=1601); a = ap.parse_args()
ns = [n for n in (21, 41, 61, 81, 101, 121, 201, 401, 801, 1201, 1601) if n <= a.nmax]; R = {}
for n in ns:
    R[n] = k1_outer_total(n)
    C = math.exp(R[n]['ln_c1_over_c0'] + 3 * math.log(n) - 2 * n * math.log(2))
    print(f"n={n:5d}  n^3 c1/(c0 4^n)={C:.6f}  Var/n={R[n]['var_over_n']:.6f}  n*corr={n*R[n]['corr']:.4f}", flush=True)
x = np.array([1 / n for n in ns if n >= 201]); A = np.vstack([np.ones(len(x)), x, x ** 2]).T
C = np.array([math.exp(R[n]['ln_c1_over_c0'] + 3 * math.log(n) - 2 * n * math.log(2)) for n in ns if n >= 201])
V = np.array([R[n]['var_over_n'] for n in ns if n >= 201])
if len(x) >= 3:
    print(f'extrapolated prefactor {np.linalg.lstsq(A, C, rcond=None)[0][0]:.7f}   (2/pi = {2/math.pi:.7f})')
    print(f'extrapolated Var/n     {np.linalg.lstsq(A, V, rcond=None)[0][0]:.7f}   (1/8)')
ex = {}
for f in os.listdir(DATA):
    if f.startswith('k1_exact_n'):
        n = int(f[10:-4]); d = np.load(os.path.join(DATA, f)); L = d['L']; ex[n] = L.max() + math.log(np.exp(L - L.max()).sum())
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5)); nn = np.array(ns)
ax[0].plot(1 / nn, [math.exp(R[n]['ln_c1_over_c0'] + 3 * math.log(n) - 2 * n * math.log(2)) for n in ns], 'o-', color=BL, ms=4, label='outer-layer formula')
if ex:
    ke = sorted(ex); ax[0].plot([1 / n for n in ke], [math.exp(ex[n] + 3 * math.log(n) - 2 * n * math.log(2)) for n in ke], 's', color=OR, label='multimodular (6 layers)')
ax[0].axhline(2 / math.pi, color='#555', ls='--'); ax[0].set_title(r'(a) $\to 2/\pi$'); ax[0].legend(frameon=False, fontsize=7)
ax[1].plot(1 / nn, [R[n]['var_over_n'] for n in ns], 'o-', color=BL, ms=4); ax[1].axhline(1 / 8, color='#555', ls='--'); ax[1].set_title(r'(b) Var$/n\to1/8$')
ax[2].plot(1 / nn, [n * R[n]['corr'] for n in ns], 'o-', color=BL, ms=4); ax[2].axhline(4, color='#555', ls='--'); ax[2].set_title(r'(c) $n\,$corr$\to4$')
for x_ in ax:
    x_.set_xlabel('$1/n$'); clean(x_)
fig.tight_layout(); save(fig, 'fig10_asymptotics')
