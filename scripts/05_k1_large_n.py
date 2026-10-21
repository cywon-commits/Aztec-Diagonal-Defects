"""Exact single-pair (k=1) weights up to large n by multimodular banded elimination (Sec. XI, Fig. 14; file fig9).

    python scripts/05_k1_large_n.py --n 21 41 61 81 [101 121]
    (n = 81: ~2 min, 101: ~5 min, 121: ~12 min)
"""
import argparse, os, math
import numpy as np
from _common import *
from aztecdefects.kasteleyn import k1_weights_exact

ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, nargs='+', default=[21, 41, 61, 81]); a = ap.parse_args()
for n in a.n:
    path = os.path.join(DATA, f'k1_exact_n{n}.npz')
    if not os.path.exists(path):
        L, cB, cW = k1_weights_exact(n); np.savez(path, L=L, cB=cB, cW=cW)
done = sorted(int(f[10:-4]) for f in os.listdir(DATA) if f.startswith('k1_exact_n'))
for n in done:
    d = np.load(os.path.join(DATA, f'k1_exact_n{n}.npz')); L, cB, cW = d['L'], d['cB'], d['cW']
    P = np.exp(L - L.max()); P /= P.sum(); pB = P.sum(1); dB = n - (np.abs(cB[:, 0]) + np.abs(cB[:, 1]))
    ne = (dB == 0) & (cB[:, 0] > 0) & (cB[:, 1] > 0); x = cB[:, 0][ne] / n; w = pB[ne] / pB[ne].sum(); mu = (w * x).sum()
    lnr = L.max() + math.log(np.exp(L - L.max()).sum())
    print(f'n={n}: ln(c1/c0)={lnr:.4f}  outer-layer {pB[dB==0].sum():.4f}  sqrt(n)*sigma(x/n)={math.sqrt(n*(w*(x-mu)**2).sum()):.4f}')
nmax = max(done); d = np.load(os.path.join(DATA, f'k1_exact_n{nmax}.npz')); L = d['L']; P = np.exp(L - L.max()); P /= P.sum()
fig, ax = plt.subplots(1, 2, figsize=(7, 3.2), gridspec_kw={'width_ratios': [1, 1.2]}); a_ = ax[0]
a_.plot([0, nmax, 0, -nmax, 0], [nmax, 0, -nmax, 0, nmax], color='#999', lw=.8); th = np.linspace(0, 2 * np.pi, 300)
a_.plot(nmax / math.sqrt(2) * np.cos(th), nmax / math.sqrt(2) * np.sin(th), ':', color='#999', lw=.8)
for c_, p_, col in [(d['cB'], P.sum(1), BL), (d['cW'], P.sum(0), OR)]:
    k = p_ > 1e-4; a_.scatter(c_[k, 0], c_[k, 1], c=col, s=4 + 400 * p_[k])
a_.text(-.62 * nmax, .62 * nmax, 'W', color=OR, ha='center', fontsize=10); a_.text(.62 * nmax, .62 * nmax, 'B', color=BL, ha='center', fontsize=10)
a_.set_aspect('equal'); a_.axis('off'); a_.set_title(f'(a) marginals, $n={nmax}$, $k=1$')
for n, col in zip([m for m in done if m >= 41][-5:], [BL, OR, AQ, VI, '#111111']):
    dd = np.load(os.path.join(DATA, f'k1_exact_n{n}.npz')); L_ = dd['L']; cB = dd['cB']; P_ = np.exp(L_ - L_.max()); P_ /= P_.sum(); pB = P_.sum(1)
    dB = n - (np.abs(cB[:, 0]) + np.abs(cB[:, 1])); ne = (dB == 0) & (cB[:, 0] > 0) & (cB[:, 1] > 0)
    z = (cB[:, 0][ne] - n / 2) / math.sqrt(n); w = pB[ne] / pB[ne].sum(); o = np.argsort(z)
    ax[1].plot(z[o], w[o] / (2 / math.sqrt(n)), color=col, lw=1.5, label=f'$n={n}$')
zz = np.linspace(-1.6, 1.6, 200); s = 1 / (2 * math.sqrt(2)); ax[1].plot(zz, np.exp(-zz ** 2 / (2 * s * s)) / (s * math.sqrt(2 * math.pi)), 'k--', lw=1.2, label=r'Gaussian $\sigma^2=1/8$')
ax[1].set_xlim(-1.6, 1.6); ax[1].set_xlabel(r'$(x-n/2)/\sqrt{n}$ along the NE edge'); ax[1].set_ylabel('density'); ax[1].legend(frameon=False, fontsize=8); clean(ax[1]); ax[1].set_title('(b) B position along its edge')
fig.tight_layout(); save(fig, 'fig9_largen')
