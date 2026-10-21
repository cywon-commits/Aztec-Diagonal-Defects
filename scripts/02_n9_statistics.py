"""Complete k-resolved statistics for n = 9 (Sec. IV, Figs. 3-5).

    python scripts/02_n9_statistics.py        (~2 minutes)
"""
import os, collections, math
import numpy as np
from _common import *
import matplotlib.colors as mc
import networkx as nx
from aztecdefects.lattice import aztec, tile_center, depth
from aztecdefects.pfaffian import local_statistics

n = 9; path = os.path.join(DATA, 'n9_local.npz')
if not os.path.exists(path):
    ls = local_statistics(aztec(n))
    np.savez(path, c=np.array([str(x) for x in ls['c']]), ce=ls['ce'].astype(str), cef=ls['cef'].astype(str),
             I=ls['I'], J=ls['J'], tiles=np.array([[t[0][0], t[0][1], t[1][0], t[1][1], t[2] == 'B'] for t in ls['tiles']]))
d = np.load(path)
c = [int(x) for x in d['c']]; kmax = max(k for k in range(len(c)) if c[k] > 0)
ce = np.vectorize(int, otypes=[object])(d['ce']); cef = np.vectorize(int, otypes=[object])(d['cef'])
tiles = [((int(r[0]), int(r[1])), (int(r[2]), int(r[3])), 'B' if r[4] else 'W') for r in d['tiles']]
I, J = d['I'], d['J']; T = len(tiles); typ = np.array([t[2] for t in tiles])
cen = np.array([tile_center(t) for t in tiles]); dep = np.array([depth(x, n) for x in cen])
rho = np.array([[float(ce[k, e]) / float(c[k]) for e in range(T)] for k in range(kmax + 1)])
dij = np.linalg.norm(cen[I] - cen[J], axis=1); kind = np.array([''.join(sorted(typ[i] + typ[j])) for i, j in zip(I, J)])
print(f'c_1 = {c[1]}  (c_1/c_0 = {c[1]/c[0]:.4f}),  k_max = {kmax}, c_kmax = {c[kmax]}')

summary = []; dist = {}
for k in range(1, kmax + 1):
    row = {'k': k, 'mean_depth_B': (rho[k][typ == 'B'] * dep[typ == 'B']).sum() / k}
    w = np.array([float(cef[k, m]) / float(c[k]) for m in range(len(I))])
    for kd in ['BW', 'BB']:
        sel = kind == kd
        if w[sel].sum() == 0:
            continue
        P = w[sel] / w[sel].sum(); Q = rho[k][I[sel]] * rho[k][J[sel]]; Q /= Q.sum()
        row['d_' + kd] = (P * dij[sel]).sum(); row['d_' + kd + '_ind'] = (Q * dij[sel]).sum()
        agg = collections.defaultdict(lambda: [0., 0.])
        for ds, pp, qq in zip(np.round(dij[sel] ** 2).astype(int), P, Q):
            agg[ds][0] += pp; agg[ds][1] += qq
        dist[(k, kd)] = agg
    summary.append(row)
r1 = summary[0]; r2 = summary[1]
print(f"k=1: outer-layer prob {rho[1][(typ=='B')&(dep==0)].sum():.3f}, <d_BW> {r1['d_BW']:.2f} vs indep {r1['d_BW_ind']:.2f}")
print(f"k=2: <d_BB> {r2['d_BB']:.2f} vs indep {r2['d_BB_ind']:.2f}")
g = lambda k, kd, d2: dist[(k, kd)][d2][0] / dist[(k, kd)][d2][1]
print(f"k=1: g_BW(sqrt2) = {g(1, 'BW', 2):.2f};  k=2: g_BB(2) = {g(2, 'BB', 4):.3f}")
print(f"k=1: mean depth of B {r1['mean_depth_B']:.2f} (uniform {dep[typ == 'B'].mean():.2f});  "
      f"depth crosses the uniform value at k = {next(r['k'] for r in summary if r['mean_depth_B'] > dep[typ == 'B'].mean())}")

# ---------------- Fig. 3 : k = 1
R = aztec(n); Bt = [k for k in range(T) if typ[k] == 'B']; Wt = [k for k in range(T) if typ[k] == 'W']
P1 = np.zeros((len(Bt), len(Wt))); bpos = {e: a for a, e in enumerate(Bt)}; wpos = {e: a for a, e in enumerate(Wt)}
for m, (i, j) in enumerate(zip(I, J)):
    if kind[m] == 'BW':
        b, w = (i, j) if typ[i] == 'B' else (j, i); P1[bpos[b], wpos[w]] = float(cef[1, m])
P1 /= P1.sum(); pB = P1.sum(1); pW = P1.sum(0)
print(f'k=1: maximal B occupation {pB.max():.3f}, most probable pair {P1.max():.3f}')
seq = mc.LinearSegmentedColormap.from_list('b', ['#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
def base(ax):
    for (i, j) in R:
        ax.add_patch(plt.Rectangle((i, j), 1, 1, fc='#eeeeee' if (i + j) % 2 == 0 else 'white', ec='#cccccc', lw=.4))
    ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-n - .5, n + .5); ax.set_ylim(-n - .5, n + .5)
fig, axs = plt.subplots(2, 2, figsize=(12, 12)); nm = mc.Normalize(0, max(pB.max(), pW.max()))
for ax, S_, p, lab in [(axs[0, 0], Bt, pB, '(a) B tile'), (axs[0, 1], Wt, pW, '(b) W tile')]:
    base(ax)
    for e, pp in zip(S_, p):
        a_, b_, _ = tiles[e]; ax.plot([a_[0] + .5, b_[0] + .5], [a_[1] + .5, b_[1] + .5], color=seq(nm(pp)), lw=7, solid_capstyle='round')
        if pp > 0.02:
            ax.text((a_[0] + b_[0]) / 2 + .5, (a_[1] + b_[1]) / 2 + 1.05, f'{pp:.3f}', ha='center', fontsize=8)
    ax.set_title(f'{lab} occupation, $n=9$, $k=1$')
cax = fig.add_axes([0.93, 0.56, 0.012, 0.36]); fig.colorbar(plt.cm.ScalarMappable(norm=nm, cmap=seq), cax=cax, label='probability')
ax = axs[1, 0]; ks = sorted(dist[(1, 'BW')]); x = np.arange(len(ks))
ax.bar(x - .2, [dist[(1, 'BW')][k][0] for k in ks], .38, color=BL, label='actual'); ax.bar(x + .2, [dist[(1, 'BW')][k][1] for k in ks], .38, color=OR, label='product of marginals')
ax.set_xticks(x); ax.set_xticklabels([f'{math.sqrt(k):.2f}' for k in ks], rotation=60, fontsize=8); ax.set_xlabel('B–W center distance $d$'); ax.set_ylabel('probability'); ax.set_title('(c) B–W distance'); ax.legend(frameon=False); clean(ax)
ax = axs[1, 1]; a_, b_ = np.unravel_index(P1.argmax(), P1.shape); be, we = tiles[Bt[a_]], tiles[Wt[b_]]
rest = R - {be[0], be[1], we[0], we[1]}; G = nx.Graph()
for (i, j) in rest:
    for q in [(i + 1, j), (i, j + 1)]:
        if q in rest:
            G.add_edge((i, j), q)
Mt = nx.bipartite.hopcroft_karp_matching(G, top_nodes=[v for v in G if sum(v) % 2 == 0]); base(ax); seen = set()
for u_, v_ in Mt.items():
    if (v_, u_) in seen:
        continue
    seen.add((u_, v_)); x0, y0 = min(u_[0], v_[0]), min(u_[1], v_[1])
    ax.add_patch(plt.Rectangle((x0 + .1, y0 + .1), abs(u_[0] - v_[0]) + .8, abs(u_[1] - v_[1]) + .8, fc='none', ec='#555', lw=1))
for t_, col in [(be, BCOL), (we, WCOL)]:
    ax.plot([t_[0][0] + .5, t_[1][0] + .5], [t_[0][1] + .5, t_[1][1] + .5], color=col, lw=8, solid_capstyle='round')
ax.set_title(f'(d) most probable pair ($P={P1.max():.3f}$)')
fig.subplots_adjust(left=.05, right=.91, top=.96, bottom=.06, wspace=.08, hspace=.12); save(fig, 'fig3_n9k1')

# ---------------- Fig. 4 : all k
fig, ax = plt.subplots(2, 2, figsize=(13, 10)); K = [r['k'] for r in summary]
ax[0, 0].semilogy(range(kmax + 1), [c[k] for k in range(kmax + 1)], 'o-', color=BL); ax[0, 0].set_xlabel('$k$'); ax[0, 0].set_ylabel('$c_k$'); ax[0, 0].set_title('(a) exact counts')
ax[0, 1].plot(K, [r['mean_depth_B'] for r in summary], 'o-', color=BL); ax[0, 1].axhline(dep[typ == 'B'].mean(), ls='--', color='#888'); ax[0, 1].set_xlabel('$k$'); ax[0, 1].set_ylabel('mean depth'); ax[0, 1].set_title('(b) mean depth of B tiles')
ax[1, 0].plot(K, [r['d_BW'] for r in summary], 'o-', color=BL, label='B–W'); ax[1, 0].plot(K, [r['d_BW_ind'] for r in summary], '--', color=BL)
K2 = [r['k'] for r in summary if 'd_BB' in r]
ax[1, 0].plot(K2, [r['d_BB'] for r in summary if 'd_BB' in r], 'o-', color=OR, label='B–B'); ax[1, 0].plot(K2, [r['d_BB_ind'] for r in summary if 'd_BB' in r], '--', color=OR)
ax[1, 0].legend(frameon=False); ax[1, 0].set_xlabel('$k$'); ax[1, 0].set_ylabel('mean distance'); ax[1, 0].set_title('(c) mean pair distance (dashed: independent)')
for k_, cc in zip([2, 5, 10, 20], [BL, OR, AQ, VI]):
    for kd, ls_ in [('BW', '-'), ('BB', ':')]:
        rows = sorted(dist[(k_, kd)].items()); rows = [(math.sqrt(a), v[0] / v[1]) for a, v in rows if v[1] > 1e-4]
        ax[1, 1].plot([r[0] for r in rows], [r[1] for r in rows], ls_, marker='o', ms=3, color=cc, label=f'$k={k_}$ {kd[0]}–{kd[1]}')
ax[1, 1].set_yscale('log'); ax[1, 1].axhline(1, color='#888'); ax[1, 1].legend(frameon=False, ncol=2, fontsize=7); ax[1, 1].set_xlabel('distance $d$'); ax[1, 1].set_ylabel('$g(d)$'); ax[1, 1].set_title('(d) pair correlation')
for a_ in ax.ravel():
    clean(a_)
fig.tight_layout(); save(fig, 'fig4_allk')

# ---------------- Fig. 5 : occupancy maps
ks = [1, 2, 4, 8, 12, 18, 27, 36]; fig, axs = plt.subplots(2, 4, figsize=(7, 4.1)); nm = mc.Normalize(0, 1)
cm = mc.LinearSegmentedColormap.from_list('b', ['#e3edfb', '#86b6ef', '#2a78d6', '#0d366b'])
for a_, k in zip(axs.ravel(), ks):
    for (i, j) in R:
        a_.add_patch(plt.Rectangle((i, j), 1, 1, fc='#f2f2f2' if (i + j) % 2 == 0 else 'white', ec='#dddddd', lw=.2))
    for e in Bt:
        p_, q_, _ = tiles[e]; a_.plot([p_[0] + .5, q_[0] + .5], [p_[1] + .5, q_[1] + .5], color=cm(nm(rho[k, e])), lw=2.6, solid_capstyle='round')
    a_.set_aspect('equal'); a_.axis('off'); a_.set_xlim(-n - .3, n + .3); a_.set_ylim(-n - .3, n + .3); a_.set_title(f'$k={k}$')
fig.subplots_adjust(left=.01, right=.99, top=.93, bottom=.12, wspace=.03, hspace=.15)
cax = fig.add_axes([0.3, 0.05, 0.4, 0.025]); fig.colorbar(plt.cm.ScalarMappable(norm=nm, cmap=cm), cax=cax, orientation='horizontal', label='occupation probability')
save(fig, 'fig5_occupancy')
