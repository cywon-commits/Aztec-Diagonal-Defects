"""Local statistics for all k of A_9 (placement A) and A_8, A_10 (placement B), Sec. XIII.

Writes data/placement_comparison_local.json, which scripts/11_placement_B.py reads.

    python scripts/13_placement_local.py      (needs data/n9_local.npz from script 02; ~15 min)
"""
import os, json, time
import numpy as np
from _common import *
from aztecdefects.lattice import aztec, tile_center, depth
from aztecdefects.pfaffian import local_statistics


def analyze(c, ce, cef, tiles, I, J, n, ox):
    c = [int(x) for x in c]; kmax = max(k for k in range(len(c)) if c[k] > 0); T = len(tiles)
    typ = np.array([t[2] for t in tiles]); cen = np.array([tile_center(t) for t in tiles]); dep = np.array([depth(x, n, ox) for x in cen])
    rho = np.array([[float(int(ce[k, e])) / c[k] for e in range(T)] for k in range(kmax + 1)])
    dij = np.linalg.norm(cen[I] - cen[J], axis=1); kind = np.array([''.join(sorted(typ[i] + typ[j])) for i, j in zip(I, J)])
    out = dict(n=n, ox=ox, kmax=kmax, c1c0=c[1] / c[0], nB=int((typ == 'B').sum()), nW=int((typ == 'W').sum()),
               nBouter=int(((typ == 'B') & (dep == 0)).sum()), nWouter=int(((typ == 'W') & (dep == 0)).sum()))
    rows = []
    for k in range(1, kmax + 1):
        w = np.array([float(int(cef[k, m])) / c[k] for m in range(len(I))]); r = {'k': k}
        r['depthB'] = (rho[k][typ == 'B'] * dep[typ == 'B']).sum() / k; r['depthW'] = (rho[k][typ == 'W'] * dep[typ == 'W']).sum() / k
        r['outerB'] = rho[k][(typ == 'B') & (dep == 0)].sum() / k; r['outerW'] = rho[k][(typ == 'W') & (dep == 0)].sum() / k
        for kd in ['BW', 'BB', 'WW']:
            sel = kind == kd
            if w[sel].sum() == 0:
                continue
            P = w[sel] / w[sel].sum(); Q = rho[k][I[sel]] * rho[k][J[sel]]; Q /= Q.sum()
            r['d_' + kd] = (P * dij[sel]).sum(); r['d_' + kd + '_ind'] = (Q * dij[sel]).sum()
            m = dij[sel] < 1.5
            if m.any():
                r['g_' + kd + '_adj'] = P[m].sum() / Q[m].sum()
            m2 = np.abs(dij[sel] - 2) < 0.01
            if m2.any():
                r['g_' + kd + '_2'] = P[m2].sum() / Q[m2].sum()
        rows.append(r)
    out['rows'] = rows; out['uniform_depth'] = float(dep[typ == 'B'].mean())
    out['top_B'] = float(rho[1][typ == 'B'].max()); out['top_W'] = float(rho[1][typ == 'W'].max())
    return out


d9 = np.load(os.path.join(DATA, 'n9_local.npz'), allow_pickle=True)
toint = np.vectorize(int, otypes=[object])
out = {'A9': analyze(d9['c'], toint(d9['ce']), toint(d9['cef']),
                     [((int(r[0]), int(r[1])), (int(r[2]), int(r[3])), 'B' if r[4] else 'W') for r in d9['tiles']], d9['I'], d9['J'], 9, 0)}
for n in (10, 8):
    t = time.time(); ls = local_statistics(aztec(n, 1, 0))
    out[f'B{n}'] = analyze(ls['c'], ls['ce'], ls['cef'], ls['tiles'], ls['I'], ls['J'], n, 1); print(f'B{n}: {time.time() - t:.0f}s', flush=True)
for key, o in out.items():
    r1, r2 = o['rows'][0], o['rows'][1]
    cross = next((r['k'] for r in o['rows'] if r['depthB'] > o['uniform_depth']), None)
    print(f"{key}: kmax {o['kmax']}  top B/W {o['top_B']:.3f}/{o['top_W']:.3f}  k=1 outer B/W {r1['outerB']:.3f}/{r1['outerW']:.3f}  "
          f"<d_BW> {r1['d_BW']:.2f} vs {r1['d_BW_ind']:.2f}  g_adj {r1['g_BW_adj']:.2f}  k=2 <d_BB> {r2['d_BB']:.2f} vs {r2['d_BB_ind']:.2f}  "
          f"g_BB(2) {r2.get('g_BB_2', float('nan')):.3f}  depth crosses uniform at k={cross}")
json.dump(out, open(os.path.join(DATA, 'placement_comparison_local.json'), 'w'), default=float)
