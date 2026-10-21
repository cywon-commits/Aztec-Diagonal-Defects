"""Independent checks of every ingredient of Theorems 1 and 2 (Secs. IX-XI).

    python scripts/10_verify_proofs.py [--nmax 10]    (~2 minutes)
"""
import argparse, itertools
from _common import *
from aztecdefects.lattice import aztec, boundary_edges
from aztecdefects.kasteleyn import count_tilings as M
from aztecdefects import boundary, shuffling

ap = argparse.ArgumentParser(); ap.add_argument('--nmax', type=int, default=10); a = ap.parse_args()
print('renewal gadget Z_old/Z_new =', shuffling.renewal_ratio(3, ((0, 0), (1, 0), (1, -1), (0, -1))), '(expected 2)')
for n in range(2, a.nmax + 1):
    R = aztec(n); E = boundary_edges(n); NE, SE = E['NE'], E['SE']; P = 2 ** (n * (n - 1) // 2)
    F = lambda nn, i, j: M(aztec(nn) - {boundary_edges(nn)['NE'][i], boundary_edges(nn)['SE'][j]})
    thm1 = all(F(n, i, j) == P * boundary.N_sum(n, i, n - 1 - j) for i in range(n) for j in range(n))
    crust = M(R - {NE[0], SE[n - 1]}) == P
    lem2 = all(M(R - {NE[0], NE[i], SE[j], SE[n - 1]}) == F(n - 1, i - 1, j) for i in range(1, n) for j in range(n - 1))
    kuo = all(M(R - {NE[0], SE[n - 1]}) * F(n, i, j) == M(R) * M(R - {NE[0], NE[i], SE[j], SE[n - 1]}) + F(n, 0, j) * F(n, i, n - 1)
              for i in range(1, n) for j in range(n - 1))
    lem4 = all(F(n, i, n - 1) == 2 * sum(F(n - 1, s, i - 1) for s in range(n - 1)) for i in range(1, n))
    print(f'n={n:2d}: Theorem1 {thm1}  crust {crust}  Lemma2 {lem2}  Kuo {kuo}  Lemma4 {lem4}', flush=True)
for n in (7, 9):
    R = aztec(n); E = boundary_edges(n); c0 = 2 ** (n * (n + 1) // 2); T = boundary.boundary_tiles(n)
    Bs = [t for t in T if t[0] in ('NE', 'SW')]; Ws = [t for t in T if t[0] in ('SE', 'NW')]; ok = True; cnt = 0
    for Bp in itertools.combinations(Bs, 2):
        for Wp in itertools.combinations(Ws, 2):
            removed = {E[e][s + d] for (e, s) in Bp + Wp for d in (0, 1)}
            ok &= M(R - removed) * 2 ** (4 * n) == c0 * boundary.det_formula(n, Bp, Wp); cnt += 1
    print(f'Theorem 2, k=2, n={n}: all {cnt} outer-layer configurations {ok}')
