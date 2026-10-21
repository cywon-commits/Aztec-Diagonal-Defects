"""Auxiliary checks quoted in the text (Secs. III, VIII, X, XII, XIII and Appendix B).

    python scripts/12_audit_checks.py        (~2 minutes)

1. Fraction of the k = 1 weight carried by tiles within D = 6 layers of the boundary
   (full double-precision inversion, n <= 51).
2. Ratio r = J2/J1 of the vertex couplings for two limiting island models:
   point dipoles at the island centres (r = 3/sqrt 2) and dumbbells whose charges
   sit at a common distance from the vertex (r = sqrt 2).
3. Half-step identities Z_{2m} = 4^m Z_{2m-1} (placement A) and
   Z^B_{2m+1} = 2^{2m+1} Z^B_{2m} for all cached polynomials.
"""
import math
import numpy as np
from _common import *
from aztecdefects.kasteleyn import k1_weights_float
from aztecdefects.lattice import depth, tile_center

for n in (21, 41, 51):
    J, Bt, Wt = k1_weights_float(n)
    dB = np.array([depth(tile_center(t), n, 0) for t in Bt]); dW = np.array([depth(tile_center(t), n, 0) for t in Wt])
    print(f'n={n}: weight within 6 layers {J[np.ix_(dB <= 6, dW <= 6)].sum() / J.sum():.6f}   outermost layer (B) {J[dB == 0].sum() / J.sum():.4f}')


def pair_energy(r1, m1, r2, m2):
    d = np.subtract(r2, r1); R = np.linalg.norm(d); e = d / R
    return (np.dot(m1, m2) - 3 * np.dot(m1, e) * np.dot(m2, e)) / R ** 3


# vertex at the origin; sigma = +1 means pointing into the vertex
J1 = pair_energy((-.5, 0), (1, 0), (.5, 0), (-1, 0))     # collinear pair, both pointing in
J2 = pair_energy((-.5, 0), (1, 0), (0, .5), (0, -1))     # perpendicular pair, both pointing in
print(f'point dipoles: J1={J1:.4f}  J2={J2:.4f}  r={J2 / J1:.4f}  (3/sqrt2={3 / math.sqrt(2):.4f})')
g = 0.1                                                   # dumbbell charges at distance g from the vertex
print(f'dumbbells:     r = (1/(sqrt2 g))/(1/(2 g)) = {(1 / (math.sqrt(2) * g)) / (1 / (2 * g)):.4f}  (sqrt2={math.sqrt(2):.4f})')

A, B = load_polys()
okA = [m for m in range(1, 20) if 2 * m in A and 2 * m - 1 in A and all(x == 4 ** m * y for x, y in zip(A[2 * m], A[2 * m - 1])) and len(A[2 * m]) == len(A[2 * m - 1])]
okB = [m for m in range(1, 20) if 2 * m + 1 in B and 2 * m in B and all(x == 2 ** (2 * m + 1) * y for x, y in zip(B[2 * m + 1], B[2 * m])) and len(B[2 * m + 1]) == len(B[2 * m])]
print('half-step A verified for m =', okA)
print('half-step B verified for m =', okB)

# 4. Theorem 1 (via Corollary 1 / Theorem 2 with k = 1) against the multimodular data:
#    every outer-layer B-W pair of every cached exact k = 1 computation
from aztecdefects import boundary
from aztecdefects.lattice import boundary_edges
import os
for f in sorted(os.listdir(DATA)):
    if not (f.startswith('k1_exact_n') or f.startswith('k1_exact_B_n')):
        continue
    ox = 1 if '_B_' in f else 0; n = int(f.split('_n')[1][:-4])
    d = np.load(os.path.join(DATA, f)); L, cB, cW = d['L'], d['cB'], d['cW']
    E = boundary_edges(n, ox); N = boundary.N_table(n)
    def centre(e, s): return tuple(np.mean([E[e][s], E[e][s + 1]], axis=0) + .5)
    idxB = {tuple(np.round(c, 3)): k for k, c in enumerate(cB)}; idxW = {tuple(np.round(c, 3)): k for k, c in enumerate(cW)}
    tiles = boundary.boundary_tiles(n, ox); dev = 0.0; cnt = 0
    for (eb, sb) in [t for t in tiles if t[0] in ('NE', 'SW')]:
        for (ew, sw) in [t for t in tiles if t[0] in ('SE', 'NW')]:
            M = [[boundary.tau(eb, ew) * boundary.nu(n, eb, sb + p, ew, sw + q, N) for p in (0, 1)] for q in (0, 1)]
            val = abs(M[0][0] * M[1][1] - M[0][1] * M[1][0])
            lnT = math.log(val) - 2 * n * math.log(2)
            kb = idxB[tuple(np.round(centre(eb, sb), 3))]; kw = idxW[tuple(np.round(centre(ew, sw), 3))]
            dev = max(dev, abs(lnT - L[kb, kw])); cnt += 1
    print(f'{f}: {cnt} outer-layer pairs, max |ln(theorem) - ln(multimodular)| = {dev:.1e}')

# 5. Double precision versus the exact multimodular weights (Appendix B, Sec. III)
from aztecdefects.lattice import depth as _depth
for n in (21, 41, 61):
    d = np.load(os.path.join(DATA, f'k1_exact_n{n}.npz')); L = d['L']
    J, Bt, Wt = k1_weights_float(n)
    sb = [k for k, t in enumerate(Bt) if _depth(tile_center(t), n) <= 6]; sw = [k for k, t in enumerate(Wt) if _depth(tile_center(t), n) <= 6]
    Jf = J[np.ix_(sb, sw)]; Je = np.exp(L)
    m = Je > 1e-6 * Je.max()
    print(f'n={n}: max relative deviation float/exact over significant pairs = {np.abs(Jf[m] / Je[m] - 1).max():.1e}')

# 6. Two pairs (Sec. XII): role of the sign tau, outer-layer share of c_2 at n = 9,
#    and the probability that both B tiles share an edge (all sectors, float 4x4 determinants)
import itertools
from math import lgamma
for n in (7, 9):
    T = boundary.boundary_tiles(n); Bs = [t for t in T if t[0] in ('NE', 'SW')]; Ws = [t for t in T if t[0] in ('SE', 'NW')]
    Nn = boundary.N_table(n); fail = 0; tot = 0; cnt = 0
    for Bp in itertools.combinations(Bs, 2):
        for Wp in itertools.combinations(Ws, 2):
            blacks = [(e, s + d) for (e, s) in Bp for d in (0, 1)]; whites = [(e, s + d) for (e, s) in Wp for d in (0, 1)]
            v = boundary.det_formula(n, Bp, Wp, Nn)
            import flint
            M0 = flint.fmpz_mat(4, 4, [boundary.nu(n, eb, ib, ew, iw, Nn) for (ew, iw) in whites for (eb, ib) in blacks])
            fail += abs(int(M0.det())) != v; tot += v; cnt += 1
    c2 = A[n][2]; c0 = A[n][0]
    print(f'n={n}: {cnt} outer-layer k=2 configurations, {fail} fail without tau; outer-layer share of c_2 = {tot * c0 / 2 ** (4 * n) / c2:.3f}')


k2_all_sectors = boundary.k2_all_sectors


for n in (41, 81, 121):
    lnc2, psame = k2_all_sectors(n)
    print(f'n={n}: all-sector n^6 c2/(c0 16^n) = {math.exp(lnc2 + 6 * math.log(n) - 4 * n * math.log(2)):.5f};  n^2 P(two B on one edge) = {n * n * psame:.3f}')

# 7. Placement B: Theorem 2 against direct counts for k <= 2 (n = 4, 6, 8) and threshold temperatures
from aztecdefects.kasteleyn import count_tilings
from aztecdefects.lattice import aztec
from aztecdefects import thermo
for n in (4, 6, 8):
    R = aztec(n, 1, 0); E = boundary_edges(n, 1); c0 = 2 ** (n * (n + 1) // 2)
    T = boundary.boundary_tiles(n, 1); Bs = [t for t in T if t[0] in ('NE', 'SW')]; Ws = [t for t in T if t[0] in ('SE', 'NW')]
    ok = True; cnt = 0
    for k in (1, 2):
        for Bp in itertools.combinations(Bs, k):
            for Wp in itertools.combinations(Ws, k):
                removed = {E[e][s + d] for (e, s) in Bp + Wp for d in (0, 1)}
                ok &= count_tilings(R - removed) * 2 ** (2 * k * n) == c0 * boundary.det_formula(n, Bp, Wp); cnt += 1
    print(f'placement B, n={n}: Theorem 2 for all {cnt} outer-layer configurations with k<=2: {ok}')
print('T_1 (units eps_d):', {f'{"B" if n % 2 == 0 else "A"}{n}': round(float(thermo.thresholds((B if n % 2 == 0 else A)[n])[1][0]), 3) for n in (10, 11, 12)})
