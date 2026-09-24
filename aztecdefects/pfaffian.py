"""Pfaffian formulation of coverings by dominoes and diagonal (B/W) tiles.

The covering graph (vertices = cells, edges = dominoes and diagonal tiles) is
planar but not bipartite.  `orientation` builds a Pfaffian (Kasteleyn)
orientation: horizontal edges left->right, vertical edges alternating by
column, and each diagonal oriented so that its triangular faces are
clockwise-odd (paper Sec. III and App. A).  Then Pf K(u) = Z(u) with weight u
on every B tile (and 1 on W tiles; since N_B = N_W the exponent of u counts the
number k of B-W pairs).
"""
import math
import numpy as np
from .lattice import partner, color, sorted_cells
from .modular import pfaff_mod, inv_mod, interpolate_mod, crt, primes_below


def orientation(cells):
    """Return (R, idx, O, tiles).

    R     : cells sorted row by row
    idx   : cell -> index
    O     : dict of oriented edges (a, b) -> 1 meaning a -> b
    tiles : list of admissible diagonal tiles (c, p, 'B'|'W')
    """
    R = sorted_cells(cells); idx = {c: k for k, c in enumerate(R)}; S = set(R); O = {}
    for (i, j) in R:
        if (i + 1, j) in S:
            O[((i, j), (i + 1, j))] = 1
        if (i, j + 1) in S:
            if i % 2 == 0:
                O[((i, j), (i, j + 1))] = 1
            else:
                O[((i, j + 1), (i, j))] = 1
    cw = lambda a, b: 1 if (a, b) in O else 0
    area = lambda x, y, z: (y[0] - x[0]) * (z[1] - x[1]) - (y[1] - x[1]) * (z[0] - x[0])
    tiles = []
    for c in R:
        p = partner(*c)
        if p in S and c < p:
            qs = [q for q in [(c[0], p[1]), (p[0], c[1])] if q in S]
            if qs:
                tri = [c, qs[0], p]
                if area(*tri) > 0:
                    tri = tri[::-1]                      # clockwise
                cyc = [(tri[k], tri[(k + 1) % 3]) for k in range(3)]
                dg = [e for e in cyc if set(e) == {c, p}][0]
                cnt = sum(cw(*e) for e in cyc if set(e) != {c, p})
                O[dg if cnt % 2 == 0 else (dg[1], dg[0])] = 1
            else:
                O[(c, p)] = 1
            tiles.append((c, p, 'B' if color(c) == 0 else 'W'))
    return R, idx, O, tiles


def _edge_arrays(idx, O, tiles):
    tt = {frozenset((c, q)): t for c, q, t in tiles}
    I = np.array([idx[a] for a, b in O]); J = np.array([idx[b] for a, b in O])
    isB = np.array([tt.get(frozenset((a, b))) == 'B' for a, b in O])
    return I, J, isB


def K_mod(R, idx, O, tiles, u, p):
    N = len(R); I, J, isB = _edge_arrays(idx, O, tiles)
    K = np.zeros((N, N), dtype=np.int64); w = np.where(isB, u % p, 1)
    K[I, J] = w; K[J, I] = (-w) % p
    return K


def K_float(R, idx, O, tiles, u):
    N = len(R); I, J, isB = _edge_arrays(idx, O, tiles)
    K = np.zeros((N, N)); w = np.where(isB, u, 1.0)
    K[I, J] = w; K[J, I] = -w
    return K


# ---------------------------------------------------------------- exact polynomials
def generating_polynomial(cells, verbose=False):
    """Exact coefficients c_k of Z(u) = sum_k c_k u^k for an arbitrary simply connected region.

    Uses det K(u) = Z(u)^2: det K(u) is evaluated modulo several primes at
    u = 0..2D with FLINT, interpolated, reconstructed by CRT, and the integer
    polynomial square root is taken.  All intermediate identities are checked.
    """
    import flint
    R, idx, O, tiles = orientation(cells); N = len(R)
    if N % 2:
        return [0]
    D = sum(1 for t in tiles if t[2] == 'B'); deg = 2 * D; pts = list(range(deg + 1))
    bits = 2 * int(1.4427 * 0.37 * N) + 80
    ps = primes_below((1 << 31) - 1, bits // 30 + 2)
    I, J, isB = _edge_arrays(idx, O, tiles); coefs = []
    for p in ps:
        vals = []
        for u in pts:
            A = np.zeros((N, N), dtype=np.int64); w = np.where(isB, u % p, 1)
            A[I, J] = w; A[J, I] = (-w) % p
            vals.append(int(flint.nmod_mat(N, N, A.ravel().tolist(), p).det()))
        coefs.append(interpolate_mod(vals, p)[:, 0])
        if verbose:
            print('  prime', p, 'done', flush=True)
    d = [int(x) for x in crt(coefs, ps)]
    a0 = math.isqrt(d[0]); assert a0 * a0 == d[0]
    a = [a0]
    for m in range(1, D + 1):
        s = d[m] - sum(a[i] * a[m - i] for i in range(1, m))
        assert s % (2 * a0) == 0; a.append(s // (2 * a0))
    for m in range(D + 1, deg + 1):
        assert d[m] == sum(a[i] * a[m - i] for i in range(max(0, m - D), min(m, D) + 1))
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def local_statistics(cells, nprimes=5):
    """Exact k-resolved single-tile and tile-pair counts.

    Returns dict with
      c   : c_k                       (k = 0..D)
      ce  : c_k(e)   [k, tile]        number of k-coverings containing tile e
      cef : c_k(e,f) [k, pair]        for all unordered tile pairs (I, J)
      tiles, I, J, R
    Method: Pf K(u) and K(u)^{-1} modulo primes at u = 0..D, the local formula
    Z_{e1..em}/Z = prod K_e (-1)^m Pf[K^{-1}]_{a1 b1 .. am bm}, interpolation, CRT.
    """
    R, idx, O, tiles = orientation(cells); T = len(tiles)
    D = sum(1 for t in tiles if t[2] == 'B'); pts = list(range(D + 1))
    bits = int(1.4427 * 0.37 * len(R)) + 40
    ps = primes_below(2147483647, max(nprimes, bits // 30 + 2))
    ea = np.array([idx[c] for c, q, t in tiles]); eb = np.array([idx[q] for c, q, t in tiles])
    I, J = np.triu_indices(T, 1); res = []
    for p in ps:
        Zs, Ze, Zef = [], [], []
        for u in pts:
            K = K_mod(R, idx, O, tiles, u, p); pf = pfaff_mod(K, p); Ki = inv_mod(K, p)
            ke = K[ea, eb]
            Ze.append((pf * ((ke * Ki[eb, ea]) % p)) % p)
            v = np.stack([ea[I], eb[I], ea[J], eb[J]], 1)
            Mx = lambda x, y: Ki[v[:, x], v[:, y]]
            pf4 = (Mx(0, 1) * Mx(2, 3) % p - Mx(0, 2) * Mx(1, 3) % p + Mx(0, 3) * Mx(1, 2) % p) % p
            Zef.append((pf * ((ke[I] * ke[J]) % p) % p * pf4) % p); Zs.append(pf)
        res.append((interpolate_mod(Zs, p)[:, 0], interpolate_mod(Ze, p), interpolate_mod(Zef, p)))
    c = crt([r[0] for r in res], ps); ce = crt([r[1] for r in res], ps); cef = crt([r[2] for r in res], ps)
    return dict(R=R, tiles=tiles, c=c, ce=ce, cef=cef, I=I, J=J)


# ---------------------------------------------------------------- float evaluations
def lnZ(cells, u):
    """ln Z(u) in double precision (reliable when the system is not deep in the frozen regime)."""
    R, idx, O, tiles = orientation(cells); K = K_float(R, idx, O, tiles, u)
    s, l = np.linalg.slogdet(K); return l / 2


def tile_probabilities(cells, u, gauge_seed=None):
    """Probability of every diagonal tile at fugacity u (double precision).

    With gauge_seed, a random diagonal gauge is applied (probabilities are gauge
    invariant, so comparing two seeds estimates the numerical error).
    """
    R, idx, O, tiles = orientation(cells); K = K_float(R, idx, O, tiles, u)
    if gauge_seed is not None:
        g = np.exp(np.random.default_rng(gauge_seed).normal(0, .4, len(R))); K = g[:, None] * K * g[None, :]
    Ki = np.linalg.inv(K)
    p = np.array([abs(K[idx[a], idx[b]] * Ki[idx[b], idx[a]]) for a, b, t in tiles])
    return tiles, p
