"""Exact boundary formulas (paper Sec. VII, IX).

Theorem 1:  M(A_n minus {b_i, w_j}) = 2^{n(n-1)/2} N_n(i, n-1-j),
            N_n(a, b) = sum_k 2^k C(n-1-k, a-k) C(n-1-k, b-k),
for b_i on the NE staircase (from the N corner) and w_j on the SE staircase
(from the E corner).

Theorem 2:  for k B tiles on NE/SW and k W tiles on SE/NW (outermost layer),
            c_k / c_0 = 2^{-2kn} |det[ tau(b, w) nu(b, w) ]|,
with tau = -1 iff b in SW and w in NW.
"""
import math
from math import comb, lgamma
import numpy as np
from .lattice import boundary_edges, partner

EDGE_ORDER = ['NE', 'SE', 'SW', 'NW']


def N_sum(n, a, b):
    """N_n(a, b) from the binomial sum (paper Eq. 8)."""
    if not (0 <= a < n and 0 <= b < n):
        return 0
    return sum(2 ** k * comb(n - 1 - k, a - k) * comb(n - 1 - k, b - k) for k in range(min(a, b) + 1))


def N_table(n):
    """N_n(a, b) for all 0 <= a, b < n via N(a,b)+N(a-1,b)+N(a,b-1)-N(a-1,b-1) = C(n,a)C(n,b)."""
    C = [comb(n, k) for k in range(n + 1)]; N = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            N[a][b] = (C[a] * C[b] - (N[a - 1][b] if a else 0) - (N[a][b - 1] if b else 0)
                       + (N[a - 1][b - 1] if a and b else 0))
    return N


def nu(n, Eb, i, Ew, j, N=None):
    """M(A_n minus {b, w}) / 2^{n(n-1)/2} for a black boundary cell (Eb, i) and white (Ew, j)."""
    k1, k2 = EDGE_ORDER.index(Eb), EDGE_ORDER.index(Ew)
    if (k2 - k1) % 4 == 1:
        a, b = i, j                 # Ew follows Eb clockwise
    elif (k1 - k2) % 4 == 1:
        a, b = j, i
    else:
        raise ValueError('cells must lie on adjacent edges')
    return N[a][n - 1 - b] if N is not None else N_sum(n, a, n - 1 - b)


def tau(Eb, Ew):
    return -1 if (Eb == 'SW' and Ew == 'NW') else 1


def boundary_tiles(n):
    """Diagonal tiles lying on the outermost layer: list of (edge, s) with cells edge[s], edge[s+1]."""
    E = boundary_edges(n); out = []
    for name in EDGE_ORDER:
        L = E[name]
        for s in range(n - 1):
            if partner(*L[s]) == L[s + 1]:
                out.append((name, s))
    return out


def det_formula(n, Btiles, Wtiles, N=None):
    """Integer |det[tau * nu]| of Theorem 2; then c_k/c_0 = value / 2^{2kn}."""
    import flint
    blacks = [(e, s + d) for (e, s) in Btiles for d in (0, 1)]
    whites = [(e, s + d) for (e, s) in Wtiles for d in (0, 1)]
    M = flint.fmpz_mat(len(whites), len(blacks),
                       [tau(eb, ew) * nu(n, eb, ib, ew, iw, N) for (ew, iw) in whites for (eb, ib) in blacks])
    return abs(int(M.det()))


def k1_outer_total(n):
    """Exact sum over the four edge combinations of c_1(B,W)/c_0 for outer-layer pairs,
    plus the moments needed for Fig. 10.  Row-by-row big-integer recursion, O(n^2)."""
    C = [comb(n, k) for k in range(n + 1)]; prev = [0] * n; r0 = None
    S0 = Sa = Sb = Saa = Sbb = Sab = 0
    for a in range(n):
        cur = [0] * n
        for b in range(n):
            cur[b] = C[a] * C[b] - prev[b] - (cur[b - 1] if b else 0) + (prev[b - 1] if b else 0)
        if a >= 1 and (a - 1) % 2 == 0:
            for b in range(2, n, 2):
                v = abs(r0[b] * cur[b - 1] - cur[b] * r0[b - 1])
                S0 += v; Sa += v * a; Sb += v * b; Saa += v * a * a; Sbb += v * b * b; Sab += v * a * b
        if a % 2 == 0:
            r0 = cur
        prev = cur
    ma, mb = Sa / S0, Sb / S0
    va = Saa / S0 - ma * ma; vb = Sbb / S0 - mb * mb; cab = Sab / S0 - ma * mb
    lnS = math.log(4 * S0) - 2 * n * math.log(2)
    return dict(ln_c1_over_c0=lnS, var_over_n=va / n, corr=cab / math.sqrt(va * vb))


def k2_split_statistics(n, chunk=3_000_000):
    """k = 2 with B tiles on NE and SW and W tiles on SE and NW (the dominant sector).

    Returns ln(c_2/c_0) of that sector, Var/n of a B position, n*corr(B_NE, B_SW),
    n*corr(B_NE, W_SE).  Float evaluation of the 4x4 determinants after factoring
    out the binomial scale of each cell.
    """
    N = N_table(n); E = boundary_edges(n)
    lf = lambda k: lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    blacks = [('NE', i) for i in range(n)] + [('SW', i) for i in range(n)]
    whites = [('SE', j) for j in range(n)] + [('NW', j) for j in range(n)]
    bi = {c: k for k, c in enumerate(blacks)}; wi = {c: k for k, c in enumerate(whites)}
    H = np.zeros((2 * n, 2 * n))
    for r, (ew, j) in enumerate(whites):
        for c, (eb, i) in enumerate(blacks):
            H[r, c] = tau(eb, ew) * math.exp(math.log(nu(n, eb, i, ew, j, N)) - lf(i) - lf(n - 1 - j))
    T = boundary_tiles(n)
    NE = [t for t in T if t[0] == 'NE']; SW = [t for t in T if t[0] == 'SW']
    SE = [t for t in T if t[0] == 'SE']; NW = [t for t in T if t[0] == 'NW']
    cols = np.array([[bi[(x[0], x[1])], bi[(x[0], x[1] + 1)], bi[(y[0], y[1])], bi[(y[0], y[1] + 1)]] for x in NE for y in SW])
    rows = np.array([[wi[(x[0], x[1])], wi[(x[0], x[1] + 1)], wi[(y[0], y[1])], wi[(y[0], y[1] + 1)]] for x in SE for y in NW])
    lb = np.array([lf(x[1]) + lf(x[1] + 1) + lf(y[1]) + lf(y[1] + 1) for x in NE for y in SW])
    lw = np.array([lf(n - 1 - x[1]) + lf(n - 2 - x[1]) + lf(n - 1 - y[1]) + lf(n - 2 - y[1]) for x in SE for y in NW])
    xB1 = np.array([x[1] for x in NE for y in SW], float); xB2 = np.array([y[1] for x in NE for y in SW], float)
    xW1 = np.array([x[1] for x in SE for y in NW], float)
    logW = np.empty((len(cols), len(rows))); step = max(1, chunk // len(rows))
    for s in range(0, len(cols), step):
        e = min(len(cols), s + step); sub = H[rows[None, :, :, None], cols[s:e, None, None, :]]
        logW[s:e] = np.log(np.abs(np.linalg.det(sub)) + 1e-300) + lb[s:e, None] + lw[None, :] - 4 * n * math.log(2)
    m = logW.max(); P = np.exp(logW - m); Z = P.sum(); P /= Z
    pb = P.sum(1); pw = P.sum(0)
    mom = lambda x, p: ((p * x).sum(), (p * (x - (p * x).sum()) ** 2).sum())
    mB1, vB1 = mom(xB1, pb); mB2, vB2 = mom(xB2, pb); mW1, vW1 = mom(xW1, pw)
    cBB = (pb * (xB1 - mB1) * (xB2 - mB2)).sum() / math.sqrt(vB1 * vB2)
    cBW = (P * (xB1 - mB1)[:, None] * (xW1 - mW1)[None, :]).sum() / math.sqrt(vB1 * vW1)
    return dict(ln_c2_over_c0=m + math.log(Z), var_over_n=vB1 / n, n_corr_BB=n * cBB, n_corr_BW=n * cBW)
