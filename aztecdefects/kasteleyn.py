"""Pure-domino (bipartite) Kasteleyn matrices and exact single-pair (k=1) weights.

For k = 1 the weight of a B tile e and a W tile f is
    c_1(e, f) / c_0 = M(A_n minus the four cells) / M(A_n) = |det 2x2 of K^{-1}|,
where K is the bipartite Kasteleyn matrix of the pure domino problem.
"""
import math, sys, time
import numpy as np
from .lattice import aztec, sorted_cells, diagonal_tiles, tile_center, depth
from .modular import primes_below


def bipartite_K(cells):
    """Real Kasteleyn matrix: rows = black cells, columns = white cells (both row-by-row)."""
    Bl = sorted_cells(c for c in cells if sum(c) % 2 == 0); Wh = sorted_cells(c for c in cells if sum(c) % 2)
    bi = {c: k for k, c in enumerate(Bl)}; wi = {c: k for k, c in enumerate(Wh)}
    K = np.zeros((len(Bl), len(Wh)), dtype=np.int64)
    for (i, j), r in bi.items():
        for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            q = (i + di, j + dj)
            if q in wi:
                K[r, wi[q]] = 1 if dj == 0 else (-1 if i % 2 else 1)
    return Bl, Wh, bi, wi, K


def count_tilings(cells):
    """Exact number of domino tilings of a simply connected region (FLINT integer determinant)."""
    import flint
    Bl, Wh, bi, wi, K = bipartite_K(set(cells))
    if len(Bl) != len(Wh):
        return 0
    if not Bl:
        return 1
    N = len(Bl)
    return abs(int(flint.fmpz_mat(N, N, [int(v) for v in K.ravel()]).det()))


def k1_weights_float(n, ox=0):
    """All k=1 pair weights J[B, W] = c_1(B,W)/c_0 in double precision (reliable for n <~ 50)."""
    cells = aztec(n, ox, 0); Bl, Wh, bi, wi, K = bipartite_K(cells)
    T = diagonal_tiles(cells); Bt = [t for t in T if t[2] == 'B']; Wt = [t for t in T if t[2] == 'W']
    G = np.linalg.inv(K.astype(float))           # G[w, b]
    b1 = np.array([bi[t[0]] for t in Bt]); b2 = np.array([bi[t[1]] for t in Bt])
    w1 = np.array([wi[t[0]] for t in Wt]); w2 = np.array([wi[t[1]] for t in Wt])
    J = np.abs(G[np.ix_(w1, b1)] * G[np.ix_(w2, b2)] - G[np.ix_(w1, b2)] * G[np.ix_(w2, b1)]).T
    return J, Bt, Wt


# ------------------------------------------------------- multimodular banded solver
def solve_banded_mod(K, rhs_rows, p):
    """X = K^{-1}[:, rhs_rows] mod p for a banded square K (partial pivoting inside the band)."""
    N = K.shape[0]; A = (K % p).astype(np.int64)
    nz = np.nonzero(K); l = int((nz[0] - nz[1]).max()); u = int((nz[1] - nz[0]).max())
    E = np.zeros((N, len(rhs_rows)), dtype=np.int64); E[rhs_rows, np.arange(len(rhs_rows))] = 1
    for k in range(N):
        r1 = min(N, k + l + 1); c1 = min(N, k + l + u + 1)
        piv = k + np.nonzero(A[k:r1, k])[0][0]
        if piv != k:
            A[[k, piv], k:c1] = A[[piv, k], k:c1]; E[[k, piv]] = E[[piv, k]]
        inv = pow(int(A[k, k]), p - 2, p)
        A[k, k:c1] = (A[k, k:c1] * inv) % p; E[k] = (E[k] * inv) % p
        if r1 > k + 1:
            f = A[k + 1:r1, k]; m = np.nonzero(f)[0]
            if len(m):
                rows = k + 1 + m; fm = f[m][:, None].copy()
                A[rows, k:c1] = (A[rows, k:c1] - fm * A[k, k:c1][None, :]) % p
                E[rows] = (E[rows] - fm * E[k][None, :]) % p
    w = l + u
    for k in range(N - 2, -1, -1):
        c1 = min(N, k + w + 1); seg = A[k, k + 1:c1]; nzs = np.nonzero(seg)[0]
        if len(nzs):
            E[k] = (E[k] - (seg[nzs][None, :] @ E[k + 1 + nzs]) % p) % p
    return E


def k1_weights_exact(n, D=6, verbose=True, ox=0):
    """Exact k=1 pair weights for all B and W tiles within `D` layers of the boundary.

    Returns (logJ, cB, cW) with logJ[B, W] = ln(c_1(B,W)/c_0) (exact up to the final
    float conversion).  c_1(e,f)/c_0 is a dyadic rational with denominator
    <= 2^{2n+O(1)}; with s = 2n + 8 the numbers 2^s c_1/c_0 are integers,
    reconstructed from primes < 2^26 (products stay < 2^63 in the banded solver).
    """
    cells = aztec(n, ox, 0); Bl, Wh, bi, wi, K = bipartite_K(cells)
    T = diagonal_tiles(cells)
    Bt = [t for t in T if t[2] == 'B' and depth(tile_center(t), n, ox) <= D]
    Wt = [t for t in T if t[2] == 'W' and depth(tile_center(t), n, ox) <= D]
    bcells = sorted({bi[c] for t in Bt for c in t[:2]}); bpos = {b: k for k, b in enumerate(bcells)}
    b1 = np.array([bpos[bi[t[0]]] for t in Bt]); b2 = np.array([bpos[bi[t[1]]] for t in Bt])
    w1 = np.array([wi[t[0]] for t in Wt]); w2 = np.array([wi[t[1]] for t in Wt])
    s = 2 * n + 8; nb = int(1.45 * 1.39 * n) + s + 60
    primes = primes_below(1 << 26, nb // 26 + 1); res = []
    for p in primes:
        t0 = time.time(); X = solve_banded_mod(K, np.array(bcells), p)
        A = X[np.ix_(w1, b1)]; B = X[np.ix_(w1, b2)]; C = X[np.ix_(w2, b1)]; Dm = X[np.ix_(w2, b2)]
        res.append(((A * Dm % p - B * C % p) % p).T.copy())
        if verbose:
            print(f'  n={n} prime {p}: {time.time()-t0:.1f}s', file=sys.stderr, flush=True)
    Mt = 1; x = np.zeros(res[0].shape, dtype=object)
    for r, p in zip(res, primes):                     # Garner
        r2 = (r.astype(object) * pow(2, s, p)) % p
        t = ((r2 - (x % p)) * pow(Mt % p, -1, p)) % p
        x = x + Mt * t; Mt *= p
    x = np.where(x > Mt // 2, x - Mt, x)
    bits = max(int(abs(v)).bit_length() for v in x.ravel())
    assert bits < Mt.bit_length() - 30, 'CRT bound too small'
    logJ = np.array([[math.log(abs(int(v))) - s * math.log(2) if v != 0 else -np.inf for v in row] for row in x])
    return logJ, np.array([tile_center(t) for t in Bt]), np.array([tile_center(t) for t in Wt])
