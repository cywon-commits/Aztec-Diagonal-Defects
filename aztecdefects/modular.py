"""Modular linear algebra helpers (numpy int64, primes < 2^31)."""
import numpy as np
from functools import reduce


def primes_below(bound, k):
    """k largest primes below `bound`."""
    import sympy
    ps = []; p = bound
    while len(ps) < k:
        p = sympy.prevprime(p); ps.append(p)
    return ps


def pfaff_mod(A, p):
    """Pfaffian of a skew-symmetric int64 matrix modulo p (p < 2^31)."""
    A = A.copy() % p; N = len(A); res = 1
    for k in range(0, N, 2):
        nz = np.nonzero(A[k, k + 1:])[0]
        if len(nz) == 0:
            return 0
        j = k + 1 + nz[0]
        if j != k + 1:
            A[[k + 1, j]] = A[[j, k + 1]]; A[:, [k + 1, j]] = A[:, [j, k + 1]]; res = (-res) % p
        a = int(A[k, k + 1]); res = res * a % p
        if k + 2 < N:
            ia = pow(a, p - 2, p); r0 = A[k, k + 2:]; r1 = A[k + 1, k + 2:]
            r0i = (r0 * ia) % p
            A[k + 2:, k + 2:] = (A[k + 2:, k + 2:] + (np.outer(r1, r0i) % p) - (np.outer(r0i, r1) % p)) % p
    return res


def inv_mod(A, p):
    """Matrix inverse modulo p by Gauss-Jordan elimination."""
    N = len(A); M = np.concatenate([A % p, np.eye(N, dtype=np.int64)], 1)
    for c in range(N):
        r = c + np.nonzero(M[c:, c])[0][0]
        if r != c:
            M[[c, r]] = M[[r, c]]
        iv = pow(int(M[c, c]), p - 2, p); M[c] = (M[c] * iv) % p
        f = M[:, c].copy(); f[c] = 0
        M = (M - (np.outer(f, M[c]) % p)) % p
    return M[:, N:]


def interpolate_mod(values, p):
    """Coefficients (mod p) of the polynomial taking `values[u]` at u = 0..len-1."""
    pts = list(range(len(values)))
    V = np.array([[pow(u, e, p) for e in pts] for u in pts], dtype=np.int64)
    Vi = inv_mod(V, p)
    Y = np.array(values, dtype=np.int64).reshape(len(pts), -1) % p
    out = np.zeros_like(Y)
    for a in range(len(pts)):
        acc = np.zeros(Y.shape[1], dtype=np.int64)
        for b in range(len(pts)):
            acc = (acc + Vi[a, b] * Y[b] % p) % p
        out[a] = acc
    return out


def crt(residues, primes, signed=True):
    """Chinese remaindering of equally-shaped integer arrays -> object array of ints."""
    M = reduce(lambda a, b: a * b, primes)
    arrs = [np.asarray(r) for r in residues]
    out = np.zeros(arrs[0].shape, dtype=object); flat = [a.ravel() for a in arrs]; o = out.ravel()
    coef = [(M // p) * pow(M // p, -1, p) for p in primes]
    for k in range(len(o)):
        x = sum(int(f[k]) * c for f, c in zip(flat, coef)) % M
        o[k] = x - M if (signed and x > M // 2) else x
    return out
