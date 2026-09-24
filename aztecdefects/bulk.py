"""Periodic (bulk) model on the torus with a 2x2 fundamental domain (paper Sec. V, VI, VIII)."""
import math
import numpy as np
from .lattice import aztec
from .pfaffian import orientation

DOMAIN = [(0, 0), (1, 0), (0, 1), (1, 1)]


def periodic_terms():
    """Fourier data of the Pfaffian-oriented covering graph restricted to one 2x2 domain.

    Returns a list of (a, b, sx, sy, sign, type) meaning K[(a,0),(b,(sx,sy))] = sign * weight,
    type 'D' (domino), 'B' or 'W'.  Extracted from the interior of a large diamond.
    """
    R, idx, O, tiles = orientation(aztec(11))
    tt = {frozenset((c, q)): t for c, q, t in tiles}; di = {v: k for k, v in enumerate(DOMAIN)}; terms = []
    for (a, b) in O:
        for (x, y, s) in [(a, b, 1), (b, a, -1)]:
            if x in di:
                rep = (y[0] % 2, y[1] % 2); sh = (y[0] // 2, y[1] // 2)
                terms.append((di[x], di[rep], sh[0], sh[1], s, tt.get(frozenset((a, b)), 'D')))
    return terms


TERMS = None


def _terms():
    global TERMS
    if TERMS is None:
        TERMS = periodic_terms()
    return TERMS


def momentum_matrices(L):
    """M0 (dominoes) and M1 (diagonal tiles, unit weight) on an L x L antiperiodic momentum grid."""
    ks = (np.arange(L) + 0.5) * 2 * np.pi / L; TH, PH = np.meshgrid(ks, ks, indexing='ij')
    M0 = np.zeros((L, L, 4, 4), complex); M1 = np.zeros((L, L, 4, 4), complex)
    for a, b, sx, sy, s, ty in _terms():
        e = s * np.exp(1j * (sx * TH + sy * PH))
        (M0 if ty == 'D' else M1)[:, :, a, b] += e
    return TH, PH, M0, M1


def spectral_min(t, grid=181):
    """min over the unit torus of |det K(e^{i th}, e^{i ph})| with diagonal weight t (u = t^2)."""
    g = np.linspace(-np.pi, np.pi, grid); best = np.inf
    terms = _terms()
    for a_ in g:
        for b_ in g:
            M = np.zeros((4, 4), complex)
            for a, b, sx, sy, s, ty in terms:
                M[a, b] += s * (1 if ty == 'D' else t) * np.exp(1j * (sx * a_ + sy * b_))
            best = min(best, abs(np.linalg.det(M)))
    return best


def bulk_pair_weight(directions=((1, 0), (1, 1)), steps=(1, 2, 3, 4, 6, 8, 11, 16, 22, 32, 45), L=500):
    """Weight of one B tile at the origin and one W tile at separation r in the infinite
    critical (u -> 0) dimer plane: |Pf of the 4x4 block of K^{-1}|."""
    TH, PH, M0, M1 = momentum_matrices(L); Mi = np.linalg.inv(M0)
    di = {v: k for k, v in enumerate(DOMAIN)}
    rep = lambda c: (di[(c[0] % 2, c[1] % 2)], (c[0] // 2, c[1] // 2))
    G = lambda b, Y, c, Z: (Mi[:, :, b, c] * np.exp(-1j * (TH * (Z[0] - Y[0]) + PH * (Z[1] - Y[1])))).mean()
    pf4 = lambda A: A[0, 1] * A[2, 3] - A[0, 2] * A[1, 3] + A[0, 3] * A[1, 2]
    Bc = [(0, 0), (1, -1)]; out = {}
    for d in directions:
        rs, vs = [], []
        for s in steps:
            m, nn = d[0] * s, d[1] * s; V = [rep(x) for x in Bc + [(2 * m, 2 * nn + 1), (2 * m - 1, 2 * nn)]]
            A = np.array([[G(V[i][0], V[i][1], V[j][0], V[j][1]) if i != j else 0 for j in range(4)] for i in range(4)])
            rs.append(math.hypot(2 * m - 1, 2 * nn + 1)); vs.append(abs(pf4(A)))
        out[d] = (np.array(rs), np.array(vs))
    return out


def thermo(u, L=500, mats=None):
    """Bulk free energy phi(u) = ln Z per cell and pair density rho = u dphi/du per cell."""
    TH, PH, M0, M1 = mats if mats is not None else momentum_matrices(L)
    t = math.sqrt(u); M = M0 + t * M1
    phi = np.log(np.abs(np.linalg.det(M))).mean() / 8
    Mi = np.linalg.inv(M)
    rho = np.real(np.einsum('xyij,xyji->xy', Mi, (t / 2) * M1)).mean() / 8
    return phi, rho


def diagonal_fraction(u, L=400, mats=None):
    """Fraction of tiles that are diagonal in the bulk, f_inf(u) = 4 rho(u)."""
    return 4 * thermo(u, L, mats)[1]
