"""Fast validation suite (a few minutes).  Run:  python -m pytest -q tests"""
import itertools, math, random
import numpy as np
import pytest
from aztecdefects.lattice import aztec, boundary_edges, diagonal_tiles
from aztecdefects import pfaffian, kasteleyn, boundary, bruteforce, bulk, shuffling, thermo


@pytest.mark.parametrize('n', [2, 3, 4, 5])
def test_polynomial_vs_bruteforce(n):
    for ox in (0, 1):
        cells = aztec(n, ox, 0)
        bf = bruteforce.k_distribution(cells)
        poly = pfaffian.generating_polynomial(cells)
        assert {k: v for k, v in enumerate(poly) if v} == bf


def test_known_polynomial_n3():
    import sympy as sp
    u = sp.symbols('u')
    c = pfaffian.generating_polynomial(aztec(3))
    assert sp.expand(sum(ci * u ** k for k, ci in enumerate(c)) - (u + 2) ** 2 * (u ** 2 + 16)) == 0


def test_half_step_recursion():
    for m in (1, 2, 3):
        a = pfaffian.generating_polynomial(aztec(2 * m)); b = pfaffian.generating_polynomial(aztec(2 * m - 1))
        assert a == [4 ** m * x for x in b]


def test_local_statistics_vs_enumeration():
    n = 3; cells = aztec(n); ls = pfaffian.local_statistics(cells)
    tid = {frozenset((c, q)): k for k, (c, q, t) in enumerate(ls['tiles'])}
    ce = {}; cef = {}
    for cov in bruteforce.enumerate_coverings(cells):
        diag = [tid[frozenset(t)] for t in cov if frozenset(t) in tid]
        k = sum(1 for d in diag if ls['tiles'][d][2] == 'B')
        for d in diag:
            ce[(k, d)] = ce.get((k, d), 0) + 1
        for a, b in itertools.combinations(sorted(diag), 2):
            cef[(k, a, b)] = cef.get((k, a, b), 0) + 1
    for k in range(ls['ce'].shape[0]):
        for d in range(len(ls['tiles'])):
            assert ls['ce'][k, d] == ce.get((k, d), 0)
        for m in range(len(ls['I'])):
            assert ls['cef'][k, m] == cef.get((k, ls['I'][m], ls['J'][m]), 0)


@pytest.mark.parametrize('n', range(2, 9))
def test_theorem1(n):
    E = boundary_edges(n); P = 2 ** (n * (n - 1) // 2); R = aztec(n)
    for i in range(n):
        for j in range(n):
            assert kasteleyn.count_tilings(R - {E['NE'][i], E['SE'][j]}) == P * boundary.N_sum(n, i, n - 1 - j)


def test_N_recursion_equals_sum():
    for n in (5, 9, 14):
        T = boundary.N_table(n)
        assert all(T[a][b] == boundary.N_sum(n, a, b) for a in range(n) for b in range(n))


@pytest.mark.parametrize('n', [5, 7])
def test_theorem2_random(n):
    random.seed(0); R = aztec(n); E = boundary_edges(n); c0 = 2 ** (n * (n + 1) // 2)
    T = boundary.boundary_tiles(n); Bs = [t for t in T if t[0] in ('NE', 'SW')]; Ws = [t for t in T if t[0] in ('SE', 'NW')]
    for k in (1, 2, 3):
        for _ in range(10):
            if len(Bs) < k:
                break
            Bsel = random.sample(Bs, k); Wsel = random.sample(Ws, k)
            removed = {E[e][s + d] for (e, s) in Bsel + Wsel for d in (0, 1)}
            assert kasteleyn.count_tilings(R - removed) * 2 ** (2 * k * n) == c0 * boundary.det_formula(n, Bsel, Wsel)


def test_k1_float_matches_exact_total():
    for n, total in [(3, 64), (4, 1024), (5, 241920)]:
        J, _, _ = kasteleyn.k1_weights_float(n)
        assert abs(J.sum() * 2 ** (n * (n + 1) // 2) - total) < 1e-6 * total


def test_k1_modular_matches_float():
    from aztecdefects.lattice import tile_center, depth
    n = 15
    logJ, cB, cW = kasteleyn.k1_weights_exact(n, D=6, verbose=False)
    J, Bt, Wt = kasteleyn.k1_weights_float(n)
    sb = [k for k, t in enumerate(Bt) if depth(tile_center(t), n) <= 6]
    sw = [k for k, t in enumerate(Wt) if depth(tile_center(t), n) <= 6]
    Jf = J[np.ix_(sb, sw)]
    assert np.allclose(np.exp(logJ), Jf, rtol=1e-9, atol=1e-12 * Jf.max())


def test_renewal_gadget():
    assert shuffling.renewal_ratio(3, ((0, 0), (1, 0), (1, -1), (0, -1))) == 2


def test_face_classes():
    assert shuffling.face_classes(4) == (9, 16) and shuffling.face_classes(5) == (25, 16)


def test_bulk_dimer_entropy():
    phi, rho = bulk.thermo(0.0, L=300)
    assert abs(phi - 0.915965594177219 / math.pi) < 1e-4


def test_thermo_thresholds_concave():
    c = pfaffian.generating_polynomial(aztec(7))
    d, T, concave = thermo.thresholds(c)
    assert concave and T[0] < T[1] < T[2]
