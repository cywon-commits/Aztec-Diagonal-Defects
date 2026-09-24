"""Generalized urban renewal for a square face carrying a diagonal chord (paper Sec. VIII).

Face (a, b, c, d) with edge weights alpha=ab, beta=bc, gamma=cd, delta=da and a
chord mu between a and c is replaced by legs a-a', b-b', c-c', d-d' (weight 1),
an inner square with weights a'b' = gamma/D, b'c' = delta/D, c'd' = alpha/D,
d'a' = beta/D and a chord b'-d' of weight mu/D, where D = alpha*gamma + beta*delta.
Then Z_old = D * Z_new.
"""
from functools import lru_cache
import sympy as sp
from .lattice import aztec, partner

u = sp.symbols('u')


def covering_graph(n):
    R = aztec(n); E = {}
    for (i, j) in R:
        for q in [(i + 1, j), (i, j + 1)]:
            if q in R:
                E[frozenset(((i, j), q))] = sp.Integer(1)
        p = partner(i, j)
        if p in R:
            E[frozenset(((i, j), p))] = u if (i + j) % 2 == 0 else sp.Integer(1)
    return R, E


def matching_polynomial(V, E):
    V = sorted(V, key=str); adj = {v: [] for v in V}
    for e, w in E.items():
        a, b = tuple(e); adj[a].append((b, w)); adj[b].append((a, w))
    order = {v: k for k, v in enumerate(V)}

    @lru_cache(None)
    def f(rem):
        if not rem:
            return sp.Integer(1)
        v = min(rem, key=lambda x: order[x]); s = 0
        for (q, w) in adj[v]:
            if q in rem:
                s += w * f(rem - {v, q})
        return sp.expand(s)
    return f(frozenset(V))


def renewal_ratio(n, face):
    """Z(original) / Z(after renewing `face`) -- should equal alpha*gamma + beta*delta = 2."""
    R, E = covering_graph(n); Z0 = matching_polynomial(R, E)
    a, b, c, d = face; E2 = dict(E)
    for x, y in [(a, b), (b, c), (c, d), (d, a), (a, c)]:
        E2.pop(frozenset((x, y)))
    P = {v: ('inner', v) for v in face}; V2 = set(R) | set(P.values())
    for v in face:
        E2[frozenset((v, P[v]))] = sp.Integer(1)
    for x, y in [(a, b), (b, c), (c, d), (d, a)]:
        E2[frozenset((P[x], P[y]))] = sp.Rational(1, 2)
    E2[frozenset((P[b], P[d]))] = u / 2
    return sp.simplify(Z0 / matching_polynomial(V2, E2))


def face_classes(n):
    """Complete faces of A_n split by the parity of the lattice corner (x + y)."""
    R = aztec(n); pts = {(i + x, j + y) for (i, j) in R for x in (0, 1) for y in (0, 1)}
    even = odd = 0
    for (x, y) in pts:
        if all(q in R for q in [(x - 1, y - 1), (x, y - 1), (x - 1, y), (x, y)]):
            if (x + y) % 2:
                odd += 1
            else:
                even += 1
    return even, odd
