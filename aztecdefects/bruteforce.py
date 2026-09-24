"""Direct enumeration of coverings (validation only; feasible for n <= 5)."""
from functools import lru_cache
from .lattice import aztec, sorted_cells, partner


def k_distribution(cells):
    """dict k -> number of coverings with k B-W pairs, by memoised enumeration."""
    cells = sorted_cells(cells); idx = {c: k for k, c in enumerate(cells)}; N = len(cells)

    @lru_cache(None)
    def f(mask):
        if mask == (1 << N) - 1:
            return {0: 1}
        k = 0
        while mask >> k & 1:
            k += 1
        i, j = cells[k]; res = {}
        opts = [((i + 1, j), 0), ((i, j + 1), 0)]
        p = partner(i, j)
        if p in idx and idx[p] > k:
            opts.append((p, 1))
        for q, w in opts:
            if q in idx and not mask >> idx[q] & 1:
                for e, c in f(mask | 1 << k | 1 << idx[q]).items():
                    res[e + w] = res.get(e + w, 0) + c
        return res
    d = f(0)
    return {k // 2: v for k, v in sorted(d.items())}   # two diagonal tiles per pair


def enumerate_coverings(cells):
    """Yield every covering as a list of tiles (pairs of cells)."""
    cells = sorted_cells(cells); S = set(cells)

    def rec(used, tiles):
        free = [c for c in cells if c not in used]
        if not free:
            yield list(tiles); return
        i, j = free[0]
        for q in [(i + 1, j), (i, j + 1), partner(i, j)]:
            if q in S and q not in used and (q[1], q[0]) > (j, i):
                yield from rec(used | {(i, j), q}, tiles + [((i, j), q)])
    yield from rec(frozenset(), [])
