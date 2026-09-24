"""Lattice conventions.

Cells are unit squares labelled by their lower-left corner (i, j).  A cell is
*black* if i + j is even and *white* otherwise.

Diagonal tiles (paper Eq. 1-2):
    B : (2m, 2n)   -- (2m+1, 2n-1)      covers two black cells, orientation '\\'
    W : (2m, 2n+1) -- (2m-1, 2n)        covers two white cells, orientation '/'
Every cell has exactly one diagonal partner.
"""

def color(c):
    """0 = black (i+j even), 1 = white."""
    return (c[0] + c[1]) % 2


def partner(i, j):
    """Diagonal partner of cell (i, j)."""
    if (i + j) % 2 == 0:                      # black -> B tile
        return (i + 1, j - 1) if i % 2 == 0 else (i - 1, j + 1)
    return (i - 1, j - 1) if i % 2 == 0 else (i + 1, j + 1)   # white -> W tile


def aztec(n, ox=0, oy=0):
    """Cells of the Aztec diamond of order n.

    Placement A (the one used in the paper) is ox = oy = 0; placement B is
    ox = 1, oy = 0.  Cell (i, j) belongs to placement A iff
    |i + 1/2| + |j + 1/2| <= n.
    """
    return {(i + ox, j + oy) for j in range(-n - 1, n + 1) for i in range(-n - 1, n + 1)
            if abs(i + .5) + abs(j + .5) <= n}


def sorted_cells(cells):
    return sorted(cells, key=lambda c: (c[1], c[0]))


def diagonal_tiles(cells):
    """All admissible diagonal tiles inside a cell set: list of (c, p, 'B'|'W')."""
    S = set(cells); T = []
    for c in sorted_cells(S):
        p = partner(*c)
        if p in S and c < p:
            T.append((c, p, 'B' if color(c) == 0 else 'W'))
    return T


def boundary_edges(n, ox=0):
    """The four boundary staircases of A_n (shifted by ox), each indexed clockwise from its corner.

    NE from the N corner, SE from the E corner, SW from the S corner, NW from the
    W corner.  Placement A, odd n / placement B (ox=1), even n: NE, SW black; SE, NW white.
    """
    E = {'NE': [(i, n - 1 - i) for i in range(n)],
         'SE': [(n - 1 - j, -1 - j) for j in range(n)],
         'SW': [(-1 - i, -n + i) for i in range(n)],
         'NW': [(-n + i, i) for i in range(n)]}
    return {k: [(x + ox, y) for x, y in v] for k, v in E.items()}


def tile_center(t):
    a, b = t[0], t[1]
    return ((a[0] + b[0]) / 2 + .5, (a[1] + b[1]) / 2 + .5)


def depth(center, n, ox=0):
    """Distance (in layers) of a tile centre from the boundary; 0 = outermost layer."""
    return n - (abs(center[0] - ox) + abs(center[1]))
