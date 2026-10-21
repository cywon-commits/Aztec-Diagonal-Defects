"""Schematic figures: model (Fig. 1), n=3 examples (Fig. 2), renewal (Fig. 12, file fig7), proof (Fig. 13, file fig8).

    python scripts/09_schematics.py          (seconds)
"""
from _common import *
from matplotlib.patches import Rectangle, Circle
from aztecdefects.lattice import aztec, partner
from aztecdefects.bruteforce import enumerate_coverings

plt.rcParams.update({'mathtext.fontset': 'cm'})
DCOL, GRAY = '#2a78d6', '#dcdcdc'


def diag(ax, c1, c2, col, lw=5):
    ax.plot([c1[0] + .5, c2[0] + .5], [c1[1] + .5, c2[1] + .5], color=col, lw=lw, solid_capstyle='round')


# ---------------- Fig. 1: model
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6)); a = axs[0]
for bx in range(3):
    for by in range(3):
        x, y = 2 * bx, 2 * by
        rs = [(x, y, 2, 1), (x, y + 1, 2, 1)] if (bx + by) % 2 == 0 else [(x, y, 1, 2), (x + 1, y, 1, 2)]
        for r in rs:
            a.add_patch(Rectangle(r[:2], r[2], r[3], fc='white', ec='k', lw=.8))
fr = lambda rx, ry, w, h, c: a.add_patch(Rectangle((rx + .08, ry + .08), w - .16, h - .16, fc=c, ec='none', alpha=.85))
fr(2, 2, 2, 1, DCOL); fr(4, 2, 1, 2, DCOL); fr(0, 4, 2, 1, '#f0a500'); fr(0, 5, 2, 1, '#f0a500')
a.text(3.5, 1.55, r'L (two $\perp$ plaquettes)', ha='center', fontsize=7, color=DCOL)
a.text(1, 6.25, r'square (two $\parallel$ plaquettes)', ha='center', fontsize=7, color='#b07800')
a.set_xlim(-.2, 6.2); a.set_ylim(-.2, 6.6); a.set_aspect('equal'); a.axis('off'); a.set_title('(a) plaquette picture', fontsize=9)
b = axs[1]
for i in range(6):
    for j in range(6):
        b.add_patch(Rectangle((i, j), 1, 1, fc=GRAY if (i + j) % 2 == 0 else 'white', ec='#aaaaaa', lw=.4))
b.add_patch(Rectangle((0.12, 4.12), 1.76, .76, fc='none', ec=DCOL, lw=2))
diag(b, (2, 2), (3, 1), BCOL); diag(b, (4, 3), (3, 2), WCOL)
b.text(1, 5.2, 'domino', ha='center', color=DCOL, fontsize=7)
b.text(2.1, 0.45, 'B tile', color=BCOL, fontsize=7); b.text(4.1, 3.95, 'W tile', color=WCOL, fontsize=7)
b.set_xlim(-.2, 6.2); b.set_ylim(-.2, 6.6); b.set_aspect('equal'); b.axis('off'); b.set_title('(b) chessboard picture', fontsize=9)
c = axs[2]; L, x0, y0 = 8, -2, -2
for i in range(x0, x0 + L):
    for j in range(y0, y0 + L):
        c.add_patch(Rectangle((i, j), 1, 1, fc=GRAY if (i + j) % 2 == 0 else 'white', ec='#bbbbbb', lw=.3))
for i in range(x0, x0 + L):
    for j in range(y0, y0 + L):
        p = partner(i, j)
        if (i, j) < p and x0 <= p[0] < x0 + L and y0 <= p[1] < y0 + L:
            black = (i + j) % 2 == 0
            diag(c, (i, j), p, BCOL if black else WCOL, 4 if black else 2.5)
c.set_xlim(x0 - .2, x0 + L + .2); c.set_ylim(y0 - .2, y0 + L + .2); c.set_aspect('equal'); c.axis('off'); c.set_title('(c) allowed positions', fontsize=9)
fig.tight_layout(); save(fig, 'fig1_model')

# ---------------- Fig. 2: one k=1 covering of A_3 in each placement (counts by enumeration)
fig, axs = plt.subplots(1, 2, figsize=(7, 3.7))
for ax, (ox, name) in zip(axs, [(0, 'A'), (1, 'B')]):
    R = aztec(3, ox, 0); count = 0; example = None
    for cov in enumerate_coverings(R):
        if sum(1 for t in cov if partner(*t[0]) == t[1]) == 2:
            count += 1; example = example or cov
    for (i, j) in R:
        ax.add_patch(Rectangle((i, j), 1, 1, fc='#d9d9d9' if (i + j) % 2 == 0 else 'white', ec='#bbbbbb', lw=.5))
    for p_, q_ in example:
        if partner(*p_) == q_:
            col = BCOL if sum(p_) % 2 == 0 else '#d62728'
            diag(ax, p_, q_, col, 9); ax.plot([p_[0] + .5, q_[0] + .5], [p_[1] + .5, q_[1] + .5], 'o', color=col, ms=14)
        else:
            xa, ya = min(p_[0], q_[0]), min(p_[1], q_[1])
            ax.add_patch(Rectangle((xa + .08, ya + .08), abs(p_[0] - q_[0]) + .84, abs(p_[1] - q_[1]) + .84, fc='none', ec='#333', lw=1.6))
    xs = [q[0] for q in R]; ys = [q[1] for q in R]
    for x in sorted(set(xs)):
        ax.text(x + .5, min(ys) - .2, str(x), ha='center', va='top', fontsize=7, color='#777')
    for y in sorted(set(ys)):
        ax.text(min(xs) - .2, y + .5, str(y), ha='right', va='center', fontsize=7, color='#777')
    ax.set_xlim(min(xs) - .8, max(xs) + 1.3); ax.set_ylim(min(ys) - .8, max(ys) + 1.3)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(f'Placement {name}: $n=3$, $k=1$  ({count} coverings)', fontsize=11)
    print(f'placement {name}, n=3: {count} coverings with k=1')
fig.text(0.5, 0.02, 'gray = black cell ($i+j$ even), white = white cell | black bar = B tile, red bar = W tile, outlined = domino', ha='center', fontsize=9)
fig.tight_layout(rect=(0, 0.05, 1, 0.95)); save(fig, 'fig2_n3')

# ---------------- Fig. 7: generalized urban renewal and face classes
fig, axs = plt.subplots(1, 3, figsize=(7, 2.5), gridspec_kw={'width_ratios': [1, 1.15, 1.25]})


def vert(ax, p, col, lab='', off=(0, 0)):
    ax.add_patch(Circle(p, .09, fc=col if col != 'w' else 'white', ec='k', lw=1, zorder=5))
    if lab:
        ax.text(p[0] + off[0], p[1] + off[1], lab, ha='center', va='center', fontsize=9)


a = axs[0]; P = {'a': (0, 1), 'b': (1, 1), 'c': (1, 0), 'd': (0, 0)}
for u, v, l, o in [('a', 'b', r'$\alpha$', (0, .12)), ('b', 'c', r'$\beta$', (.12, 0)), ('c', 'd', r'$\gamma$', (0, -.13)), ('d', 'a', r'$\delta$', (-.13, 0))]:
    a.plot([P[u][0], P[v][0]], [P[u][1], P[v][1]], 'k', lw=1.2)
    a.text((P[u][0] + P[v][0]) / 2 + o[0], (P[u][1] + P[v][1]) / 2 + o[1], l, ha='center', va='center')
a.plot([0, 1], [1, 0], color=DCOL, lw=2); a.text(.62, .62, r'$\mu$', color=DCOL)
for k, col, off in [('a', BCOL, (-.17, .12)), ('b', 'w', (.17, .12)), ('c', BCOL, (.17, -.12)), ('d', 'w', (-.17, -.12))]:
    vert(a, P[k], col, '$%s$' % k, off)
a.set_xlim(-.4, 1.4); a.set_ylim(-.4, 1.4); a.set_aspect('equal'); a.axis('off'); a.set_title('(a) face with a chord')
b = axs[1]; S = 1.6
I_ = {'a': (.45, 1.15), 'b': (1.15, 1.15), 'c': (1.15, .45), 'd': (.45, .45)}; Q = {'a': (0, S), 'b': (S, S), 'c': (S, 0), 'd': (0, 0)}
for k in 'abcd':
    b.plot([Q[k][0], I_[k][0]], [Q[k][1], I_[k][1]], 'k', lw=1)
for u, v, l, o in [('a', 'b', r'$\gamma/\Delta$', (0, .13)), ('b', 'c', r'$\delta/\Delta$', (.25, 0)), ('c', 'd', r'$\alpha/\Delta$', (0, -.13)), ('d', 'a', r'$\beta/\Delta$', (-.25, 0))]:
    b.plot([I_[u][0], I_[v][0]], [I_[u][1], I_[v][1]], 'k', lw=1.2)
    b.text((I_[u][0] + I_[v][0]) / 2 + o[0], (I_[u][1] + I_[v][1]) / 2 + o[1], l, ha='center', va='center', fontsize=8)
b.plot([I_['b'][0], I_['d'][0]], [I_['b'][1], I_['d'][1]], color=DCOL, lw=2); b.text(.86, .72, r'$\mu/\Delta$', color=DCOL, fontsize=8)
for k, col in [('a', BCOL), ('b', 'w'), ('c', BCOL), ('d', 'w')]:
    vert(b, Q[k], col); vert(b, I_[k], 'w' if col == BCOL else BCOL)
for k, off in [('a', (-.15, .12)), ('b', (.15, .12)), ('c', (.15, -.12)), ('d', (-.15, -.12))]:
    b.text(Q[k][0] + off[0], Q[k][1] + off[1], '$%s$' % k, ha='center')
b.text(S / 2, -.42, r'$Z=\Delta\,Z^{\prime},\ \ \Delta=\alpha\gamma+\beta\delta$', ha='center', fontsize=8)
b.set_xlim(-.4, S + .4); b.set_ylim(-.55, S + .35); b.set_aspect('equal'); b.axis('off'); b.set_title('(b) generalized urban renewal')
c = axs[2]; n = 4; R = aztec(n)
for (i, j) in R:
    c.add_patch(Rectangle((i, j), 1, 1, fc='#f0f0f0' if (i + j) % 2 == 0 else 'white', ec='#cccccc', lw=.4))
pts = {(i + x, j + y) for (i, j) in R for x in (0, 1) for y in (0, 1)}
for (x, y) in pts:
    if all(q in R for q in [(x - 1, y - 1), (x, y - 1), (x - 1, y), (x, y)]):
        if (x + y) % 2:
            c.plot(x, y, 'o', color=DCOL, ms=3.5)
        else:
            c.plot(x, y, 's', color='#aaaaaa', ms=3)
for cc in R:
    p = partner(*cc)
    if p in R and cc < p:
        c.plot([cc[0] + .5, p[0] + .5], [cc[1] + .5, p[1] + .5], color=BCOL if sum(cc) % 2 == 0 else WCOL, lw=1.1, alpha=.8)
c.set_aspect('equal'); c.axis('off'); c.set_title(r'(c) face classes, $\mathcal{A}_4$', fontsize=9)
c.text(0, -5.6, 'blue: $n^2$ complete faces with a chord\ngray: chord-free faces', ha='center', fontsize=7)
fig.tight_layout(); save(fig, 'fig7_renewal')

# ---------------- Fig. 8: proof ingredients (n = 5)
fig, axs = plt.subplots(1, 3, figsize=(7, 2.6))
n = 5; R = aztec(n); NE = [(i, n - 1 - i) for i in range(n)]; SE = [(n - 1 - j, -1 - j) for j in range(n)]


def base(ax):
    for (i, j) in R:
        ax.add_patch(Rectangle((i, j), 1, 1, fc='#ececec' if (i + j) % 2 == 0 else 'white', ec='#bbbbbb', lw=.4))
    ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-n - .5, n + .5); ax.set_ylim(-n - .8, n + .5)


def mark(ax, c_, col, lab, dx, dy):
    ax.add_patch(Rectangle(c_, 1, 1, fc=col, ec='k', lw=.8)); ax.text(c_[0] + .5 + dx, c_[1] + .5 + dy, lab, fontsize=8, ha='center')


def dom(ax, p, q):
    xa, ya = min(p[0], q[0]), min(p[1], q[1]); ax.add_patch(Rectangle((xa + .1, ya + .1), .8, 1.8, fc='#9fc3ef', ec='#1c5cab', lw=1))


a = axs[0]; base(a); i, j = 2, 1
mark(a, NE[0], DCOL, r'$b_0$', 0, .95); mark(a, NE[i], DCOL, r'$b_i$', .9, .6)
mark(a, SE[j], OR, r'$w_j$', 1.0, -.2); mark(a, SE[n - 1], OR, r'$w_{n-1}$', .2, -1.0)
a.set_title('(a) four boundary cells (Kuo)')
b = axs[1]; base(b)
crust = [((-k, n - k), (-k, n - k - 1)) for k in range(1, n)] + [((-k, -n + k - 1), (-k, -n + k)) for k in range(1, n + 1)]
for p, q in crust:
    dom(b, p, q)
for c_ in [NE[0], SE[n - 1]]:
    b.add_patch(Rectangle(c_, 1, 1, fc='k', ec='k'))
Rm = aztec(n - 1, 1, 0)
assert R - {NE[0], SE[n - 1]} - {x for t in crust for x in t} == Rm   # Lemma 2: the remainder is A_{n-1} + (1,0)
for (x, y) in Rm:
    for dx, dy, seg in [(1, 0, ((x + 1, y), (x + 1, y + 1))), (-1, 0, ((x, y), (x, y + 1))), (0, 1, ((x, y + 1), (x + 1, y + 1))), (0, -1, ((x, y), (x + 1, y)))]:
        if (x + dx, y + dy) not in Rm:
            b.plot([seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]], color=OR, lw=1.4)
b.set_title(r'(b) forced crust $\Rightarrow\mathcal{A}_{n-1}$')
c = axs[2]; base(c)
for p, q in crust[n - 1:]:
    dom(c, p, q)
c.add_patch(Rectangle(SE[n - 1], 1, 1, fc='k'))
Qc = [(0, n - 1)] + [cc for k in range(1, n) for cc in [(-k, n - k), (-k, n - k - 1)]]
for q in Qc:
    c.add_patch(Rectangle(q, 1, 1, fc='#f7c9a8', ec=OR, lw=1))
c.plot([q[0] + .5 for q in Qc], [q[1] + .5 for q in Qc], '-', color='#8a3210', lw=1)
for k in range(1, n):
    c.add_patch(Rectangle((-k + 1, n - k - 1), 1, 1, fc='none', ec=AQ, lw=1.6))
c.text(-2.5, n - .2, r'path $Q$', color='#8a3210', fontsize=8); c.text(1.2, 1.3, r'$u_s$', color=AQ, fontsize=8)
c.set_title(r'(c) path $Q$ and cells $u_s$')
fig.tight_layout(); save(fig, 'fig8_proof')
