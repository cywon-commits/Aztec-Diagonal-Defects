"""Schematic figures: model (Fig. 1), n=3 examples (Fig. 2), renewal (Fig. 7), proof (Fig. 8)."""
import numpy as np
from _common import *
from matplotlib.patches import Rectangle, Circle
from aztecdefects.lattice import aztec, partner
from aztecdefects.bruteforce import enumerate_coverings

def diag(ax, c1, c2, col, lw=5):
    ax.plot([c1[0] + .5, c2[0] + .5], [c1[1] + .5, c2[1] + .5], color=col, lw=lw, solid_capstyle='round')

# Fig. 1
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6)); a = axs[0]
for bx in range(3):
    for by in range(3):
        x, y = 2 * bx, 2 * by
        rs = [(x, y, 2, 1), (x, y + 1, 2, 1)] if (bx + by) % 2 == 0 else [(x, y, 1, 2), (x + 1, y, 1, 2)]
        for r in rs:
            a.add_patch(Rectangle(r[:2], r[2], r[3], fc='white', ec='k', lw=.8))
fr = lambda rx, ry, w, h, c: a.add_patch(Rectangle((rx + .08, ry + .08), w - .16, h - .16, fc=c, alpha=.85))
fr(2, 2, 2, 1, BL); fr(4, 2, 1, 2, BL); fr(0, 4, 2, 1, YE); fr(0, 5, 2, 1, YE)
a.set_xlim(-.2, 6.2); a.set_ylim(-.2, 6.6); a.set_aspect('equal'); a.axis('off'); a.set_title('(a) plaquette picture')
b = axs[1]
for i in range(6):
    for j in range(6):
        b.add_patch(Rectangle((i, j), 1, 1, fc='#dcdcdc' if (i + j) % 2 == 0 else 'white', ec='#aaa', lw=.4))
b.add_patch(Rectangle((0.12, 4.12), 1.76, .76, fc='none', ec=BL, lw=2)); diag(b, (2, 2), (3, 1), BCOL); diag(b, (4, 3), (3, 2), WCOL)
b.set_xlim(-.2, 6.2); b.set_ylim(-.2, 6.6); b.set_aspect('equal'); b.axis('off'); b.set_title('(b) chessboard picture')
c = axs[2]
for i in range(-2, 6):
    for j in range(-2, 6):
        c.add_patch(Rectangle((i, j), 1, 1, fc='#dcdcdc' if (i + j) % 2 == 0 else 'white', ec='#bbb', lw=.3))
        p = partner(i, j)
        if (i, j) < p and -2 <= p[0] < 6 and -2 <= p[1] < 6:
            diag(c, (i, j), p, BCOL if (i + j) % 2 == 0 else WCOL, 3)
c.set_xlim(-2.2, 6.2); c.set_ylim(-2.2, 6.2); c.set_aspect('equal'); c.axis('off'); c.set_title('(c) allowed positions')
fig.tight_layout(); save(fig, 'fig1_model')

# Fig. 2: one k=1 covering of A_3 in each placement
fig, axs = plt.subplots(1, 2, figsize=(7, 3.7))
for ax, (ox, name) in zip(axs, [(0, 'A'), (1, 'B')]):
    R = aztec(3, ox, 0); count = 0; example = None
    for cov in enumerate_coverings(R):
        nd = sum(1 for t in cov if partner(*t[0]) == t[1])
        if nd == 2:
            count += 1; example = example or cov
    for (i, j) in R:
        ax.add_patch(Rectangle((i, j), 1, 1, fc='#d9d9d9' if (i + j) % 2 == 0 else 'white', ec='#bbb', lw=.5))
    for p_, q_ in example:
        if partner(*p_) == q_:
            diag(ax, p_, q_, BCOL if sum(p_) % 2 == 0 else WCOL, 8)
        else:
            x0, y0 = min(p_[0], q_[0]), min(p_[1], q_[1])
            ax.add_patch(Rectangle((x0 + .08, y0 + .08), abs(p_[0] - q_[0]) + .84, abs(p_[1] - q_[1]) + .84, fc='none', ec='#333', lw=1.6))
    ax.set_aspect('equal'); ax.axis('off'); ax.autoscale(); ax.set_title(f'placement {name}: {count} coverings with k=1')
fig.tight_layout(); save(fig, 'fig2_n3')

# Fig. 7: renewal gadget and face classes
fig, axs = plt.subplots(1, 2, figsize=(5.5, 2.6)); a = axs[0]
P = {'a': (0, 1), 'b': (1, 1), 'c': (1, 0), 'd': (0, 0)}
for u_, v_ in [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'a')]:
    a.plot([P[u_][0], P[v_][0]], [P[u_][1], P[v_][1]], 'k')
a.plot([0, 1], [1, 0], color=BL, lw=2)
S = 1.6; I_ = {'a': (.45, 1.15), 'b': (1.15, 1.15), 'c': (1.15, .45), 'd': (.45, .45)}; Q = {'a': (0, S), 'b': (S, S), 'c': (S, 0), 'd': (0, 0)}
for k in 'abcd':
    a.plot([Q[k][0] + 2.3, I_[k][0] + 2.3], [Q[k][1] - .3, I_[k][1] - .3], 'k', lw=.8)
for u_, v_ in [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'a')]:
    a.plot([I_[u_][0] + 2.3, I_[v_][0] + 2.3], [I_[u_][1] - .3, I_[v_][1] - .3], 'k')
a.plot([I_['b'][0] + 2.3, I_['d'][0] + 2.3], [I_['b'][1] - .3, I_['d'][1] - .3], color=BL, lw=2)
a.set_aspect('equal'); a.axis('off'); a.set_title('(a) face with chord $\\to$ gadget', fontsize=8)
c = axs[1]; R = aztec(4)
for (i, j) in R:
    c.add_patch(Rectangle((i, j), 1, 1, fc='#f0f0f0' if (i + j) % 2 == 0 else 'white', ec='#ccc', lw=.4))
    p = partner(i, j)
    if p in R and (i, j) < p:
        diag(c, (i, j), p, BCOL if (i + j) % 2 == 0 else WCOL, 1.1)
c.set_aspect('equal'); c.axis('off'); c.autoscale(); c.set_title('(b) chords in $\\mathcal{A}_4$', fontsize=8)
fig.tight_layout(); save(fig, 'fig7_renewal')

# Fig. 8: proof ingredients (n = 5)
n = 5; R = aztec(n); NE = [(i, n - 1 - i) for i in range(n)]; SE = [(n - 1 - j, -1 - j) for j in range(n)]
fig, axs = plt.subplots(1, 2, figsize=(5.5, 2.8))
for ax in axs:
    for (i, j) in R:
        ax.add_patch(Rectangle((i, j), 1, 1, fc='#ececec' if (i + j) % 2 == 0 else 'white', ec='#bbb', lw=.4))
    ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-n - .5, n + .5); ax.set_ylim(-n - .8, n + .5)
for cc, col in [(NE[0], BL), (NE[2], BL), (SE[1], OR), (SE[n - 1], OR)]:
    axs[0].add_patch(Rectangle(cc, 1, 1, fc=col, ec='k'))
axs[0].set_title('(a) $b_0, b_i, w_j, w_{n-1}$', fontsize=8)
crust = [((-k, n - k), (-k, n - k - 1)) for k in range(1, n)] + [((-k, -n + k - 1), (-k, -n + k)) for k in range(1, n + 1)]
for p_, q_ in crust:
    axs[1].add_patch(Rectangle((min(p_[0], q_[0]) + .1, min(p_[1], q_[1]) + .1), .8, 1.8, fc='#9fc3ef', ec='#1c5cab'))
for cc in [NE[0], SE[n - 1]]:
    axs[1].add_patch(Rectangle(cc, 1, 1, fc='k'))
axs[1].set_title('(b) forced crust', fontsize=8)
fig.tight_layout(); save(fig, 'fig8_proof')
