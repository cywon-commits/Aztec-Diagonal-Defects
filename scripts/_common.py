import os, sys, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
DATA = os.path.join(ROOT, 'data'); FIG = os.path.join(ROOT, 'figures')
os.makedirs(DATA, exist_ok=True); os.makedirs(FIG, exist_ok=True)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'serif', 'font.size': 9, 'mathtext.fontset': 'cm'})
BL, OR, AQ, YE, VI = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#4a3aa7'
BCOL, WCOL = '#111111', '#e34948'


def clean(ax):
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + '.pdf')); fig.savefig(os.path.join(FIG, name + '.png'), dpi=130)
    print('wrote figures/' + name + '.pdf')


def load_polys():
    p = os.path.join(DATA, 'polynomials.json')
    if not os.path.exists(p):
        raise SystemExit('run scripts/01_polynomials.py first')
    d = json.load(open(p)); return {int(k): [int(x) for x in v] for k, v in d['A'].items()}, \
        {int(k): [int(x) for x in v] for k, v in d.get('B', {}).items()}
