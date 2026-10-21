import os as _os, sys as _sys
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.path.dirname(_os.path.dirname(_HERE))
_sys.path.insert(0,_HERE); _sys.path.insert(0,_ROOT)
import sys, math

from aztec_sample import AztecSample
from opentiling import tiling_open
from ortools.sat.python import cp_model
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Circle, FancyArrow
plt.rcParams.update({'font.family':'serif','font.size':9,'mathtext.fontset':'cm'})
n=3; s=AztecSample(n,0,2,3)
m,x,E,Q=s.model(); m.Minimize(E); sv=cp_model.CpSolver(); sv.parameters.num_workers=2; sv.Solve(m); E0=int(sv.ObjectiveValue())
# pick a k=1 state with the B-W pair at two adjacent tangency points (most probable pair)
m,x,E,Q=s.model(); m.Add(sum(Q.values())==2); m.Add(E==E0+4)
target=[v for v in Q if v in [(13,13),(15,13)]]
for v in target: m.Add(Q[v]==1)
sv=cp_model.CpSolver(); sv.parameters.num_workers=2; st=sv.Solve(m); print(sv.StatusName(st))
S=[1 if sv.Value(v) else -1 for v in x]; T=tiling_open(s,S)
TL=sorted(T,key=lambda p:sorted(p)); CI={pr:k for k,pr in enumerate(TL)}
import matplotlib.cm as cm
cols=['#aec7e8','#ffbb78','#98df8a','#c5b0d5','#c49c94','#f7b6d2','#dbdb8d','#9edae5','#ffe28a','#b5cf6b','#e7ba52','#9c9ede']
fig,ax=plt.subplots(1,2,figsize=(7.0,3.5),gridspec_kw={'width_ratios':[1.15,1]})
a=ax[0]
pl=s.pl; import itertools
for pr in TL:
    k=CI[pr]
    cs=sorted(pr); diag=abs(cs[0][0]-cs[1][0])==1 and abs(cs[0][1]-cs[1][1])==1
    for c in cs:
        p=pl[c]; (x0,y0),(x1,y1)=p[1],p[2]
        fc=('#555555' if (c[0]+c[1])%2==0 else '#e34948') if diag else cols[k%len(cols)]
        a.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fc=fc,ec='white',lw=1.2,alpha=0.6 if diag else 0.9))
for (p,q),Si in zip(s.edges,S):
    # arrow along island: +1 means p->q ; S index sign convention: spin +1 points a->b (out of a)
    (x0,y0),(x1,y1)=p,q
    if Si<0: (x0,y0),(x1,y1)=q,p
    a.annotate('',xy=(x0+0.8*(x1-x0),y0+0.8*(y1-y0)),xytext=(x0+0.2*(x1-x0),y0+0.2*(y1-y0)),arrowprops=dict(arrowstyle='-|>',lw=0.9,color='#222',mutation_scale=7))
for v,(c1,c2) in s.sites.items():
    arms=s.VA[v]; q=sum(sg*S[i] for d,i,sg in arms)
    if q: a.add_patch(Circle(v,0.28,fc='#d62728' if q>0 else '#1f77b4',ec='k',lw=0.5,zorder=5)); a.text(v[0],v[1],'%+d'%q,ha='center',va='center',fontsize=6,color='w',zorder=6)
X=[p[0] for e in s.edges for p in e]; Y=[p[1] for e in s.edges for p in e]; a.set_xlim(min(X)-.5,max(X)+.5); a.set_ylim(min(Y)-.5,max(Y)+.5); a.set_aspect(1); a.axis('off'); a.set_title('(a) Shakti sample, one B–W pair')
b=ax[1]
for c in s.cells:
    b.add_patch(Rectangle(c,1,1,fc='#dddddd' if (c[0]+c[1])%2==0 else 'white',ec='#bbbbbb',lw=0.4))
for pr in TL:
    k=CI[pr]
    cs=sorted(pr); diag=abs(cs[0][0]-cs[1][0])==1 and abs(cs[0][1]-cs[1][1])==1
    if diag:
        col='#111111' if (cs[0][0]+cs[0][1])%2==0 else '#e34948'
        b.plot([cs[0][0]+.5,cs[1][0]+.5],[cs[0][1]+.5,cs[1][1]+.5],color=col,lw=4,solid_capstyle='round')
    else:
        x0=min(c[0] for c in cs); y0=min(c[1] for c in cs); w=abs(cs[0][0]-cs[1][0])+1; h=abs(cs[0][1]-cs[1][1])+1
        b.add_patch(Rectangle((x0+.08,y0+.08),w-.16,h-.16,fc=cols[k%len(cols)],ec='#444',lw=0.8))
b.set_xlim(-4.3,3.3); b.set_ylim(-4.3,3.3); b.set_aspect(1); b.axis('off'); b.set_title(r'(b) the same covering of $\mathcal{A}_3$')
fig.subplots_adjust(left=0.01,right=0.99,top=0.9,bottom=0.02,wspace=0.05)
for ext in ('pdf','png'): fig.savefig(_os.path.join(_ROOT,'figures',f'fig17_shakti_aztec.{ext}'),dpi=200)
print('ok',len(T))
