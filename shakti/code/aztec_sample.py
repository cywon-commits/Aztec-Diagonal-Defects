"""Open Shakti sample whose plaquettes are exactly the Aztec diamond AD_n (placement ox)."""
import os as _os, sys as _sys
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.path.dirname(_os.path.dirname(_HERE))
_sys.path.insert(0,_HERE); _sys.path.insert(0,_ROOT)

import sys, itertools, math, json

from openmult import plaquettes, touch, same_block
from aztecdefects.lattice import aztec
from ortools.sat.python import cp_model
PS=plaquettes(24)
def cell(p):
    (x0,y0),(x1,y1)=p[1],p[2]; X=(x0+x1)/2-12; Y=(y0+y1)/2-12
    return (math.floor((X+Y)/2), math.floor((Y-X)/2))
CELL={cell(p):p for p in PS}
def unit_edges(p):
    (x0,y0),(x1,y1)=p[1],p[2]; E=[]
    for x in range(x0,x1): E+=[((x,y0),(x+1,y0)),((x,y1),(x+1,y1))]
    for y in range(y0,y1): E+=[((x0,y),(x0,y+1)),((x1,y),(x1,y+1))]
    return E
class AztecSample:
    def __init__(s,n,ox,J1,J2):
        s.cells=aztec(n,ox); s.pl={c:CELL[c] for c in s.cells}; s.J1,s.J2=J1,J2
        s.edges=sorted({e for p in s.pl.values() for e in unit_edges(p)})
        s.VA={}
        for i,(a,b) in enumerate(s.edges):
            s.VA.setdefault(a,[]).append(((b[0]-a[0],b[1]-a[1]),i,-1)); s.VA.setdefault(b,[]).append(((a[0]-b[0],a[1]-b[1]),i,+1))
        s.V=list(s.VA); s.nE=len(s.edges)
        # block centres with both plaquettes present = diagonal sites
        s.sites={}
        for c,p in s.pl.items():
            for d,q in s.pl.items():
                if c<d and same_block(p,q):
                    (x0,y0),(x1,y1)=p[1],p[2]; bx,by=(x0//2)*2+1,(y0//2)*2+1
                    s.sites[(bx,by)]=(c,d)
    def ev(s,v,sp):
        E=0
        for ((d1,i1,s1),t1),((d2,i2,s2),t2) in itertools.combinations(zip(s.VA[v],sp),2):
            E+=(s.J1 if (d1[0]==-d2[0] and d1[1]==-d2[1]) else s.J2)*s1*s2*t1*t2
        return E
    def model(s):
        m=cp_model.CpModel(); x=[m.NewBoolVar('') for i in range(s.nE)]; Es=[]; Q={}
        for v in s.V:
            arms=s.VA[v]; rows=[]
            for bits in itertools.product([0,1],repeat=len(arms)):
                sp=[1 if b else -1 for b in bits]; q=sum(sg*t for (d,i,sg),t in zip(arms,sp))
                rows.append(bits+(s.ev(v,sp),int(v in s.sites and q!=0)))
            e=m.NewIntVar(min(r[-2] for r in rows),max(r[-2] for r in rows),''); ch=m.NewBoolVar('')
            m.AddAllowedAssignments([x[i] for d,i,sg in arms]+[e,ch],rows); Es.append(e)
            if v in s.sites: Q[v]=ch
        return m,x,sum(Es),Q
def spectrum(n,ox,J1,J2,levels=2,tl=900):
    tl=int(_os.environ.get("ENUM_TL",tl))
    s=AztecSample(n,ox,J1,J2); out=[]
    m,x,E,Q=s.model(); sv=cp_model.CpSolver(); sv.parameters.num_workers=8; sv.parameters.max_time_in_seconds=120
    m.Minimize(E); sv.Solve(m); E0=int(sv.ObjectiveValue())
    for kk in range(levels):
        # count states with exactly 2kk charged diagonal sites at the lowest energy for that k
        m,x,E,Q=s.model(); m.Add(sum(Q.values())==2*kk); m.Minimize(E)
        sv=cp_model.CpSolver(); sv.parameters.num_workers=8; sv.parameters.max_time_in_seconds=120
        st=sv.Solve(m)
        if st not in (cp_model.OPTIMAL,cp_model.FEASIBLE): out.append((kk,None,0)); continue
        Ek=int(sv.ObjectiveValue())
        m,x,E,Q=s.model(); m.Add(sum(Q.values())==2*kk); m.Add(E==Ek)
        class C(cp_model.CpSolverSolutionCallback):
            def __init__(c): super().__init__(); c.n=0; c.pos={}
            def on_solution_callback(c):
                c.n+=1; key=tuple(sorted(v for v,q in Q.items() if c.Value(q)))
                c.pos[key]=c.pos.get(key,0)+1
        sv=cp_model.CpSolver(); sv.parameters.enumerate_all_solutions=True; sv.parameters.max_time_in_seconds=tl
        c=C(); st=sv.Solve(m,c); out.append((kk,Ek-E0,c.n,sv.StatusName(st),c.pos))
    return E0,out,s
if __name__=='__main__':
    n,ox=int(sys.argv[1]),int(sys.argv[2]); J1,J2=int(sys.argv[3]),int(sys.argv[4])
    E0,out,s=spectrum(n,ox,J1,J2,levels=int(sys.argv[5]) if len(sys.argv)>5 else 2)
    print(f'AD_{n} ox={ox} J1={J1} J2={J2}: plaquettes {len(s.cells)}, islands {s.nE}, diagonal sites {len(s.sites)}, E0={E0}')
    for kk,dE,N,*rest in out:
        print(f'  k={kk}: lowest dE={dE} (units of J1/{J1}), #states={N}', rest[0] if rest else '')
        if rest and kk==1: print('   per-position counts (first 12):',sorted(rest[1].items(),key=lambda z:-z[1])[:12])
