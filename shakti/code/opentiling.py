import os as _os, sys as _sys
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.path.dirname(_os.path.dirname(_HERE))
_sys.path.insert(0,_HERE); _sys.path.insert(0,_ROOT)
import sys, math

from aztec_sample import AztecSample, cell, CELL
from ortools.sat.python import cp_model
import itertools
def vinfo(s,v,S):
    arms=s.VA[v]; sp={d:sg*S[i] for d,i,sg in arms}; return arms,sp
def tiling_open(s,S):
    """Pairs of cells. U vertex: z=3, collinear pair head-to-tail."""
    inv={}
    for c,p in s.pl.items(): inv[((p[1][0]+p[2][0])/2,(p[1][1]+p[2][1])/2)]=c
    pairs=[]
    for v in s.V:
        arms,sp=vinfo(s,v,S); z=len(arms); q=sum(sp.values())
        if v in s.sites and q!=0: pairs.append(frozenset(s.sites[v])); continue
        if z!=3: continue
        dirs={d:(i,sg) for d,i,sg in arms}
        perp=[d for d in dirs if (-d[0],-d[1]) not in dirs]
        if len(perp)!=1: continue
        perp=perp[0]; col=[d for d in dirs if d!=perp]
        if sp[col[0]]==sp[col[1]]: continue      # happy
        i,sg=dirs[perp]; sin=sg*S[i]
        n=(-perp[1],perp[0]); ic,sc=dirs[col[0]]; f=(-col[0][0]*sc*S[ic], -col[0][1]*sc*S[ic])
        k=-sin*(1 if f[0]*n[0]+f[1]*n[1]>0 else -1)
        a=inv.get((v[0]-0.5*perp[0], v[1]-0.5*perp[1])); b=inv.get((v[0]+perp[0]+0.5*k*n[0], v[1]+perp[1]+0.5*k*n[1]))
        if a is not None and b is not None: pairs.append(frozenset([a,b]))
    return pairs
if __name__=='__main__':
    n,ox=int(sys.argv[1]),int(sys.argv[2]); s=AztecSample(n,ox,2,3)
    m,x,E,Q=s.model(); sv=cp_model.CpSolver(); m.Minimize(E); sv.Solve(m); E0=int(sv.ObjectiveValue())
    for kk,dE in [(0,0),(1,4)]:
        m,x,E,Q=s.model(); m.Add(sum(Q.values())==2*kk); m.Add(E==E0+dE)
        res={'ok':0,'bad':0}; tilings=set()
        class C(cp_model.CpSolverSolutionCallback):
            def __init__(c): super().__init__()
            def on_solution_callback(c):
                S=[1 if c.Value(v) else -1 for v in x]; T=tiling_open(s,S)
                cov={a for p in T for a in p}
                good= len(cov)==2*len(T)==len(s.cells)
                res['ok' if good else 'bad']+=1
                if good: tilings.add(frozenset(T))
        sv=cp_model.CpSolver(); sv.parameters.enumerate_all_solutions=True; sv.Solve(m,C())
        print(f'AD_{n} ox={ox} k={kk}: states mapped to perfect coverings {res}, distinct coverings {len(tilings)}')
        from aztecdefects.bruteforce import enumerate_coverings
        from aztecdefects.lattice import partner as _partner
        exact = {frozenset(frozenset(t) for t in cov) for cov in enumerate_coverings(s.cells)
                 if sum(1 for t in cov if _partner(*t[0]) == t[1]) == 2 * kk}
        print(f'   set of coverings equals the exact enumeration ({len(exact)} coverings): {tilings == exact}')
