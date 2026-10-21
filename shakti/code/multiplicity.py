"""Number of island configurations per domino tiling in the J1-J2 Shakti vertex model."""
import os as _os, sys as _sys
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.path.dirname(_os.path.dirname(_HERE))
_sys.path.insert(0,_HERE); _sys.path.insert(0,_ROOT)

import sys, random, itertools
P=int(sys.argv[1]); sys.argv=['cp2.py',str(P),'{}','/dev/null']
import os; os.chdir(_HERE)
exec(open('cp2.py').read().split("if __name__=='__main__':")[0])
from ortools.sat.python import cp_model
def md(p): return (p[0]%P, p[1]%P)
def types(S):
    return {v: vtype({d:s*S[i] for d,i,s in VA[v]}) for v in V}
def tiling(S):
    pairs=[]; tX=types(S)
    for v in V:
        t=tX[v]; arms=VA[v]
        if t=='U':
            dirs={d:(i,s) for d,i,s in arms}
            perp=[d for d in dirs if (-d[0],-d[1]) not in dirs][0]
            i,s=dirs[perp]; sin=s*S[i]
            flat=md((v[0]-0.5*perp[0], v[1]-0.5*perp[1]))
            n=(-perp[1],perp[0]); col=[d for d in dirs if d!=perp][0]; ic,sc=dirs[col]
            f=(-col[0]*sc*S[ic], -col[1]*sc*S[ic]); k=-sin*(1 if f[0]*n[0]+f[1]*n[1]>0 else -1)
            pairs.append(frozenset([flat, md((v[0]+perp[0]+0.5*k*n[0], v[1]+perp[1]+0.5*k*n[1]))]))
        elif t in('Op','Om'):
            d=arms[0][0]; n=(-d[1],d[0])
            pairs.append(frozenset([md((v[0]+0.5*n[0],v[1]+0.5*n[1])), md((v[0]-0.5*n[0],v[1]-0.5*n[1]))]))
    return frozenset(pairs)
def same_Uset_model(S):
    tX=types(S); m=cp_model.CpModel(); x=[m.NewBoolVar('') for i in range(nE)]
    for v in V:
        tU,tp,tm=m.NewBoolVar(''),m.NewBoolVar(''),m.NewBoolVar('')
        m.AddAllowedAssignments([x[i] for d,i,s in VA[v]]+[tU,tp,tm],TABLE[v])
        m.Add(tU==int(tX[v]=='U'))
        if len(VA[v])==2: m.Add(tp+tm==int(tX[v] in ('Op','Om')))   # charge sign free
    return m,x
def multiplicity(S):
    T=tiling(S); m,x=same_Uset_model(S)
    class C(cp_model.CpSolverSolutionCallback):
        def __init__(s): super().__init__(); s.n=0; s.nt=0
        def on_solution_callback(s):
            Z=[1 if s.Value(v) else -1 for v in x]; s.n+=1; s.nt+=(tiling(Z)==T)
    sv=cp_model.CpSolver(); sv.parameters.enumerate_all_solutions=True; sv.parameters.max_time_in_seconds=600
    c=C(); sv.Solve(m,c); return c.n,c.nt
def random_state(spec,seed):
    m=cp_model.CpModel(); x=build(m,'x',spec); rng=random.Random(seed)
    m.Maximize(sum(rng.choice([-3,-2,-1,1,2,3])*v for v in x))
    sv=cp_model.CpSolver(); sv.parameters.max_time_in_seconds=20; sv.parameters.random_seed=seed; sv.parameters.num_workers=8
    sv.Solve(m); return [1 if sv.Value(v) else -1 for v in x]
if __name__=='__main__':
    nb=P//2
    specs=[({},'k=0'),({(1,2):0,(2,2):0},'k=1 adjacent'),({(1,2):0,(nb-1,nb-1):0},'k=1 far'),({(1,2):0,(2,2):0,(nb-2,nb-1):0,(nb-1,nb-1):0},'k=2')]
    for spec,lab in specs:
        out=[]
        for seed in range(5):
            S=random_state(spec,seed); T=tiling(S)
            from collections import Counter
            ok=len({p for pr in T for p in pr})==2*len(T)==nb*nb*2
            n,nt=multiplicity(S); out.append((ok,n,nt))
        print(f'P={P} {lab}: (perfect matching, #configs with same U-set, #configs with same domino tiling) =',out,flush=True)
