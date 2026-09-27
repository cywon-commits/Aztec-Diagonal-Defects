import os as _os, sys as _sys
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.path.dirname(_os.path.dirname(_HERE))
_sys.path.insert(0,_HERE); _sys.path.insert(0,_ROOT)
import sys, itertools

exec(open(_os.path.join(_HERE,'cp_open.py')).read().split("    def minflip")[0])
from ortools.sat.python import cp_model
from functools import lru_cache
def plaquettes(L):
    Ps=[]
    for a in range(L//2):
        for b in range(L//2):
            if horiz(a,b): Ps+= [('H',(2*a,2*b),(2*a+2,2*b+1)),('H',(2*a,2*b+1),(2*a+2,2*b+2))]
            else: Ps+= [('V',(2*a,2*b),(2*a+1,2*b+2)),('V',(2*a+1,2*b),(2*a+2,2*b+2))]
    return Ps
def touch(p,q):
    (x0,y0),(x1,y1)=p[1],p[2]; (u0,v0),(u1,v1)=q[1],q[2]
    ox=min(x1,u1)-max(x0,u0); oy=min(y1,v1)-max(y0,v0)
    return (ox==0 and oy>0) or (oy==0 and ox>0)
def n_matchings(Ps, allow_parallel=False):
    n=len(Ps); adj=[[j for j in range(n) if j!=i and touch(Ps[i],Ps[j]) and (Ps[i][0]!=Ps[j][0] or (allow_parallel and same_block(Ps[i],Ps[j])))] for i in range(n)]
    @lru_cache(None)
    def f(mask):
        if mask==(1<<n)-1: return 1
        i=next(k for k in range(n) if not mask>>k&1)
        return sum(f(mask|1<<i|1<<j) for j in adj[i] if not mask>>j&1)
    return f(0)
def same_block(p,q): return (p[1][0]//2,p[1][1]//2)==(q[1][0]//2,q[1][1]//2)
def count_ground(L,J1,J2,spec={}):
    s=Sample(0,0,L,L,J1,J2); E0,st=s.minE(spec)
    m=cp_model.CpModel(); x,E=s.build(m,'x',spec); m.Add(E==E0)
    class C(cp_model.CpSolverSolutionCallback):
        def __init__(c): super().__init__(); c.n=0
        def on_solution_callback(c): c.n+=1
    so=cp_model.CpSolver(); so.parameters.enumerate_all_solutions=True; so.parameters.max_time_in_seconds=600
    c=C(); so.Solve(m,c); return E0,c.n
if __name__=="__main__":
   for L in (4,6,8):
      Ps=plaquettes(L)
      for J1,J2 in [(5,6),(2,3),(5,9)]:
          E0,n=count_ground(L,J1,J2)
          print(f'L={L} r={J2/J1:.2f}: ground states {n}, domino tilings of plaquette region {n_matchings(Ps)}, E0={E0}',flush=True)
