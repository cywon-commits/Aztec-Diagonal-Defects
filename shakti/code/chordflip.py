import os as _os, sys as _sys
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.path.dirname(_os.path.dirname(_HERE))
_sys.path.insert(0,_HERE); _sys.path.insert(0,_ROOT)
import sys, random
sys.argv=['x','12']; 
from multiplicity import *
for seed in range(3):
    S=random_state({},seed); T=tiling(S)
    m,x=same_Uset_model(S)
    d=[m.NewBoolVar('') for i in range(nE)]
    for i in range(nE): m.AddBoolXOr([x[i], m.NewConstant(1 if S[i]>0 else 0), d[i].Not()])
    m.Add(sum(d)>=1); m.Minimize(sum(d)); sv=cp_model.CpSolver(); sv.parameters.num_workers=2; sv.Solve(m)
    Y=[1 if sv.Value(v) else -1 for v in x]; TY=tiling(Y)
    print('seed',seed,'flips',sum(a!=b for a,b in zip(S,Y)),'covering changed:',T!=TY,'tiles changed',len(T-TY),
          'barrier (dE,barrier) per r',{r:barrier(S,Y,r) for r in (1.2,1.5,1.8)},flush=True)
