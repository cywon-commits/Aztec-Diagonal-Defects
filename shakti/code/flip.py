import sys
sys.argv=['x','8','{}','/dev/null']
exec(open('cp2.py').read().split("if __name__")[0])
from ortools.sat.python import cp_model
def buildU(m,tag):
    x=[m.NewBoolVar(f'{tag}{i}') for i in range(nE)]; U={}
    for v in V:
        tU,tp,tm=m.NewBoolVar(''),m.NewBoolVar(''),m.NewBoolVar('')
        m.AddAllowedAssignments([x[i] for d,i,s in VA[v]]+[tU,tp,tm],TABLE[v]); U[v]=tU
        if len(VA[v])==2: m.Add(tp+tm==0)
    m.Add(sum(U.values())==nplaq//2); return x,U
m=cp_model.CpModel(); X,UX=buildU(m,'x'); Y,UY=buildU(m,'y')
d=[m.NewBoolVar('') for i in range(nE)]
for i in range(nE): m.AddBoolXOr([X[i],Y[i],d[i].Not()])
dv=[m.NewBoolVar('') for v in V]
for b,v in zip(dv,V): m.AddBoolXOr([UX[v],UY[v],b.Not()])
m.Add(sum(dv)>=1)
m.Minimize(sum(d)); s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=200; s.parameters.num_workers=8
r=s.Solve(m); Xs=[1 if s.Value(v) else -1 for v in X]; Ys=[1 if s.Value(v) else -1 for v in Y]
print(s.StatusName(r),'min flips',s.ObjectiveValue(),'unhappy changed',sum(s.Value(b) for b in dv))
print({rr:barrier(Xs,Ys,rr) for rr in (1.2,1.5,1.8)})
