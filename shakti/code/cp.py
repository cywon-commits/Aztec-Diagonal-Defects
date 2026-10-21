import sys, itertools
from ortools.sat.python import cp_model
src=open('shakti.py').read().split("J1=1.0; J2=")[0]
ARGS=sys.argv[:]; P=int(ARGS[1]); sys.argv=[ARGS[0],str(P)]
exec(src)
nE=len(edges)
def vtype(v,sp):  # sp dict arm-> in(+1)/out(-1)
    z=len(sp); q=sum(sp.values())
    if z==2: return 'ok' if q==0 else 'O'
    if z==3:
        if abs(q)==3: return 'bad'
        c=[d for d in sp if (-d[0],-d[1]) in sp]; return 'ok' if sp[c[0]]==sp[c[1]] else 'U'
    if q: return 'bad'
    return 'ok' if sp[(1,0)]==sp[(-1,0)] else 'bad'
def build(m,tag,nO_H,nO_V,nU):
    x=[m.NewBoolVar(f'{tag}{i}') for i in range(nE)]   # 1 means spin = S0 reversed? no: 1 -> +1 along edge
    Ulist=[];OH=[];OV=[]
    for v in V:
        arms=VA[v]
        allowed={'U':[],'O':[],'ok':[]}
        for bits in itertools.product([0,1],repeat=len(arms)):
            sp={d:s*(1 if b else -1) for (d,i,s),b in zip(arms,bits)}
            t=vtype(v,sp)
            if t!='bad': allowed[t].append(bits)
        vars_=[x[i] for d,i,s in arms]
        tU=m.NewBoolVar('');tO=m.NewBoolVar('')
        # table with type flags
        rows=[b+(1,0) for b in allowed['U']]+[b+(0,1) for b in allowed['O']]+[b+(0,0) for b in allowed['ok']]
        m.AddAllowedAssignments(vars_+[tU,tO],rows)
        Ulist.append(tU)
        if len(arms)==2:
            a,b=v[0]//2,v[1]//2
            (OH if horiz(a,b) else OV).append(tO)
    m.Add(sum(OH)==nO_H); m.Add(sum(OV)==nO_V); m.Add(sum(Ulist)==nU)
    return x
nplaq=len([v for v in V if len(VA[v])==2])*2
m=cp_model.CpModel()
kH,kV=int(ARGS[2]),int(ARGS[3])
X=build(m,'x',0,0,nplaq//2)
Y=build(m,'y',kH,kV,(nplaq-2*(kH+kV))//2)
d=[m.NewBoolVar('') for i in range(nE)]
for i in range(nE): m.AddBoolXOr([X[i],Y[i],d[i].Not()])
m.Minimize(sum(d))
s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=float(ARGS[4]) if len(ARGS)>4 else 120; s.parameters.num_workers=8
r=s.Solve(m); print(s.StatusName(r),'min flips',s.ObjectiveValue(),'bound',s.BestObjectiveBound())
import pickle
pickle.dump(dict(X=[s.Value(v) for v in X],Y=[s.Value(v) for v in Y],edges=edges),open(f'sol_{kH}{kV}_{P}.pkl','wb'))
