import sys, itertools, pickle, json, heapq
from ortools.sat.python import cp_model
ARGS=sys.argv[:]; P=int(ARGS[1]); sys.argv=[ARGS[0],str(P)]
exec(open('shakti.py').read().split("J1=1.0; J2=")[0])
nE=len(edges)
def vtype(sp):
    z=len(sp); q=sum(sp.values())
    if z==2: return 'ok' if q==0 else ('Op' if q>0 else 'Om')
    if z==3:
        if abs(q)==3: return 'bad'
        c=[d for d in sp if (-d[0],-d[1]) in sp]; return 'ok' if sp[c[0]]==sp[c[1]] else 'U'
    if q: return 'bad'
    return 'ok' if sp[(1,0)]==sp[(-1,0)] else 'bad'
TABLE={}
for v in V:
    arms=VA[v]; rows=[]
    for bits in itertools.product([0,1],repeat=len(arms)):
        sp={d:s*(1 if b else -1) for (d,i,s),b in zip(arms,bits)}
        t=vtype(sp)
        if t=='bad': continue
        rows.append(bits+(int(t=='U'),int(t=='Op'),int(t=='Om')))
    TABLE[v]=rows
nplaq=sum(2 for v in V if len(VA[v])==2)
def center(b): return ((2*b[0]+1)%P,(2*b[1]+1)%P)
def build(m,tag,spec):
    """spec: dict block->sign (+1,-1, or 0 = either); all other centers uncharged"""
    x=[m.NewBoolVar(f'{tag}{i}') for i in range(nE)]
    U=[]
    cen={center(b):s for b,s in spec.items()}
    for v in V:
        tU,tp,tm=m.NewBoolVar(''),m.NewBoolVar(''),m.NewBoolVar('')
        m.AddAllowedAssignments([x[i] for d,i,s in VA[v]]+[tU,tp,tm],TABLE[v])
        U.append(tU)
        if len(VA[v])==2:
            s=cen.get(v,None)
            if s is None: m.Add(tp+tm==0)
            elif s==0: m.Add(tp+tm==1)
            elif s>0: m.Add(tp==1)
            else: m.Add(tm==1)
    m.Add(sum(U)==(nplaq-2*len(spec))//2)
    return x
def solve(specX,specY,tl=300):
    m=cp_model.CpModel(); X=build(m,'x',specX); Y=build(m,'y',specY)
    d=[m.NewBoolVar('') for i in range(nE)]
    for i in range(nE): m.AddBoolXOr([X[i],Y[i],d[i].Not()])
    m.Minimize(sum(d))
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=tl; s.parameters.num_workers=8
    r=s.Solve(m)
    if r not in (cp_model.OPTIMAL,cp_model.FEASIBLE): return None
    return dict(status=s.StatusName(r),nflip=int(s.ObjectiveValue()),bound=s.BestObjectiveBound(),
                X=[1 if s.Value(v) else -1 for v in X],Y=[1 if s.Value(v) else -1 for v in Y])
def Evert(v,S,J1,J2):
    t=0
    for (d1,i1,s1),(d2,i2,s2) in itertools.combinations(VA[v],2):
        t+=(J1 if d1[0]==-d2[0] and d1[1]==-d2[1] else J2)*s1*s2*S[i1]*S[i2]
    return t
def barrier(X,Y,r):
    fl=[i for i in range(nE) if X[i]!=Y[i]]; n=len(fl)
    aff=sorted({v for v in V for d,i,s in VA[v] if i in fl})
    Es=[]
    for mk in range(2**n):
        S=list(X)
        for k,i in enumerate(fl):
            if mk>>k&1: S[i]=-S[i]
        Es.append(sum(Evert(v,S,1.0,r) for v in aff))
    E0=Es[0]; Es=[e-E0 for e in Es]
    best={0:0.0};pq=[(0.0,0)]
    while pq:
        b,mk=heapq.heappop(pq)
        if b>best[mk]:continue
        for k in range(n):
            nn=mk^(1<<k);nb=max(b,Es[nn])
            if nb<best.get(nn,1e9): best[nn]=nb;heapq.heappush(pq,(nb,nn))
    T=2**n-1
    return round(Es[T],3),round(best[T],3)
if __name__=='__main__':
    queries=json.loads(ARGS[2]); out=ARGS[3]
    res={}
    for name,(sx,sy) in queries.items():
        SX={tuple(map(int,k.split(','))):v for k,v in sx.items()}
        SY={tuple(map(int,k.split(','))):v for k,v in sy.items()}
        R=solve(SX,SY)
        if R is None: res[name]=None; print(name,'INFEASIBLE',flush=True); continue
        R['barrier']={r:barrier(R['X'],R['Y'],r) for r in (1.2,1.5,1.8)} if R['nflip']<=22 else None
        res[name]=R
        print(name,R['status'],R['nflip'],R['barrier'],flush=True)
    pickle.dump(res,open(out,'wb'))
