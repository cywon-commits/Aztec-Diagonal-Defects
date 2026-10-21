import sys, itertools, heapq, json, pickle
from ortools.sat.python import cp_model
# open Shakti sample [ox,ox+L] x [oy,oy+L] cut from the infinite lattice; block (a,b) horizontal iff (a+b) odd
def horiz(a,b): return (a+b)%2==1
def present(e):
    (x1,y1),(x2,y2)=e; a,b=min(x1,x2)//2,min(y1,y2)//2
    if x1==x2: return not (horiz(a,b) and x1%2==1)
    return not ((not horiz(a,b)) and y1%2==1)
class Sample:
    def __init__(s,x0,y0,x1,y1,J1,J2):
        s.J1,s.J2=J1,J2
        s.edges=[]
        for x in range(x0,x1+1):
            for y in range(y0,y1+1):
                if x<x1 and present(((x,y),(x+1,y))): s.edges.append(((x,y),(x+1,y)))
                if y<y1 and present(((x,y),(x,y+1))): s.edges.append(((x,y),(x,y+1)))
        s.idx={e:i for i,e in enumerate(s.edges)}
        s.VA={}
        for (a,b),i in s.idx.items():
            s.VA.setdefault(a,[]).append(((b[0]-a[0],b[1]-a[1]),i,-1))   # +1 spin points a->b : out of a
            s.VA.setdefault(b,[]).append(((a[0]-b[0],a[1]-b[1]),i,+1))
        s.V=list(s.VA)
        # block centers fully inside
        s.centers={}
        for v in s.V:
            if v[0]%2==1 and v[1]%2==1 and len(s.VA[v])==2 and x0<v[0]<x1 and y0<v[1]<y1:
                s.centers[v]=(v[0]//2,v[1]//2)
        s.nE=len(s.edges)
    def ev(s,v,bits):
        a=s.VA[v]; E=0
        for (k1,(d1,i1,s1)),(k2,(d2,i2,s2)) in itertools.combinations(enumerate(a),2):
            J=s.J1 if (d1[0]==-d2[0] and d1[1]==-d2[1]) else s.J2
            E+=J*s1*s2*bits[k1]*bits[k2]
        return E
    def E(s,S): return sum(s.ev(v,[S[i] for d,i,sg in s.VA[v]]) for v in s.V)
    def build(s,m,tag,spec):
        x=[m.NewBoolVar(f'{tag}{i}') for i in range(s.nE)]; Es=[]
        for v in s.V:
            arms=s.VA[v]; rows=[]
            for bits in itertools.product([0,1],repeat=len(arms)):
                sp=[1 if b else -1 for b in bits]
                q=sum(sg*t for (d,i,sg),t in zip(arms,sp))
                if v in s.centers:
                    want=spec.get(s.centers[v],None)
                    if want is None and q!=0: continue
                    if want is not None and q==0: continue
                rows.append(bits+(s.ev(v,sp),))
            lo=min(r[-1] for r in rows); hi=max(r[-1] for r in rows)
            e=m.NewIntVar(lo,hi,''); m.AddAllowedAssignments([x[i] for d,i,sg in arms]+[e],rows); Es.append(e)
        return x,sum(Es)
    def minE(s,spec,tl=120):
        m=cp_model.CpModel(); x,E=s.build(m,'x',spec); m.Minimize(E)
        so=cp_model.CpSolver(); so.parameters.max_time_in_seconds=tl; so.parameters.num_workers=8
        r=so.Solve(m); return int(so.ObjectiveValue()), so.StatusName(r)
    def minflip(s,specX,EX,specY,EY,tl=300):
        m=cp_model.CpModel(); X,ex=s.build(m,'x',specX); Y,ey=s.build(m,'y',specY)
        m.Add(ex==EX); m.Add(ey==EY)
        d=[m.NewBoolVar('') for i in range(s.nE)]
        for i in range(s.nE): m.AddBoolXOr([X[i],Y[i],d[i].Not()])
        m.Minimize(sum(d)); so=cp_model.CpSolver(); so.parameters.max_time_in_seconds=tl; so.parameters.num_workers=8
        r=so.Solve(m)
        if r not in (cp_model.OPTIMAL,cp_model.FEASIBLE): return None
        return so.StatusName(r),int(so.ObjectiveValue()),[1 if so.Value(v) else -1 for v in X],[1 if so.Value(v) else -1 for v in Y]
    def barrier(s,X,Y):
        fl=[i for i in range(s.nE) if X[i]!=Y[i]]; n=len(fl)
        aff=sorted({v for v in s.V for d,i,sg in s.VA[v] if i in fl})
        def Eaff(S): return sum(s.ev(v,[S[i] for d,i,sg in s.VA[v]]) for v in aff)
        Es=[]
        for mk in range(2**n):
            S=list(X)
            for k,i in enumerate(fl):
                if mk>>k&1: S[i]=-S[i]
            Es.append(Eaff(S))
        E0=Es[0]; Es=[e-E0 for e in Es]
        best={0:0};pq=[(0,0)]
        while pq:
            b,mk=heapq.heappop(pq)
            if b>best[mk]:continue
            for k in range(n):
                nn=mk^(1<<k);nb=max(b,Es[nn])
                if nb<best.get(nn,1e9): best[nn]=nb;heapq.heappush(pq,(nb,nn))
        return Es[-1],best[2**n-1]
if __name__=='__main__':
    x0,y0,x1,y1,J1,J2=map(int,sys.argv[1:7]); pairs=json.loads(sys.argv[7]); out=sys.argv[8]
    S=Sample(x0,y0,x1,y1,J1,J2)
    E0,st=S.minE({}); print('E0',E0,st,flush=True)
    res={}
    for name,blocks in pairs.items():
        spec={tuple(b):0 for b in blocks}
        E1,st1=S.minE(spec)
        R=S.minflip({},E0,spec,E1)
        if R is None: print(name,'none',flush=True); continue
        st2,nf,X,Y=R
        dE,bar=S.barrier(X,Y) if nf<=20 else (E1-E0,None)
        res[name]=dict(E1=E1,nflip=nf,X=X,Y=Y,dE=dE,bar=bar)
        print(f'{name}: dE={(E1-E0)/J1:.2f} J1  flips={nf} ({st2})  barrier={None if bar is None else bar/J1:.2f} J1  annih={None if bar is None else (bar-dE)/J1:.2f} J1',flush=True)
    pickle.dump(dict(res=res,edges=S.edges,box=(x0,y0,x1,y1)),open(out,'wb'))
