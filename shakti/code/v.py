# islands: key=((x1,y1),(x2,y2)) with x1<=x2,y1<=y2 ; value +1 means pointing toward (x2,y2)
L={}
def H(y,x,d): L[((x,y),(x+1,y))]=1 if d=='>' else -1
def V(x,y,d): L[((x,y),(x,y+1))]=1 if d=='^' else -1
for x,d in zip(range(4),'><><'): H(4,x,d); H(2,x,d); H(0,x,d)
V(0,3,'^');V(0,2,'v'); H(3,0,'<');H(3,1,'>'); V(2,3,'^');V(2,2,'v')
V(3,3,'v');V(3,2,'^'); V(4,3,'^');V(4,2,'v')
V(0,1,'^');V(0,0,'v'); V(1,1,'v');V(1,0,'^'); V(2,1,'^');V(2,0,'v')
H(1,2,'<');H(1,3,'>'); V(4,1,'^');V(4,0,'v')
R=dict(L)
for e in [((1,3),(2,3)),((2,2),(2,3)),((3,2),(3,3)),((1,2),(2,2)),((2,2),(3,2)),((2,1),(2,2)),((1,1),(1,2)),((2,1),(3,1))]:
    R[e]=-R[e]
def vert(S,v):
    arms={}
    for (a,b),s in S.items():
        if v==a: arms[(b[0]-a[0],b[1]-a[1])]=-s   # +1 = pointing into v? s=+1 points to b -> out of a
        if v==b: arms[(a[0]-b[0],a[1]-b[1])]=s    # s=+1 points to b = into v
    return arms
def classify(arms):
    z=len(arms); q=sum(arms.values())
    if z==2:
        return ('z2', q, 'ok' if q==0 else 'charged')
    if z==3:
        dirs=list(arms)
        col=[d for d in dirs if (-d[0],-d[1]) in arms]
        if abs(q)==3: return ('z3',q,'q3')
        s=sum(arms[d] for d in col)
        return ('z3',q,'ok' if s==0 else 'unhappy')
    if z==4:
        if q!=0: return ('z4',q,'charged')
        # type I: each collinear pair head-to-tail
        if arms[(1,0)]+arms[(-1,0)]==0: return ('z4',0,'typeI')
        return ('z4',0,'typeII')
for v in [(x,y) for y in (3,2,1) for x in (1,2,3)]:
    a,b=classify(vert(L,v)),classify(vert(R,v))
    print(v,a,'->',b)
import itertools, heapq
flips=[((1,3),(2,3)),((2,2),(2,3)),((3,2),(3,3)),((1,2),(2,2)),((2,2),(3,2)),((2,1),(2,2)),((1,1),(1,2)),((2,1),(3,1))]
verts=[(x,y) for y in (3,2,1) for x in (1,2,3)]
def Evert(arms,J1,J2):
    d=list(arms);E=0
    for a,b in itertools.combinations(d,2):
        col=(a[0]==-b[0] and a[1]==-b[1])
        E+=(J1 if col else J2)*arms[a]*arms[b]
    return E
def Etot(S,J1,J2): return sum(Evert(vert(S,v),J1,J2) for v in verts)
def state(mask):
    S=dict(L)
    for i,e in enumerate(flips):
        if mask>>i&1: S[e]=-S[e]
    return S
for r in [1.2,1.5,1.8,2.0,2.2]:
    J1,J2=1.0,r
    E=[Etot(state(m),J1,J2) for m in range(256)]
    # bottleneck path 0 -> 255 via single flips
    best={0:E[0]};pq=[(E[0],0)];prev={}
    while pq:
        b,m=heapq.heappop(pq)
        if b>best.get(m,1e9):continue
        for i in range(8):
            n=m^(1<<i);nb=max(b,E[n])
            if nb<best.get(n,1e9): best[n]=nb;prev[n]=m;heapq.heappush(pq,(nb,n))
    path=[255]
    while path[-1]!=0: path.append(prev[path[-1]])
    path=path[::-1]
    print(f"J2/J1={r}: E_O-E_L={E[0]-E[255]:+.2f}  barrier from O={best[255]-E[0]:.2f} from L={best[255]-E[255]:.2f}  path E:",[round(E[m]-E[0],2) for m in path])
