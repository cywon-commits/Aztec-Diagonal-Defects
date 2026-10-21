import numpy as np, itertools, sys
exec(open('v.py').read().split('import itertools, heapq')[0].split('for v in [(x,y)')[0])
# R = L-ground-state unit cell (4x4) from user's figure
P=int(sys.argv[1]) if len(sys.argv)>1 else 12
def horiz(a,b): return (a+b)%2==1
def edge_present(e):
    (x1,y1),(x2,y2)=e
    a,b=(min(x1,x2)%P)//2,(min(y1,y2)%P)//2
    if x1==x2: # vertical edge
        return not (horiz(a,b) and x1%2==1)
    else:
        return not ((not horiz(a,b)) and y1%2==1)
spin={}
for x in range(P):
    for y in range(P):
        for e in [((x,y),(x+1,y)),((x,y),(x,y+1))]:
            if not edge_present(e): continue
            ex=((e[0][0]%4,e[0][1]%4),(e[1][0]-e[0][0]+e[0][0]%4, e[1][1]-e[0][1]+e[0][1]%4))
            spin[e]=R[ex]
edges=list(spin)
idx={e:i for i,e in enumerate(edges)}
def arms_of(v):
    x,y=v; out=[]
    for d,e,s in [((1,0),((x,y),(x+1,y)),-1),((0,1),((x,y),(x,y+1)),-1),
                  ((-1,0),(((x-1)%P,y),((x-1)%P+1,y)),1),((0,-1),((x,(y-1)%P),(x,(y-1)%P+1)),1)]:
        if e in idx: out.append((d,idx[e],s))
    return out
V=[(x,y) for x in range(P) for y in range(P)]
VA={v:arms_of(v) for v in V}
S0=np.array([spin[e] for e in edges])
def evert(v,S,J1,J2):
    a=VA[v];E=0
    for (d1,i1,s1),(d2,i2,s2) in itertools.combinations(a,2):
        col=d1[0]==-d2[0] and d1[1]==-d2[1]
        E=E+(J1 if col else J2)*(s1*S[...,i1])*(s2*S[...,i2])
    return E
def classify(v,S):
    a=VA[v]; z=len(a); q=sum(s*S[i] for d,i,s in a)
    return z,q
# sanity: ground state check
from collections import Counter
c=Counter()
for v in V:
    z,q=classify(v,S0)
    c[(z,int(q))]+=1
print('vertex (z,q) counts',dict(c))
J1=1.0; J2=float(sys.argv[2]); x0,y0,x1,y1=map(int,sys.argv[3:7])
win=[i for i,e in enumerate(edges) if all(x0<=p[0]<=x1 and y0<=p[1]<=y1 for p in e)]
N=len(win); print('N free islands',N)
M=np.arange(2**N,dtype=np.int64)
S=np.tile(S0,(2**N,1)).astype(np.int8)
for k,i in enumerate(win):
    S[:,i]=np.where((M>>k)&1,-S0[i],S0[i])
aff=sorted({v for v in V for d,i,s in VA[v] if i in win})
E=np.zeros(2**N)
for v in aff: E+=evert(v,S.astype(float),J1,J2)
E-=E[0]
# centers of blocks (z=2 collinear) charges
centers=[v for v in aff if len(VA[v])==2]
Q2=np.zeros(2**N,dtype=int); qs={}
for v in centers:
    q=sum(s*S[:,i].astype(int) for d,i,s in VA[v]); qs[v]=q; Q2+= (np.abs(q)==2)
# other bad vertices: charged z4, |q|=3 z3
# local minima
nbrE=np.stack([E[M^(1<<k)] for k in range(N)],1)
locmin=np.all(nbrE>E[:,None]-1e-9,1)
# bottleneck barrier from ground via union-find in energy order
order=np.argsort(E,kind='stable'); par=np.arange(2**N); added=np.zeros(2**N,bool)
def f(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
barrier=np.full(2**N,np.inf)
comp_members={}
for m in order:
    added[m]=True
    for k in range(N):
        n=m^(1<<k)
        if added[n]:
            ra,rb=f(m),f(n)
            if ra!=rb:
                par[ra]=rb
                # merging: if one contains ground
                pass
    # record lazily below
# second pass: barrier = min level connecting to ground -> do threshold sweep on unique levels
levels=np.unique(np.round(E,6))
barrier=np.full(2**N,np.inf)
for lv in levels:
    ok=E<=lv+1e-9
    # BFS from ground within ok
    seen=np.zeros(2**N,bool); seen[0]=True; frontier=np.array([0])
    while frontier.size:
        nxt=(frontier[:,None]^(1<<np.arange(N))[None,:]).ravel()
        nxt=nxt[ok[nxt]&~seen[nxt]]; nxt=np.unique(nxt); seen[nxt]=True; frontier=nxt
    newly=seen&np.isinf(barrier); barrier[newly]=lv
    if not np.isinf(barrier[locmin&(Q2>0)]).any(): break
sel=np.where(locmin&(Q2>0))[0]
rows=sorted({(int(Q2[m]),round(E[m],3),round(barrier[m],3),bin(m).count('1')) for m in sel})
print('nO  E_defect  barrier  #flips')
for r in rows[:25]: print(r)
