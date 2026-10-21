import pickle,sys,itertools,heapq
ARGS=sys.argv[:]; sys.argv=[ARGS[0],ARGS[1]]
exec(open('shakti.py').read().split("J1=1.0; J2=")[0])
D=pickle.load(open(ARGS[2],'rb'))
X=[1 if b else -1 for b in D['X']]; Y=[1 if b else -1 for b in D['Y']]
fl=[i for i in range(len(edges)) if X[i]!=Y[i]]
print('flipped islands:',[edges[i] for i in fl])
aff=sorted({v for v in V for d,i,s in VA[v] if i in fl})
def vstate(v,S):
    a=VA[v]; ar={d:s*S[i] for d,i,s in a}; z=len(a); q=sum(ar.values())
    if z==2: return 'z2q%+d'%q if q else 'z2ok'
    if z==3:
        if abs(q)==3: return 'z3q3'
        c=[d for d in ar if (-d[0],-d[1]) in ar]; return 'z3hap' if ar[c[0]]==ar[c[1]] else 'z3UNH'
    if q: return 'z4q%+d'%q
    return 'z4I' if ar[(1,0)]==ar[(-1,0)] else 'z4II'
for v in aff:
    a,b=vstate(v,X),vstate(v,Y)
    if a!=b: print('  ',v,a,'->',b)
def E(S,J1,J2):
    t=0
    for v in aff:
        for (d1,i1,s1),(d2,i2,s2) in itertools.combinations(VA[v],2):
            t+=(J1 if d1[0]==-d2[0] and d1[1]==-d2[1] else J2)*s1*s2*S[i1]*S[i2]
    return t
n=len(fl)
for r in [1.2,1.5,1.8]:
    Es=[]
    for m in range(2**n):
        S=list(X)
        for k,i in enumerate(fl):
            if m>>k&1: S[i]=-S[i]
        Es.append(E(S,1,r))
    E0=Es[0]; Es=[e-E0 for e in Es]
    best={0:0};pq=[(0,0)];prev={}
    while pq:
        b,m=heapq.heappop(pq)
        if b>best[m]:continue
        for k in range(n):
            nn=m^(1<<k);nb=max(b,Es[nn])
            if nb<best.get(nn,1e9): best[nn]=nb;prev[nn]=m;heapq.heappush(pq,(nb,nn))
    T=2**n-1; p=[T]
    while p[-1]: p.append(prev[p[-1]])
    print(f'r={r}: dE={Es[T]:.2f} barrier(create)={best[T]:.2f} (annihilate {best[T]-Es[T]:.2f}) path',[round(Es[q],2) for q in p[::-1]])
