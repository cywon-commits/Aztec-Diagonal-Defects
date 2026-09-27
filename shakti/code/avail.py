import numpy as np
exec(open('cyc.py').read().split('dfs([start,b2]')[0])
dfs([start,b2],False,6)
cyc=res[6]
Nq=1024
th=(np.arange(Nq)+0.5)*2*np.pi/Nq
T,F=np.meshgrid(th,th,indexing='ij')
den=2*np.cos(T)+2j*np.cos(F)
cache={}
def C(x,y):  # K^{-1}(w-b=(x,y))
    if (x,y) not in cache: cache[(x,y)]=np.mean(np.exp(1j*(x*T+y*F))/den)
    return cache[(x,y)]
def K(b,w):
    dx,dy=w[0]-b[0],w[1]-b[1]; return 1 if dy==0 else 1j
def prob(doms):
    bs=[d[0] if sum(d[0])%2==0 else d[1] for d in doms]; ws=[d[1] if sum(d[0])%2==0 else d[0] for d in doms]
    M=np.array([[K(bs[i],ws[i])*C(ws[j][0]-bs[i][0],ws[j][1]-bs[i][1]) for j in range(len(bs))] for i in range(len(bs))])
    return abs(np.linalg.det(M))
print('P(single horizontal domino)=',prob([((0,0),(1,0))]),' vertical',prob([((0,0),(0,1))]))
tot=0
for c in cyc:
    L=len(c); old=[(c[k],c[(k+1)%L]) for k in range(1,L,2)]
    p=prob(old); tot+=p; print(c,'old dominoes',old,'P=%.5f'%p)
print('sum over 8 moves per B site =',tot)
print('--- pinwheels')
for ox in (0,1):
  for oy in (0,1):
    c=[(ox,oy),(ox-1,oy),(ox-1,oy-1),(ox,oy-1)]
    for chi in (1,-1):
        arms=[(0,1),(-1,0),(0,-1),(1,0)] if chi==1 else [(1,0),(0,1),(-1,0),(0,-1)]
        doms=[(c[k],(c[k][0]+arms[k][0],c[k][1]+arms[k][1])) for k in range(4)]
        cells={x for d in doms for x in d}
        ok=all(diag(x) in cells for x in cells)
        if ok: print('vertex',(ox,oy),'chi',chi,'P=%.6f'%prob(doms))
