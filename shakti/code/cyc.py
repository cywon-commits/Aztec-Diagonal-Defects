# cells (i,j); black if i+j even. B: (2m,2n)-(2m+1,2n-1); W: (2m,2n+1)-(2m-1,2n)
def diag(c):
    i,j=c
    if (i+j)%2==0:  # black -> B partner
        return (i+1,j-1) if i%2==0 else (i-1,j+1)
    else:
        return (i-1,j-1) if i%2==0 else (i+1,j+1)
def nb(c):
    i,j=c; return [(i+1,j),(i-1,j),(i,j+1),(i,j-1)]
# alternating cycle: new edges at even positions; exactly one B and one W among new edges.
start=(0,0); b2=diag(start)   # B edge as first new edge
res={}
def dfs(path,usedW,maxlen):
    L=len(path)  # path of cells; edges between consecutive; edge index L-2
    cur=path[-1]; k=L-1  # next edge index
    if L>=4 and k%2==1 and False: pass
    if L-1>=maxlen: return
    new = (k%2==0)
    cands=[(n,'d') for n in nb(cur)]
    if new and not usedW and (sum(cur)%2==1): cands.append((diag(cur),'W'))
    for n,t in cands:
        if n==start and k%2==1 and usedW:
            res.setdefault(L,[]).append(list(path)); continue
        if n in path: continue
        dfs(path+[n],usedW or t=='W',maxlen)
dfs([start,b2],False,12)
for L in sorted(res): print('cycle cells',L,'count',len(res[L])); print(res[L][0])
