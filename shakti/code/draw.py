import pickle,sys,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ARGS=sys.argv[:]; sys.argv=[ARGS[0],ARGS[1]]
exec(open('shakti.py').read().split("J1=1.0; J2=")[0])
D=pickle.load(open(ARGS[2],'rb')); x0,y0,x1,y1=map(int,ARGS[3:7]); out=ARGS[7]
X=[1 if b else -1 for b in D['X']]; Y=[1 if b else -1 for b in D['Y']]
fig,axs=plt.subplots(1,2,figsize=(10,5.2))
for ax,S,title in [(axs[0],X,ARGS[9] if len(ARGS)>9 else 'before (ground state, L tiles)'),(axs[1],Y,ARGS[8] if len(ARGS)>8 else 'after (O + O)')]:
    for i,e in enumerate(edges):
        (a,b),(c,d)=e
        if not (x0<=a<=x1 and y0<=b<=y1 and x0<=c<=x1 and y0<=d<=y1): continue
        ax.plot([a,c],[b,d],color='0.75',lw=1,zorder=0)
        sx,sy=(c-a)*S[i],(d-b)*S[i]; mx,my=(a+c)/2,(b+d)/2
        fl=X[i]!=Y[i]
        ax.annotate('',xy=(mx+0.32*sx,my+0.32*sy),xytext=(mx-0.32*sx,my-0.32*sy),
            arrowprops=dict(arrowstyle='-|>',lw=2.6 if fl else 1.4,color='k' if fl else '0.35',mutation_scale=14))
    for v in V:
        x,y=v
        if not(x0<=x<=x1 and y0<=y<=y1): continue
        a=VA[v]; ar={dd:s*S[i] for dd,i,s in a}; z=len(a); q=sum(ar.values())
        if z==2 and q: ax.add_patch(plt.Circle(v,0.16,color='tab:red' if q>0 else 'tab:blue',zorder=3)); ax.text(x,y,'%+d'%q,ha='center',va='center',color='w',fontsize=8,zorder=4)
        if z==3:
            c=[dd for dd in ar if (-dd[0],-dd[1]) in ar]
            if ar[c[0]]!=ar[c[1]]: ax.plot(x,y,'o',ms=9,mfc='none',mec='tab:orange',mew=2,zorder=3)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title)
    for s in ax.spines.values(): s.set_visible(False)
fig.text(0.5,0.02,'bold arrows: flipped islands | orange ring: unhappy z=3 vertex | red/blue disk: charged z=2 block center',ha='center',fontsize=9)
plt.tight_layout(rect=(0,0.04,1,0.95)); plt.savefig(out,dpi=130)
