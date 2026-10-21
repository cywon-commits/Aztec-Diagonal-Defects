import sys
from cp_open import Sample
L=int(sys.argv[1]); J1,J2=int(sys.argv[2]),int(sys.argv[3]); nb=L//2
S=Sample(0,0,L,L,J1,J2)
E0,_=S.minE({},200); print('E0',E0,flush=True)
c=nb//2-1
refH=(c,c+1) if (2*c+1)%2==1 else (c,c)   # a horizontal block near center
refV=(c+1,c+1) if (2*c+2)%2==0 else (c,c)
refH=[b for b in [(c,c),(c,c+1),(c+1,c)] if (b[0]+b[1])%2==1][0]
refV=[b for b in [(c,c),(c+1,c+1)] if (b[0]+b[1])%2==0][0]
E,_=S.minE({refH:0,refV:0},200) if abs(refH[0]-refV[0])+abs(refH[1]-refV[1])>0 else (None,None)
print('center pair',refH,refV,'dE/J1=',(E-E0)/J1,flush=True)
eps=2*(2*J1-J2)  # bulk single-O cost
for a in range(nb):
    for b in range(nb):
        if 0<a<nb-1 and 0<b<nb-1: continue
        h=(a+b)%2==1; ref=refV if h else refH
        E,st=S.minE({(a,b):0,ref:0},120)
        print(f'block {(a,b)} {"H" if h else "V"}  site cost (J1) = {(E-E0-eps)/J1:.2f}  [{st}]',flush=True)
