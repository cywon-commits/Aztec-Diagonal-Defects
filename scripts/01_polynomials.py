"""Exact generating polynomials Z_n(u) = sum_k c_k u^k (Table I, Secs. III, VI, VII).

    python scripts/01_polynomials.py [--nmax 13] [--nmaxB 8]

n <= 13 takes a few minutes; n = 15 and 17 take ~8 and ~20 minutes.
"""
import argparse, json, os, time
from _common import DATA
from aztecdefects.lattice import aztec
from aztecdefects.pfaffian import generating_polynomial
import sympy as sp

ap = argparse.ArgumentParser(); ap.add_argument('--nmax', type=int, default=13); ap.add_argument('--nmaxB', type=int, default=8)
a = ap.parse_args()
path = os.path.join(DATA, 'polynomials.json')
out = json.load(open(path)) if os.path.exists(path) else {'A': {}, 'B': {}}
for n in range(1, a.nmax + 1):
    if str(n) in out['A']:
        continue
    t = time.time(); c = generating_polynomial(aztec(n)); out['A'][str(n)] = [str(x) for x in c]
    json.dump(out, open(path, 'w')); print(f'A n={n}: degree {len(c)-1}  ({time.time()-t:.0f}s)', flush=True)
for n in range(1, a.nmaxB + 1):
    if str(n) not in out['B']:
        out['B'][str(n)] = [str(x) for x in generating_polynomial(aztec(n, 1, 0))]
json.dump(out, open(path, 'w'))
u = sp.symbols('u')
print('\nTable I (placement A):')
for n in range(1, min(a.nmax, 11) + 1):
    c = [int(x) for x in out['A'][str(n)]]
    r = sp.Rational(c[1], c[0]) if len(c) > 1 else 0
    print(f'n={n:2d}  c1/c0 = {r}   Z = {sp.factor(sum(ci*u**k for k, ci in enumerate(c))) if n <= 4 else "degree %d" % (len(c)-1)}')
print('\nhalf-step check Z_2m = 4^m Z_2m-1:',
      all([int(x) for x in out['A'][str(2*m)]] == [4**m*int(x) for x in out['A'][str(2*m-1)]] for m in range(1, a.nmax // 2 + 1)))
