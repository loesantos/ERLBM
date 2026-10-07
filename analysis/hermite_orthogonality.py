#!/usr/bin/env python3
"""Ortogonalidade ( produto ponderado ) entre o polinomio de Hermite H_xxyy e os de segunda ordem,
em D2Q9, D3Q19 e D3Q27, com aritmetica exata.  No D3Q19, sum w H_xxyy H_zz = -1/27 ( Sec. 2.2 do artigo ).

    python3 analysis/hermite_orthogonality.py
"""
import itertools, numpy as np
from fractions import Fraction as F
def lat(name):
    if name=='D2Q9':
        c=[(x,y) for x in (-1,0,1) for y in (-1,0,1)]
        w=[{0:F(4,9),1:F(1,9),2:F(1,36)}[x*x+y*y] for x,y in c]
    elif name=='D3Q19':
        c=[v for v in itertools.product((-1,0,1),repeat=3) if sum(a*a for a in v)<=2]
        w=[{0:F(1,3),1:F(1,18),2:F(1,36)}[sum(a*a for a in v)] for v in c]
    else:
        c=list(itertools.product((-1,0,1),repeat=3))
        w=[{0:F(8,27),1:F(2,27),2:F(1,54),3:F(1,216)}[sum(a*a for a in v)] for v in c]
    return c,w
cs2=F(1,3)
for name in ['D2Q9','D3Q19','D3Q27']:
    c,w=lat(name); d=len(c[0])
    H2=lambda v,a,b: v[a]*v[b]-(cs2 if a==b else 0)
    Hxxyy=lambda v: (v[0]**2-cs2)*(v[1]**2-cs2)
    out=[]
    for a in range(d):
        for b in range(a,d):
            s=sum(wi*Hxxyy(v)*H2(v,a,b) for v,wi in zip(c,w))
            out.append(((a,b),s))
    # also third-order Hermite x second order
    H3=lambda v: v[0]**2*v[1]-cs2*v[1]
    s3=[sum(wi*H3(v)*H2(v,a,b) for v,wi in zip(c,w)) for a in range(d) for b in range(d)]
    print(name,[(k,str(s)) for k,s in out if s!=0], 'H3xH2 nonzero:',any(x!=0 for x in s3))
