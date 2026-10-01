"""Audit the uniform cubic construction using the shared field builders.
The proof is in the paper. This checks chosen representatives at m=3,4,5,6,
including even m and non-unit parameters, and directly checks cubic APN.
"""
from pathlib import Path
import json
import numpy as np
from compare_all import BinaryField, bivariate_clv, degree

def parity_table(size):
    return np.array([i.bit_count() & 1 for i in range(size)], dtype=np.uint8)

def radical_basis(h):
    n=(len(h)-1).bit_length()
    rows=[sum(int(h[(1<<i)^(1<<j)]^h[1<<i]^h[1<<j]^h[0])<<j
              for j in range(n)) for i in range(n)]
    rank=0; pivots=[]
    for c in range(n):
        pivot=next((i for i in range(rank,n) if (rows[i]>>c)&1),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        for i in range(n):
            if i!=rank and ((rows[i]>>c)&1):rows[i]^=rows[rank]
        pivots.append(c);rank+=1
    basis=[]
    for c in range(n):
        if c in pivots:continue
        v=1<<c
        for i,p in enumerate(pivots):
            if (rows[i]>>c)&1:v^=1<<p
        basis.append(v)
    return rank,basis

def main():
    results=[]
    for m,modulus in [(3,0xb),(4,0x13),(5,0x25),(6,0x43)]:
        field=BinaryField(m,modulus);elems=np.arange(field.size)
        admissible=[a for a in range(1,field.size)
                    if not np.any((field.power(elems,3)^elems^a)==0)]
        selected=admissible if m<=4 else ([1,6,7] if m==5 else [6,14,34,40])
        for alpha in selected:
            _,g=bivariate_clv(field,alpha)
            size=len(g);n=2*m;z=np.arange(size);parity=parity_table(size)
            for beta in range(1,size):
                h=parity[g & beta]
                rank,basis=radical_basis(h)
                if 4<=rank<=n-2:break
            else:raise AssertionError('required component absent')
            rad=[0]
            for e in basis:rad += [r^e for r in rad]
            v=next(r for r in rad if r and not h[r])
            assert np.all(h[z^v]==h)
            phi=z ^ np.where(h,v,0)
            assert np.all(phi[phi]==z)
            cubic=g[phi]
            polar=g[z^v]^g[v]^g^g[0]
            formula=g ^ np.where(h,polar^g[v],0).astype(np.uint32)
            assert np.array_equal(cubic,formula)
            assert degree(cubic)==3
            assert np.array_equal(np.bincount(g,minlength=size),
                                  np.bincount(cubic,minlength=size))
            max_diff=0
            for a in range(1,size):
                row=np.bincount(cubic[z^a]^cubic,minlength=size)
                max_diff=max(max_diff,int(row.max()))
            assert max_diff==2
            results.append(dict(m=m,modulus=modulus,alpha=alpha,
                binary_dot_product_mask=beta,component_polar_rank=rank,
                radical_dimension=len(basis),direction=v,degree=3,
                differential_uniformity=max_diff,image_fibers_preserved=True))
            print(results[-1],flush=True)
    Path('cubic_extension_audit_results.json').write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':main()
