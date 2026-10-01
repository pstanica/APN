"""Exhaustive switching ranks and direct quotient checks using shared builders."""
import json,re,subprocess,tempfile
from pathlib import Path
import numpy as np
from compare_all import BinaryField,bivariate_clv,ROOT


def rank(rows):
    piv={}
    for r in rows:
        r=int(r)
        while r:
            h=r.bit_length()-1
            if h in piv:r^=piv[h]
            else:piv[h]=r;break
    return len(piv)


def component_forms(table,n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    coeff=[int(table[(1<<i)^(1<<j)]^table[1<<i]^table[1<<j]^table[0]) for i,j in pairs]
    return [sum(((w&beta).bit_count()%2)<<p for p,w in enumerate(coeff))
            for beta in range(1<<n)]


subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'switchn.cpp'),'-o',str(ROOT/'switchn')],check=True)
records=[]
with tempfile.TemporaryDirectory() as tmp:
    tmp=Path(tmp)
    for m,reps in [(4,[1,9]),(5,[1]),(6,[6,14,34,40])]:
        k=BinaryField(m,{4:0x13,5:0x25,6:0x43}[m]); n=2*m
        for c in reps:
            f,g=bivariate_clv(k,c)
            inv=int(k.power(c,k.size-2))
            z=np.arange(1<<n);x=z%k.size;y=z//k.size
            h=k.trace(k.mul(c,k.mul(x,y)))
            traceform=sum(int(h[(1<<i)^(1<<j)]^h[1<<i]^h[1<<j]^h[0])<<p
                          for p,(i,j) in enumerate(( (i,j) for i in range(n) for j in range(i+1,n))))
            for label,vals in [('F',f),('G',g)]:
                path=tmp/'input.raw';np.asarray(vals,dtype='<u4').tofile(path)
                prefix=tmp/f'{label}{m}a{c}'
                run=subprocess.run([str(ROOT/'switchn'),str(n),str(path),str(prefix)],capture_output=True,text=True,check=True)
                forms=component_forms(vals,n)
                cforms=[forms[beta] for beta in range(1<<n) if (beta&inv).bit_count()%2==0]
                assert rank(cforms)==n-1
                assert rank(cforms+[traceform])==n
                expected={n-1:(1<<n)-2,n:1} if m>=5 else {n-1:(1<<n)-4,n:3}
                actual={int(a):int(b) for a,b in re.findall(r'(\d+):(\d+)',run.stderr)}
                assert actual==expected,(m,c,label,actual)
                files=list(tmp.glob(prefix.name+'_v*_s*.raw'))
                directions=sorted(int(p.name.split('_v')[1].split('_')[0]) for p in files)
                if m>=5:assert directions==[inv]
                assert inv in directions
                records.append(dict(m=m,alpha=c,map=label,exceptional_directions=directions,
                                    component_form_rank=n-1,rank_with_trace=n,nullity_histogram=actual))
                print(label,m,c,'quotient/rank checks pass',flush=True)
(ROOT/'switching_audit_results.json').write_text(json.dumps(records,indent=2)+'\n')
