"""One dimension-12 comparison driver, one C++ invariant engine.

python compare_all.py
Records every exact parameter tuple, binary field basis, table SHA256, ANF
degree, APN check, and ortho-derivative histogram in comparison_manifest.json.
Equal histograms are recorded as unresolved, NEVER as equivalence proofs.
"""
from pathlib import Path
import argparse, hashlib, json, re, subprocess, tempfile
from collections import Counter
import numpy as np
from itertools import chain

ROOT = Path(__file__).resolve().parent
TARGET = [[0,10173780],[2,5090166],[4,1265283],[6,213372],[8,27441],
          [10,2727],[12,234],[14,63],[16,18],[18,36]]


class BinaryField:
    def __init__(self, degree, modulus):
        self.degree, self.size, self.modulus = degree, 1 << degree, modulus
        a = np.arange(self.size, dtype=np.int64)
        assert np.all(self.power(a, self.size) == a)
        assert np.all(self.power(a[1:], self.size-1) == 1)

    def mul(self, a, b):
        a, b = np.broadcast_arrays(np.asarray(a,dtype=np.int64),
                                   np.asarray(b,dtype=np.int64))
        a, b = a.copy(), b.copy()
        r = np.zeros_like(a)
        for _ in range(self.degree):
            r ^= np.where(b & 1, a, 0)
            b >>= 1
            a <<= 1
            a ^= np.where(a & self.size, self.modulus, 0)
        return r

    def power(self, a, exponent):
        a = np.asarray(a,dtype=np.int64).copy()
        r = np.ones_like(a)
        while exponent:
            if exponent & 1:
                r = self.mul(r,a)
            a = self.mul(a,a)
            exponent >>= 1
        return r

    def trace(self, a, subdegree=1):
        r = np.zeros_like(np.asarray(a,dtype=np.int64))
        t = np.asarray(a,dtype=np.int64)
        for _ in range(self.degree//subdegree):
            r ^= t
            t = self.power(t,1 << subdegree)
        return r


def degree(vals):
    c = vals.astype(np.uint32).copy()
    n = (len(c)-1).bit_length()
    for i in range(n):
        step = 1 << i
        block = c.reshape(-1,2*step)
        block[:,step:] ^= block[:,:step]
    return max((int(z).bit_count() for z in np.flatnonzero(c)),default=0)


def packed(coords, m):
    result = np.zeros_like(coords[0],dtype=np.uint32)
    for i,c in enumerate(coords):
        result ^= np.asarray(c,dtype=np.uint32) << (i*m)
    return result


def bivariate_clv(k,c):
    m=k.degree; zz=np.arange(k.size*k.size);x=zz%k.size;y=zz//k.size
    M,P=k.mul,k.power; t=M(x,y)
    p=P(x,3)^M(x,P(y,2))^M(c,P(y,3))^t
    q=P(x,5)^M(c,M(P(x,4),y))^M(P(c^1,2),M(x,P(y,4)))^M(c,P(y,5))^t^M(c,P(t,2))
    g=p^M(P(c,k.size-2),k.trace(M(c,t)))
    return packed([p,q],m),packed([g,q],m)


def cases():
    k=BinaryField(6,0x43);M,P=k.mul,k.power
    zz=np.arange(4096);x=zz%64;y=zz//64
    source='Calderini-Li-Villa, arXiv:2204.07462v2 Theorem 4'
    for c in range(1,64):
        if np.any(P(np.arange(64),3)^np.arange(64)^c == 0):continue
        f,g=bivariate_clv(k,c)
        yield 'CLV',dict(alpha=c),k,source,f
        yield 'Trace_CLV',dict(alpha=c),k,'Trace_Bivariate_APN_revised.tex',g
    source='Gologlu-Kolsch, arXiv:2111.04197 Table 2.1 and Theorem 3.3'
    noncubes=[a for a in range(1,64) if int(P(a,21))!=1]
    sub8=[a for a in range(1,64) if int(P(a,8))==a]
    for i in [1,5]:
        q=1<<i
        xq,yq=P(x,q),P(y,q)
        for d in range(1,64):
            if not np.any(P(np.arange(64),q+1)^np.arange(64)^d == 0):
                # (1,0,1,d)_q, (0,0,1,0)_(q^2)
                f=packed([M(xq,x)^M(x,yq)^M(d,M(yq,y)),M(x,P(y,q*q))],6)
                yield 'Taniguchi_normalized',dict(i=i,d=d),k,source,f
        for d in noncubes:
            for j in [0,2,4]:
                f=packed([M(xq,x)^M(d,M(yq,y)),M(x,P(y,1<<j))],6)
                yield 'Zhou_Pott_even_j',dict(i=i,j=j,d=d),k,source,f
        r=1<<(i+3)
        xr,yr=P(x,r),P(y,r)
        for b in noncubes:
            for a in sub8:
                if int(P(b,q+r))==int(P(a,q+1)):continue
                f=packed([M(xq,x)^M(b,M(yq,y)),M(xr,y)^M(M(a,P(b,62)),M(x,yr))],6)
                yield 'Gologlu_Kolsch_F4',dict(i=i,B=b,a=a),k,source,f
    k=BinaryField(4,0x13);M,P=k.mul,k.power
    zz=np.arange(4096);x=zz%16;y=(zz//16)%16;z=zz//256
    source='Bartoli-Stanica arXiv:2603.15146 coordinate definitions; finite APN test only at m=4'
    for i in [1,3]:
        q=1<<i;xq,yq,zq=P(x,q),P(y,q),P(z,q)
        for a in range(1,16):
            f=packed([M(xq,x)^M(a,M(xq,z))^M(y,zq),M(xq,z)^M(yq,y),M(x,yq)^M(a,M(yq,z))^M(zq,z)],4)
            h=packed([M(xq,x)^M(a,M(x,yq))^M(y,zq),M(x,yq)^M(zq,z),M(xq,z)^M(yq,y)^M(a,M(yq,z))],4)
            for name,vals in [('Bartoli_Stanica_G',f),('Bartoli_Stanica_H',h)]:
                yield name,dict(i=i,a=a,Li_Kaleyski_specialization=(a==1)),k,source,vals
    k=BinaryField(3,0x0b);M,P=k.mul,k.power
    zz=np.arange(4096);coords=[(zz>>(3*j))&7 for j in range(4)]
    for i in [1,2]:
        q=1<<i;outs=[]
        for shift in range(4):
            x,y,z,w=[coords[(j+shift)%4] for j in range(4)]
            outs.append(M(P(x,q),y)^M(P(x,q),w)^P(y,q+1)^M(P(y,q),z)
                        ^M(P(z,q),x)^M(P(z,q),w)^M(P(w,q),x))
        yield 'Cyclic_quadrivariate_uploaded',dict(i=i),k,'Uploaded Quadrivariate_APN_family_CCZ(1).tex first-coordinate formula',packed(outs,3)
    k=BinaryField(12,0x1053);M,P=k.mul,k.power;x=np.arange(4096)
    for i in [1,5,7,11]:
        yield 'Gold',dict(i=i),k,'Gold power family, gcd(i,12)=1',P(x,(1<<i)+1)
    yield 'Cubic_plus_trace',{},k,'Li-Kaleyski ePrint 2022/1522 Table II F4, a=1',P(x,3)^k.trace(P(x,9))
    primitive=next(a for a in range(2,4096) if all(int(P(a,4095//p))!=1 for p in [3,5,7,13]))
    for t in [0,1,2]:
        a=int(P(primitive,t)); inv=int(P(a,4094))
        for shift in [0,1]:
            inside=M(P(a,3*(1<<shift)),P(x,9*(1<<shift)))^M(P(a,6*(1<<shift)),P(x,18*(1<<shift)))
            vals=P(x,3)^M(inv,k.trace(inside,3))
            yield 'Budaghyan_Carlet_trace3',dict(a=a,shift=shift,cube_class=t),k,'Li-Kaleyski ePrint 2022/1522 Table II F5,F6',vals
    sub16=[a for a in range(4096) if int(P(a,16))==a]
    u=primitive
    for s in [5,11]:
        base=M(u,P(x,(1<<s)+1))^M(P(u,16),P(x,(1<<8)+(1<<((4+s)%12))))
        vterm=P(x,(1<<8)+1)
        wterm=M(P(u,17),P(x,(1<<s)+(1<<((4+s)%12))))
        for v in sub16:
            for w in sub16:
                if int(M(v,w))==1:continue
                vals=base^M(v,vterm)^M(w,wterm)
                yield 'BBMM_fixed_primitive',dict(s=s,u=u,v=v,w=w),k,'Li-Kaleyski ePrint 2022/1522 Table II F7-F9; fixed primitive u only',vals


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--families',nargs='*',help='Optional exact registry family names')
    parser.add_argument('--external-manifest',type=Path,help='JSON list of exact raw-table descriptors, for additional comparisons')
    parser.add_argument('--output-prefix',default='comparison',help='Prefix of result files; useful for separate deterministic batches')
    args=parser.parse_args()
    engine=ROOT/'invn_canonical'
    subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'invn_canonical.cpp'),'-o',str(engine)],check=True)
    entries=[]
    with tempfile.TemporaryDirectory() as tmp:
        raw=Path(tmp)/'table.raw'
        def external_cases():
            if not args.external_manifest:return
            for item in json.loads(args.external_manifest.read_text()):
                k=BinaryField(item['field_degree'],item['field_modulus'])
                path=args.external_manifest.parent/item['raw_table']
                vals=np.fromfile(path,dtype='<u4')
                yield item['family'],item['parameters'],k,item['source'],vals
        for family,params,k,source,vals in chain(cases(),external_cases()):
            if args.families and family not in args.families:continue
            vals=np.asarray(vals,dtype='<u4')
            assert len(vals)==4096 and np.all(vals<4096)
            d=degree(vals);assert d<=2,(family,params,d)
            raw.write_bytes(vals.tobytes())
            run=subprocess.run([str(engine),'12',str(raw),'0'],capture_output=True,text=True,check=True)
            apn='APN=1' in run.stdout
            hist=[]
            if apn:
                body=re.search(r'ods\{([^}]*)\}',run.stdout).group(1)
                hist=[[int(j),int(count)] for j,count in re.findall(r'(\d+):(\d+)',body)]
                assert sum(count for _,count in hist)==4095*4096
                assert sum(j*count for j,count in hist)==4095*4096
            item=dict(family=family,parameters=params,field_degree=k.degree,
                      field_modulus=k.modulus,binary_encoding='polynomial coefficients; packed coordinates least significant first',
                      source=source,table_sha256=hashlib.sha256(vals.tobytes()).hexdigest(),
                      algebraic_degree=d,APN=apn,ortho_histogram=hist,
                      comparison=('same_invariant_unresolved' if hist==TARGET else 'separated') if apn else 'not_APN',
                      engine_output=run.stdout.strip())
            entries.append(item)
            if len(entries)%50==0:print('checked',len(entries),'tables;',family,flush=True)
    manifest=dict(dimension=12,target_histogram=TARGET,engine='invn_canonical.cpp',
                  scope_notes=['BBMM covers all subfield v,w at one primitive u, not all primitive u.',
                               'Budaghyan-Carlet covers three cube classes for each of the two trace formulas.',
                               'Bartoli-Stanica formulas are evaluated at even base degree 4; no odd-degree theorem is invoked.',
                               'Gologlu F1,F2 have no admissible exponent at bivariate base degree 6.',
                               'The eight cyclic quadrivariate types and the six Gold-switching tables were not supplied; pending exact builders or tables.'],
                  entries=entries)
    path=ROOT/(args.output_prefix+'_manifest.json')
    path.write_text(json.dumps(manifest,indent=2)+'\n')
    summary={}
    for family in sorted({v['family'] for v in entries}):
        group=[v for v in entries if v['family']==family]
        summary[family]=dict(tables=len(group),APN=sum(v['APN'] for v in group),
                            separated=sum(v['comparison']=='separated' for v in group),
                            same_invariant=sum(v['comparison']=='same_invariant_unresolved' for v in group))
    (ROOT/(args.output_prefix+'_summary.json')).write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
