"""Reproduce every stated permutation exclusion with the shared builders/engine."""
from pathlib import Path
import json,subprocess,tempfile
import numpy as np
from compare_all import ROOT,BinaryField,bivariate_clv

subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'invn_canonical.cpp'),'-o',str(ROOT/'invn_canonical')],check=True)
results=[]
with tempfile.TemporaryDirectory() as tmp:
    path=Path(tmp)/'table.raw'
    for m,reps,mod in [(4,[1,9],0x13),(5,[1,6,7],0x25),
                       (6,[6,14,34,40],0x43),(7,[1],0x83)]:
        k=BinaryField(m,mod)
        for a in reps:
            f,g=bivariate_clv(k,a)
            for label,vals in [('F',f),('G',g)]:
                np.asarray(vals,dtype='<u4').tofile(path)
                run=subprocess.run([str(ROOT/'invn_canonical'),str(2*m),str(path),'1'],capture_output=True,text=True,check=True)
                expect=3 if m==4 and label=='F' else 2 if m==4 else 4
                assert 'APN=1' in run.stdout and f'maxsub={expect} ' in run.stdout
                results.append(dict(m=m,alpha=a,map=label,field_modulus=mod,result=run.stdout.strip()))
                print(label,m,a,'maxsub',expect,flush=True)
(ROOT/'invariant_audit_results.json').write_text(json.dumps(results,indent=2)+'\n')
