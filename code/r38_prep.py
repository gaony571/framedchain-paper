
"""(prereg round36 7c) 준비: 시험 90 도메인 서열, ATLAS train 과의 서열 동일성, AlphaFlow 입력 CSV, ColabFold MSA.
results/r38_prep.json, derived/r38/alphaflow_input.csv, derived/r38/msa/<dom>/a3m/<dom>.a3m"""
import os, sys, json, zipfile, time, subprocess, shutil, collections
for k_ in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k_]='1'
import pandas as pd, numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
R='/content/drive/MyDrive/atlas_framed_curves'
def backup(p):
    b=R+'/results/_backup/'+time.strftime('%Y-%m-%d'); os.makedirs(b,exist_ok=True); shutil.copy(p,b+'/'+time.strftime('%H%M%S_')+os.path.basename(p))
N=json.load(open(R+'/raw/asamt_test90_needs.json')); T=N['test90']
EXC=['1i2tA00','1nijA01','4qbuA03']
AA=dict(ALA='A',ARG='R',ASN='N',ASP='D',CYS='C',GLN='Q',GLU='E',GLY='G',HIS='H',HSD='H',HSE='H',HSP='H',HIE='H',HID='H',HIP='H',ILE='I',LEU='L',LYS='K',MET='M',PHE='F',PRO='P',SER='S',THR='T',TRP='W',TYR='Y',VAL='V')
z=zipfile.ZipFile(R+'/raw/models/asam/input_mdcath.zip')
seq={}
for t in T:
    txt=z.read(f'input_mdcath/{t}.320.pdb').decode(); s=[]; seen=set()
    for l in txt.splitlines():
        if l.startswith('ATOM') and l[12:16].strip()=='CA':
            k=(l[21],l[22:27]); 
            if k not in seen: seen.add(k); s.append(AA.get(l[17:21].strip(),'X'))
    seq[t]=''.join(s)
tr=pd.read_csv(R+'/code/ext/alphaflow/splits/atlas_train.csv')
al=Align.PairwiseAligner(); al.mode='local'; al.substitution_matrix=substitution_matrices.load('BLOSUM62'); al.open_gap_score=-10; al.extend_gap_score=-0.5
km=lambda s:set(s[i:i+3] for i in range(len(s)-2))
TK=[(n,s,km(s)) for n,s in zip(tr['name'],tr['seqres'])]
ident={}
for t in T:
    q=seq[t]; qk=km(q); cand=sorted(TK,key=lambda x:-len(qk&x[2]))[:25]; best=(0.0,None)
    for n,s,_ in cand:
        a=al.align(q,s)[0]; idn=sum(1 for (i0,i1),(j0,j1) in zip(*a.aligned) for i,j in zip(range(i0,i1),range(j0,j1)) if q[i]==s[j])
        best=max(best,(idn/len(q),n))
    ident[t]=dict(len=len(q),max_id=round(best[0],3),hit=best[1])
    print(t,ident[t],flush=True)
use=[t for t in T if t not in EXC]
os.makedirs(R+'/derived/r38',exist_ok=True)
pd.DataFrame({'name':use,'seqres':[seq[t] for t in use]}).to_csv(R+'/derived/r38/alphaflow_input.csv',index=False)
out=dict(n_test=len(T),excluded_pdbid=EXC,n_use=len(use),identity=ident,sensitivity_drop=[t for t in use if ident[t]['max_id']>=0.4],
         len_over_320=[t for t in T if ident[t]['len']>320],written=time.strftime('%F %T'))
P=R+'/results/r38_prep.json'; json.dump(out,open(P,'w'),indent=1); backup(P); print('PREP saved',flush=True)
msa=R+'/derived/r38/msa'
r=subprocess.run(f"cd {R}/code/ext/alphaflow && python -m scripts.mmseqs_query --split {R}/derived/r38/alphaflow_input.csv --outdir {msa}",shell=True)
got=[t for t in use if os.path.exists(f'{msa}/{t}/a3m/{t}.a3m')]
out['msa_done']=len(got); out['msa_rc']=r.returncode; json.dump(out,open(P,'w'),indent=1); backup(P)
print('MSA',len(got),'/',len(use),flush=True); open(R+'/results/logs/r38_prep.done','w').write(time.strftime('%F %T'))
