"""(prereg round36 7c, 공통 계산 코드) 모델·MD 앙상블 비교. 고친 통계: 모든 실행 쌍, C·B_P·B_Q·실행 안 평균 거리 저장,
Δ² = 2C − B_P − B_Q, 비 = √((C/w̄)² − 1), Δμ² = 교차 − (바닥_P + 바닥_Q)/2 (실행 평균), 로그 퍼짐 비(r37 정의).
좌표: r37_ff 와 같다 — 잔기마다 (CA, C, O), Å, 서열은 block() 으로 맞춤, 320 K 기준 평균 점(derived/r38/md/<dom>.npz 의 ref)에 맞춘 좌표.
접힌 프레임: 첫 320 K 프레임의 천연 접촉(npz 의 ii, jj, dn), Q ≥ q (q = 0.75 주, 0.7·0.8 민감도). 'all' = 거르지 않음.
MD 실행 = 반복(창 안 접힌 프레임에서 N = 100 무작위, 못 채우면 그 반복 제외), 칸은 실행 둘 이상. 창: full(주; aSAMt 학습 창), half(뒤 절반, 시간 기준).
모델 실행 = 표본 무작위 반쪽(20 번 나눠 평균). 반쪽마다 ≥ 40.
칸: MD T 대 MD 320 (온도 반응), aSAMt T 대 MD 320, aSAMt T 대 MD T, AlphaFlow 대 MD 320, MD 320 반복 나눔(2 대 3, 10 가지; 귀무).
사용: python r38_compare.py <ROOT> <도메인 목록 파일 또는 쉼표 목록>
출력: derived/r38/compare/<dom>.json (수치), <dom>.npz (실행 평균·변별 퍼짐; E5 용)"""
import os
for k_ in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k_]='1'
import sys, json, glob, time, itertools, zlib, shutil
ROOT=sys.argv[1]; arg=sys.argv[2]; sys.argv=[sys.argv[0],ROOT,'__none__']
exec(open(ROOT+'/code/c1_core.py').read().split('doms = ONLY')[0])
src=open(ROOT+'/code/r37_ff.py').read(); exec(src[src.index('AA = dict('):src.index('def stats')])   # AA, TT, md_atoms, at_atoms, block
import h5py, mdtraj as md, difflib   # difflib: r37_ff 의 block() 이 씀(그 파일 맨 위에서 가져옴)
doms=[x.strip() for x in (open(arg).read().split() if os.path.exists(arg) else arg.split(',')) if x.strip()]
OUT=ROOT+'/derived/r38/compare'; os.makedirs(OUT,exist_ok=True)
N=100; MINH=40; NSPLIT=20; QS=(0.75,0.7,0.8)

def pdist_mean(A,B=None):
    """B 없으면 A 안 서로 다른 프레임 쌍의 평균 거리, 있으면 A–B 모든 쌍의 평균 거리."""
    if B is None:
        D=np.sqrt(cost(A,A)); n=len(A); return float((D.sum()-np.trace(D))/(n*(n-1)))
    return float(np.sqrt(cost(A,B)).mean())
def compare(P,Q,ne):
    """P, Q: 실행 목록(각 (n, 9ne) 맞춘 좌표). 모든 실행 쌍. 반환: 수치 dict."""
    C=float(np.mean([pdist_mean(a,b) for a in P for b in Q]))
    BP=float(np.mean([pdist_mean(a,b) for a,b in itertools.combinations(P,2)])); BQ=float(np.mean([pdist_mean(a,b) for a,b in itertools.combinations(Q,2)]))
    aP=float(np.mean([pdist_mean(a) for a in P])); aQ=float(np.mean([pdist_mean(a) for a in Q]))
    mP=[a.mean(0) for a in P]; mQ=[b.mean(0) for b in Q]; m2=lambda x,y: float(((x-y)**2).sum())
    cm=np.mean([m2(x,y) for x in mP for y in mQ]); fP=np.mean([m2(x,y) for x,y in itertools.combinations(mP,2)]); fQ=np.mean([m2(x,y) for x,y in itertools.combinations(mQ,2)])
    sp=lambda rs: np.mean([((r-r.mean(0))**2).sum(1).mean() for r in rs])
    D2=2*C-BP-BQ; w=(BP+BQ)/2
    return dict(C=C,B_P=BP,B_Q=BQ,a_P=aP,a_Q=aQ,D2=D2,d2n=D2/np.sqrt(ne),wbar=w,ratio=float(np.sqrt(max((C/w)**2-1,0))),
                dmu=float((cm-(fP+fQ)/2)/ne),Dmu2=float(cm-(fP+fQ)/2),logspread=float(np.log(sp(Q)/sp(P))),nP=[len(a) for a in P],nQ=[len(b) for b in Q])
def avg(ds):
    out={}
    for k in ds[0]:
        out[k]=ds[0][k] if isinstance(ds[0][k],list) else float(np.mean([d[k] for d in ds]))
    return out
def edge_profile(runs,ne):
    """실행별 평균(9ne)과 변별 퍼짐(ne; 변의 9 성분 제곱 편차 합의 평균)."""
    return np.stack([r.mean(0) for r in runs]), np.stack([((r-r.mean(0))**2).reshape(len(r),ne,9).sum(2).mean(0) for r in runs])

for d in doms:
    if os.path.exists(f'{OUT}/{d}.json'): continue
    t0=time.time(); rng=np.random.default_rng(zlib.crc32(f'r38_compare|{d}'.encode()))
    try:
        z=np.load(f'{ROOT}/derived/r38/md/{d}.npz'); ref=torch.as_tensor(z['ref']); ii,jj,dn=z['ii'],z['jj'],z['dn']
        h=h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5','r'); g=h[d]; smd,imd=md_atoms(g['pdbProteinAtoms'][()].decode()); L=len(smd)
        MD={T:[g[T][r]['coords'][:][:,imd.ravel()].reshape(-1,L,3,3) for r in sorted(g[T].keys(),key=int)] for T in TT if T in g}; h.close()
        assert np.array_equal(imd, z['imd']), 'imd 다름'
        ne=L-1; Qf=lambda X:(np.linalg.norm(X[:,ii,0]-X[:,jj,0],axis=-1)<1.2*dn).mean(1)
        A_=lambda X: fs.align_to_reference(frames(X),ref).numpy().reshape(len(X),-1)
        res=dict(dom=d,L=int(L),ne=int(ne),seq=smd,check={},cells={},nruns={}); prof={}
        # MD: 창·Q 마다 실행(반복)
        mdrun={}
        for T in MD:
            for W_ in ('full','half'):
                reps=[X if W_=='full' else X[len(X)//2:] for X in MD[T]]
                Qv=[Qf(X) for X in reps]; Ar=[A_(X) for X in reps]
                for q in QS+('all',):
                    runs=[]
                    for X,qq in zip(Ar,Qv):
                        idx=np.arange(len(X)) if q=='all' else np.where(qq>=q)[0]
                        if len(idx)>=N: runs.append(X[np.sort(rng.choice(idx,N,replace=False))])
                    mdrun[(T,W_,q)]=runs; res['nruns'][f'MD{T}|{W_}|{q}']=len(runs)
        def model_runs(Xm,q):
            Am=A_(Xm); qq=Qf(Xm); idx=np.arange(len(Xm)) if q=='all' else np.where(qq>=q)[0]
            return Am[idx], float((qq>=0.75).mean())
        def model_cells(tag,Am,ref_key):
            P=mdrun[ref_key]
            if len(P)<2 or len(Am)<2*MINH: return None
            outs=[];
            for s in range(NSPLIT):
                perm=rng.permutation(len(Am)); hf=len(Am)//2; Qr=[Am[perm[:hf]],Am[perm[hf:2*hf]]]; outs.append(compare(P,Qr,ne))
            return avg(outs)
        # MD 온도 반응과 320 귀무
        for W_ in ('full','half'):
            for q in QS+('all',):
                P=mdrun[('320',W_,q)]
                if len(P)>=2:
                    mm,ss=edge_profile(P,ne); prof[f'MD320|{W_}|{q}|mean']=mm; prof[f'MD320|{W_}|{q}|spr']=ss
                for T in TT[1:]:
                    Qr=mdrun.get((T,W_,q),[])
                    if len(P)>=2 and len(Qr)>=2:
                        res['cells'][f'MD{T}-MD320|{W_}|{q}']=compare(P,Qr,ne)
                        mm,ss=edge_profile(Qr,ne); prof[f'MD{T}|{W_}|{q}|mean']=mm; prof[f'MD{T}|{W_}|{q}|spr']=ss
                if len(P)==5:
                    nl=[compare([P[i] for i in a],[P[i] for i in range(5) if i not in a],ne) for a in itertools.combinations(range(5),2)]
                    res['cells'][f'null320|{W_}|{q}']=[x['Dmu2'] for x in nl]; res['cells'][f'null320d2n|{W_}|{q}']=[x['d2n'] for x in nl]   # 2차 1 의 귀무는 Δ²/√변
        # 모델
        files={}
        for T in TT:
            f=f'{ROOT}/derived/r38/asamt/{d}.{T}'
            if os.path.exists(f+'.traj.dcd'): files[f'aSAMt{T}']=(f+'.traj.dcd',f+'.top.pdb',T)
        f=f'{ROOT}/derived/r38/alphaflow/{d}.pdb'
        if os.path.exists(f): files['AlphaFlow']=(f,None,'320')
        for name,(fx,ftop,T) in files.items():
            t=md.load(fx,top=ftop) if ftop else md.load(fx)
            sm,im=at_atoms(t.topology); a0,b0,Lb=block(smd,sm)
            res['check'][name]=dict(L_model=len(sm),L_md=L,block=int(Lb),same_seq=bool(sm==smd),n_frames=int(t.n_frames))
            if Lb!=L or sm!=smd: res['check'][name]['note']='서열 또는 길이 다름 — 칸을 내지 않음'; continue
            Xm=t.xyz[:,im.ravel()].reshape(t.n_frames,L,3,3)*10
            for W_ in ('full','half'):
                for q in QS+('all',):
                    Am,ffrac=model_runs(Xm,q); res['check'][name]['folded_frac_q075']=ffrac
                    c1=model_cells(name,Am,('320',W_,q))
                    if c1: res['cells'][f'{name}-MD320|{W_}|{q}']=c1
                    if name.startswith('aSAMt') and T!='320':
                        c2=model_cells(name,Am,(T,W_,q))
                        if c2: res['cells'][f'{name}-MD{T}|{W_}|{q}']=c2
                    if W_=='full' and len(Am)>=2*MINH:
                        hf=len(Am)//2; perm=rng.permutation(len(Am)); mm,ss=edge_profile([Am[perm[:hf]],Am[perm[hf:2*hf]]],ne)
                        prof[f'{name}|{q}|mean']=mm; prof[f'{name}|{q}|spr']=ss
        res['sec']=round(time.time()-t0,1)
        np.savez_compressed(f'{OUT}/{d}.npz',**prof); json.dump(res,open(f'{OUT}/{d}.json','w'))
        c=res['cells']; k='full|0.75'
        print(time.strftime('%H:%M:%S'),d,L,'check',res['check'],'| aSAMt320 dmu',round(c.get(f'aSAMt320-MD320|{k}',{}).get('dmu',np.nan),5),
              'MD379',round(c.get(f'MD379-MD320|{k}',{}).get('dmu',np.nan),5),'AF',round(c.get(f'AlphaFlow-MD320|{k}',{}).get('dmu',np.nan),5),res['sec'],'s',flush=True)
    except Exception as e:
        import traceback; print('ERR',d,repr(e)[:300],traceback.format_exc()[-800:],flush=True)
print('COMPARE_DONE',flush=True)
