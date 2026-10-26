"""(prereg round36 7c 2차 3) R3′ 역장 차 d_ff 를 공통 계산 코드(모든 실행 쌍, C·B 저장)로 다시 낸다 — 단백질 사이 분포(크기 맞춤 몬테카를로의 재료).
자료 읽기·좌표·접힌 프레임·기준 평균 점은 r37_ff.py 와 같다(ATLAS stride 20, Q ≥ 0.75, 320 K 접힌 프레임 600 개의 평균 점).
창: 양쪽 전체(aSAMt 학습 창이 각 궤적 전체로 확인됐으므로 같은 창 규칙의 값). 실행은 N = 100 을 채우는 반복만(7c 와 같음).
출력: derived/r38/dff/<사슬>.json, results/r38_dff.json(사슬·도메인·길이·dmu·d2n·logspread 목록)"""
import os
for k_ in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k_]='1'
import sys, json, glob, time, itertools, zipfile, shutil, difflib, zlib
ROOT=sys.argv[1]; sys.argv=[sys.argv[0],ROOT,'__none__']
exec(open(ROOT+'/code/c1_core.py').read().split('doms = ONLY')[0])
src=open(ROOT+'/code/r37_ff.py').read(); exec(src[src.index('AA = dict('):src.index('def stats')])
cs=open(ROOT+'/code/r38_compare.py').read(); exec(cs[cs.index('def pdist_mean'):cs.index('def avg')])
import h5py, mdtraj as md
OUT=ROOT+'/derived/r38/dff'; os.makedirs(OUT,exist_ok=True); TMP='/content/dff_tmp'; os.makedirs(TMP,exist_ok=True); N=100
for x in [x for x in json.load(open(ROOT+'/results/overlap_mdcath_atlas.json'))['list_cover90'] if x.get('atlas_local')]:
    c,d=x['atlas'],x['mdcath']
    if os.path.exists(f'{OUT}/{c}.json'): continue
    try:
        t0=time.time(); rng=np.random.default_rng(zlib.crc32(f'r38_dff|{c}'.encode()) if 'zlib' in globals() else 38)
        h=h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5','r'); g=h[d]; smd,imd=md_atoms(g['pdbProteinAtoms'][()].decode())
        z=zipfile.ZipFile(f'{ROOT}/raw/atlas/{c}_protein.zip'); names=z.namelist(); pdbn=[n for n in names if n.endswith('.pdb')][0]; z.extract(pdbn,TMP); top=md.load(f'{TMP}/{pdbn}').topology; sat,iat=at_atoms(top)
        a0,b0,L=block(smd,sat)
        if L<0.9*min(len(smd),len(sat)): raise RuntimeError(f'짧은 일치 {L}/{len(smd)},{len(sat)}')
        imd,iat=imd[a0:a0+L],iat[b0:b0+L]
        MD=[g['320'][r]['coords'][:][:,imd.ravel()].reshape(-1,L,3,3) for r in sorted(g['320'].keys(),key=int)]; h.close()
        AT=[]
        for xt in sorted(n for n in names if n.endswith('.xtc')):
            z.extract(xt,TMP); t=md.load(f'{TMP}/{xt}',top=top,atom_indices=iat.ravel(),stride=20); os.remove(f'{TMP}/{xt}'); AT.append(t.xyz.reshape(t.n_frames,L,3,3)*10)
        ca0=MD[0][0,:,0]; ii,jj=np.triu_indices(L,4); dn=np.linalg.norm(ca0[ii]-ca0[jj],axis=-1); m_=dn<8.0; ii,jj,dn=ii[m_],jj[m_],dn[m_]
        Qf=lambda X:(np.linalg.norm(X[:,ii,0]-X[:,jj,0],axis=-1)<1.2*dn).mean(1)
        fMD=[X[Qf(X)>=QT] for X in MD]; fAT=[X[Qf(X)>=QT] for X in AT]
        allf=np.concatenate([X for X in fMD if len(X)])[:600]; ref=fs.mean_point(frames(allf))
        A_=lambda X: fs.align_to_reference(frames(X),ref).numpy().reshape(len(X),-1)
        pick=lambda rs:[A_(X[np.sort(rng.choice(len(X),N,replace=False))]) for X in rs if len(X)>=N]
        P=pick(fMD); Q=pick(fAT); ne=L-1
        res=dict(chain=c,dom=d,L=int(L),ne=int(ne),n_runs_md=len(P),n_runs_at=len(Q))
        if len(P)>=2 and len(Q)>=2: res['ff']=compare(P,Q,ne)
        res['sec']=round(time.time()-t0,1); json.dump(res,open(f'{OUT}/{c}.json','w')); print(time.strftime('%H:%M:%S'),c,L,res.get('ff',{}).get('dmu'),res['sec'],flush=True)
    except Exception as e:
        import traceback; print('ERR',c,repr(e)[:200],traceback.format_exc()[-300:],flush=True)
    for f in glob.glob(TMP+'/*'): os.remove(f) if os.path.isfile(f) else shutil.rmtree(f,ignore_errors=True)
J=[json.load(open(q)) for q in glob.glob(OUT+'/*.json')]
S=dict(written=time.strftime('%F %T'),n=len(J),n_ff=sum('ff' in j for j in J),rows=[dict(chain=j['chain'],dom=j['dom'],L=j['L'],**{k:round(j['ff'][k],6) for k in ('dmu','d2n','logspread','ratio','C','B_P','B_Q')}) for j in J if 'ff' in j])
P_=ROOT+'/results/r38_dff.json'; json.dump(S,open(P_,'w'),indent=1)
b=ROOT+'/results/_backup/'+time.strftime('%Y-%m-%d'); os.makedirs(b,exist_ok=True); shutil.copy(P_,b+'/'+time.strftime('%H%M%S_')+'r38_dff.json')
import numpy as _np; print('DFF_DONE', S['n_ff'], '중앙 dmu', _np.median([r['dmu'] for r in S['rows']]) if S['rows'] else None, flush=True)
