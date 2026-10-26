
"""(prereg round36 7c, 90 판 셋) mdCATH 쪽 준비 — r38_md_ref.py 사본: 시험 87 도메인의 기준 평균 점(320 K 접힌 프레임), 천연 접촉, 원자 색인,
같은 도메인 온도 곡선(320 대 348/379/413/450 K; 전체·뒤 절반; r37_ff 와 같은 통계). derived/r38/md/<dom>.{npz,json}, results/r38_md_ref.json"""
import os
for k_ in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k_]='1'
import sys, json, glob, time, itertools, shutil
ROOT=sys.argv[1]; sys.argv=[sys.argv[0],ROOT,'__none__']
exec(open(ROOT+'/code/c1_core.py').read().split('doms = ONLY')[0])
src=open(ROOT+'/code/r37_ff.py').read()
exec(src[src.index('AA = dict('):src.index('for x in [x for x')])   # AA, TT, md_atoms, at_atoms, block, stats (r37_ff 그대로)
import h5py
def backup(p):
    b=ROOT+'/results/_backup/'+time.strftime('%Y-%m-%d'); os.makedirs(b,exist_ok=True); shutil.copy(p,b+'/'+time.strftime('%H%M%S_')+os.path.basename(p))
OUT=ROOT+'/derived/r38/md'; os.makedirs(OUT,exist_ok=True)
prep=json.load(open(ROOT+'/results/r38_prep.json')); EXC=set(prep['excluded_pdbid'])
doms=['1i2tA00','1nijA01','4qbuA03']   # 90 판(2026-10-19): 처음 준비에서 PDB ID 겹침으로 뺀 셋. 기준 정의는 r38_md_ref.py 그대로
for d in doms:
    if os.path.exists(f'{OUT}/{d}.json'): continue
    try:
        t0=time.time(); rng=np.random.default_rng(38)
        h=h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5','r'); g=h[d]; smd,imd=md_atoms(g['pdbProteinAtoms'][()].decode()); L=len(smd)
        MD={T:[g[T][r]['coords'][:][:,imd.ravel()].reshape(-1,L,3,3) for r in sorted(g[T].keys(),key=int)] for T in TT if T in g}; h.close()
        ca0=MD['320'][0][0,:,0]; ii,jj=np.triu_indices(L,4); dn=np.linalg.norm(ca0[ii]-ca0[jj],axis=-1); m_=dn<8.0; ii,jj,dn=ii[m_],jj[m_],dn[m_]
        Qf=lambda X:(np.linalg.norm(X[:,ii,0]-X[:,jj,0],axis=-1)<1.2*dn).mean(1)
        foldMD={T:[X[Qf(X)>=QT] for X in MD[T]] for T in MD}
        allf=np.concatenate([X for X in foldMD['320'] if len(X)])[:600]; ref=fs.mean_point(frames(allf))
        A_=lambda X: fs.align_to_reference(frames(X),ref).numpy().reshape(len(X),-1)
        ok=lambda rs:[A_(X) for X in rs if len(X)>=40]
        R320=ok(foldMD['320']); ne=L-1
        res=dict(dom=d,L=int(L),ne=int(ne),seq=smd,nfold={T:[int(len(X)) for X in foldMD[T]] for T in foldMD},nall={T:[int(len(X)) for X in MD[T]] for T in MD},temp={},temp_half={})
        H320=[X[len(X)//2:] for X in R320 if len(X)//2>=40]
        for T in TT[1:]:
            if T not in foldMD: continue
            RT=ok(foldMD[T])
            if len(R320)>=2 and len(RT)>=2: res['temp'][T]=stats(R320,RT,ne,rng)
            HT=[X[len(X)//2:] for X in RT if len(X)//2>=40]
            if len(H320)>=2 and len(HT)>=2: res['temp_half'][T]=stats(H320,HT,ne,rng)
        np.savez(f'{OUT}/{d}.npz',ref=np.asarray(ref),ii=ii,jj=jj,dn=dn,imd=imd)
        res['sec']=round(time.time()-t0,1); json.dump(res,open(f'{OUT}/{d}.json','w')); print(time.strftime('%H:%M:%S'),d,L,res['temp_half'].get('450'),res['sec'],flush=True)
    except Exception as e:
        import traceback; print('ERR',d,repr(e)[:200],traceback.format_exc()[-300:],flush=True)
