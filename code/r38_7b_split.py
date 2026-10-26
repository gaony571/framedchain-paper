"""(prereg round36 7b 판정용) 교환 가능한 10 ns 블록 반쪽 × 20 나눔으로 표준 지표 대 정답 순서 대조.
README 판 환경(/content/afenv)에서 돈다. 원 분석 코드(AlphaFlow 0408d7c scripts/analyze_ensembles.py)의 함수 정의를 그대로 가져와
main() 과 같은 순서로 계산하되 SASA(노출 잔기 J·MI)는 뺀다(prereg 7b "계산 범위": 시험 예측 1–4 는 SASA 를 쓰지 않음).
기준 쪽 중첩·기준 PCA 는 나눔마다 한 번, 후보 프레임 수에 달린 무작위 색인(RAND1/RAND2/RAND1K)은 원 코드처럼 대조마다 seed 137 로 다시 뽑는다.
나눔: 반복(R1–R3)마다 1001 프레임을 100 프레임(10 ns) 블록 10 개로 자르고(끝 1 프레임 버림), 반복마다 5 개를 무작위로 A(기준)에, 나머지를 B 에.
후보 = B 에서 250 프레임 무작위 → 대조 8: true, shrink0.8, shrink0.6, inflate1.2, trim20, trim40, scramble, xtalnoise (r38_controls.py 의 make 와 같은 정의).
에너지 거리(공통 계산 코드): 맞춘 좌표(CA, C, O 프레임, Å), 기준 평균 점 = A 프레임(최대 600)의 평균 점. 기준 실행 = 반복별 A 블록 프레임 100 개씩(3),
후보 실행 = 기준과 대칭으로 원래 프레임이 나온 반복별 묶음(true·줄임·부풂·잘라내기; 2026-10-19 표적 1·나눔 1 시험 뒤 고침 — 무작위 반쪽이면 이웃 프레임이 닮아 후보 쪽 바닥이 작아져 A·B 교환 가능성이 추정량에서 깨짐), 프레임 출처가 없는 섞음·결정+잡음은 무작위 셋. 교차 C 도 서로 다른 반복끼리의 쌍만 쓴다(A·B 가 같은 반복에서 나오므로 같은 반복 쌍은 독립 실행이 아님 — 세 항 모두 독립 실행 쌍만; 진짜 후보의 기댓값 0). C·B_P·B_Q 저장, Δ² = 2C − B_P − B_Q, d2n = Δ²/√변.
사용: /content/afenv/bin/python r38_7b_split.py <ROOT> <작업자 수> [나눔 수=20] [표적 쉼표 목록]
출력: derived/r38/b7/<표적>__s<k>.json (이어 돌리기)"""
import os, sys
for k_ in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k_]='1'
import json, time, zlib, itertools
import numpy as np, mdtraj
ROOT=sys.argv[1]; NW=int(sys.argv[2]); NS=int(sys.argv[3]) if len(sys.argv)>3 else 20
AF=ROOT+'/code/ext/alphaflow'; AD='/content/r38_bar/atlas'; OUT=ROOT+'/derived/r38/b7'; os.makedirs(OUT,exist_ok=True)
src=open(AF+'/scripts/analyze_ensembles.py').read()
from sklearn.decomposition import PCA
from scipy.optimize import linear_sum_assignment
exec(src[src.index('def get_pca'):src.index('def condense_sidechain_sasas')])      # get_pca, get_rmsds
exec(src[src.index('def get_mean_covar'):src.index('def main(name)')])             # get_mean_covar, sqrtm, get_wasserstein, align_tops (그리고 sasa_mi 등)
sys.path[:0]=[ROOT+'/framespace_pkg',ROOT+'/framespace_wip']
import torch, framespace as fs; torch.set_num_threads(1)
CT=['true','shrink0.8','shrink0.6','inflate1.2','trim20','trim40','scramble','xtalnoise']
def frames_t(X): return fs.coordinates_from_backbone(torch.as_tensor(X[:,:,0]),torch.as_tensor(X[:,:,1]),torch.as_tensor(X[:,:,2]))
def cost(A,B): return np.maximum((A*A).sum(1)[:,None]+(B*B).sum(1)[None,:]-2*A@B.T,0)
def pmean(A,B=None):
    if B is None:
        D=np.sqrt(cost(A,A)); n=len(A); return float((D.sum()-np.trace(D))/(n*(n-1)))
    return float(np.sqrt(cost(A,B)).mean())
def edist_rep(Pd,Qd,ne):
    """Pd, Qd: {반복 번호: 맞춘 좌표}. 교차·바닥 모두 서로 다른 반복끼리의 쌍만."""
    C=float(np.mean([pmean(Pd[a],Qd[b]) for a in Pd for b in Qd if a!=b])); BP=float(np.mean([pmean(Pd[a],Pd[b]) for a,b in itertools.combinations(sorted(Pd),2)])); BQ=float(np.mean([pmean(Qd[a],Qd[b]) for a,b in itertools.combinations(sorted(Qd),2)]))
    D2=2*C-BP-BQ; w=(BP+BQ)/2
    return dict(C=C,B_P=BP,B_Q=BQ,D2=D2,d2n=D2/np.sqrt(ne),ratio=float(np.sqrt(max((C/w)**2-1,0))))
def edist_runs(P,Q,ne):
    C=float(np.mean([pmean(a,b) for a in P for b in Q])); BP=float(np.mean([pmean(a,b) for a,b in itertools.combinations(P,2)])); BQ=float(np.mean([pmean(a,b) for a,b in itertools.combinations(Q,2)]))
    D2=2*C-BP-BQ; w=(BP+BQ)/2
    return dict(C=C,B_P=BP,B_Q=BQ,D2=D2,d2n=D2/np.sqrt(ne),ratio=float(np.sqrt(max((C/w)**2-1,0))))
def bb_idx(top):
    keep=[r for r in top.residues if r.is_protein and all(any(a.name==n for a in r.atoms) for n in ('CA','C','O'))]
    return np.array([[next(a.index for a in r.atoms if a.name==n) for n in ('CA','C','O')] for r in keep])

def metrics(out):
    """print_analysis.py 의 정의(SASA 제외)."""
    m={'md_pairwise':float(out['ref_mean_pairwise_rmsd']),'af_pairwise':float(out['af_mean_pairwise_rmsd']),'cosine_sim':float(abs(out['cosine_sim'])),
       'emd_mean':float(np.square(out['emd_mean']).mean()**0.5),'emd_var':float(np.square(np.real(out['emd_var'])).mean()**0.5),
       'rmsf_pearson':float(np.corrcoef(out['af_rmsf'],out['ref_rmsf'])[0,1]),
       'refemd':float(out['EMD,ref']['ref|af']),'jointemd':float(out['EMD,joint']['ref|af']),'bures_fallback':int(out['bures_fallback'])}
    m['rmwd']=float(np.sqrt(m['emd_mean']**2+m['emd_var']**2))
    cm=out['crystal_distmat']<0.8
    rt=(~cm)&(out['ref_contact_prob']>0.1); at=(~cm)&(out['af_contact_prob']>0.1); rw=cm&(out['ref_contact_prob']<0.9); aw=cm&(out['af_contact_prob']<0.9)
    m['weak_contacts_iou']=float((rw&aw).sum()/(rw|aw).sum()) if (rw|aw).sum() else float('nan')
    m['transient_contacts_iou']=float((rt&at).sum()/(rt|at).sum()) if (rt|at).sum() else float('nan')
    return m

def run(task):
    name,s=task; f=f'{OUT}/{name}__s{s}.json'
    if os.path.exists(f): return name,s,'skip'
    t0=time.time(); rng=np.random.default_rng(zlib.crc32(f'7b|{name}|{s}'.encode()))
    pdb=f'{AD}/{name}/{name}.pdb'; xt=mdtraj.load(pdb); top=xt.topology
    reps=[mdtraj.load(f'{AD}/{name}/{name}_prod_R{k}_fit.xtc',top=pdb) for k in (1,2,3)]
    A_idx=[]; B_idx=[]
    for t in reps:
        blocks=[np.arange(b*100,(b+1)*100) for b in range(t.n_frames//100)]; pick=rng.permutation(len(blocks))
        A_idx.append(np.sort(np.concatenate([blocks[i] for i in pick[:len(blocks)//2]]))); B_idx.append(np.sort(np.concatenate([blocks[i] for i in pick[len(blocks)//2:]])))
    refA=[t[a] for t,a in zip(reps,A_idx)]; poolB=mdtraj.join([t[b] for t,b in zip(reps,B_idx)])
    labB=np.concatenate([np.full(len(b),k) for k,b in enumerate(B_idx)])
    # 후보 250 프레임과 대조(r38_controls.py make 와 같은 정의)
    heavy=[a.index for a in top.atoms if a.element.symbol!='H']
    sel=np.sort(rng.choice(poolB.n_frames,min(250,poolB.n_frames),replace=False)); cand=poolB[sel]; lab=labB[sel]; cand.superpose(xt,atom_indices=heavy)
    X=cand.xyz; mu=X.mean(0); rm=np.sqrt(((X-mu)**2).sum(-1).mean(1)); per_atom=np.sqrt(((X-mu)**2).sum(-1).mean(0))
    C_={}; LAB={}
    for c in CT:
        if c=='true': Y=X
        elif c.startswith('shrink') or c.startswith('inflate'): lam=float(c[6:] if c.startswith('shrink') else c[7:]); Y=mu+lam*(X-mu)
        elif c.startswith('trim'): q=int(c[4:])/100; Y=X[np.argsort(rm)[:int(round(len(X)*(1-q)))]]
        elif c=='scramble':
            Y=np.empty_like(X)
            for r in top.residues:
                ai=[a.index for a in r.atoms]; p=rng.permutation(len(X)); Y[:,ai]=mu[ai]+(X[p][:,ai]-mu[ai])
        elif c=='xtalnoise': Y=xt.xyz[0][None]+rng.normal(size=X.shape)*(per_atom/np.sqrt(3))[None,:,None]
        C_[c]=Y.astype(np.float32)
        LAB[c]=lab if c in ('true','shrink0.8','shrink0.6','inflate1.2') else (lab[np.argsort(rm)[:len(Y)]] if c.startswith('trim') else rng.integers(0,3,len(Y)))
    # ---- 원 코드 main() 의 기준 쪽(대조와 무관한 부분) ----
    traj_aa=mdtraj.join(refA); ref_aa=mdtraj.load(pdb)
    traj_aa.atom_slice([a.index for a in traj_aa.top.atoms if a.element.symbol!='H'],True)
    ref_aa.atom_slice([a.index for a in ref_aa.top.atoms if a.element.symbol!='H'],True)
    proto=mdtraj.Trajectory(C_['true'],top); proto.atom_slice([a.index for a in proto.top.atoms if a.element.symbol!='H'],True)
    refmask,afmask=align_tops(traj_aa.top,proto.top)
    traj_aa.atom_slice(refmask,True); ref_aa.atom_slice(refmask,True)
    traj_aa.superpose(ref_aa)
    ca_mask=[a.index for a in traj_aa.top.atoms if a.name=='CA']
    traj=traj_aa.atom_slice(ca_mask,False); ref=ref_aa.atom_slice(ca_mask,False); traj.superpose(ref)
    n_atoms=traj.n_atoms
    ref_pca,ref_coords=get_pca(traj.xyz); seed_ref_pca=ref_pca.transform(ref.xyz.reshape(1,-1))
    ref_rmsf=mdtraj.rmsf(traj_aa,ref_aa)*10
    crystal_distmat=np.linalg.norm(ref.xyz[0,None,:]-ref.xyz[0,:,None],axis=-1)
    # 에너지 거리 준비(맞춘 좌표)
    bi=bb_idx(top); L=len(bi); ne=L-1
    RX=[t.xyz[:,bi.ravel()].reshape(t.n_frames,L,3,3)*10 for t in refA]
    refpt=fs.mean_point(frames_t(np.concatenate([r[::max(1,len(r)//200)] for r in RX])[:600]))
    A_=lambda Y: fs.align_to_reference(frames_t(Y),refpt).numpy().reshape(len(Y),-1)
    Pd={k:A_(r[np.sort(rng.choice(len(r),min(100,len(r)),replace=False))]) for k,r in enumerate(RX)}
    res={'name':name,'split':s,'n_ref':int(traj_aa.n_frames),'L':int(L),'controls':{}}
    from types import SimpleNamespace
    def get_emd(ref_coords1, ref_coords2, af_coords, seed_coords, K=None):
        exec_ns={}
    for c in CT:
        Y=C_[c]; aft=mdtraj.Trajectory(Y,top); aft.atom_slice([a.index for a in aft.top.atoms if a.element.symbol!='H'],True); aft.atom_slice(afmask,True)
        np.random.seed(137)
        RAND1=np.random.randint(0,traj_aa.n_frames,aft.n_frames); RAND2=np.random.randint(0,traj_aa.n_frames,aft.n_frames); RAND1K=np.random.randint(0,traj_aa.n_frames,1000)
        aft.superpose(ref_aa); aftraj=aft.atom_slice(ca_mask,False); aftraj.superpose(ref)
        out={}
        af_coords_ref_pca=ref_pca.transform(aftraj.xyz.reshape(aftraj.n_frames,-1))
        af_pca,af_coords=get_pca(aftraj.xyz)
        joint_pca,_=get_pca(np.concatenate([traj[RAND1].xyz,aftraj.xyz]))
        af_coords_joint_pca=joint_pca.transform(aftraj.xyz.reshape(aftraj.n_frames,-1)); ref_coords_joint_pca=joint_pca.transform(traj.xyz.reshape(traj.n_frames,-1)); seed_joint=joint_pca.transform(ref.xyz.reshape(1,-1))
        out['af_rmsf']=mdtraj.rmsf(aft,ref_aa)*10; out['ref_rmsf']=ref_rmsf
        ref_mean,ref_covar=get_mean_covar(traj_aa[RAND1K].xyz); af_mean,af_covar=get_mean_covar(aft.xyz)
        out['emd_mean']=(np.square(ref_mean-af_mean).sum(-1)**0.5)*10
        try:
            out['emd_var']=(np.trace(ref_covar+af_covar-2*sqrtm(ref_covar@af_covar),axis1=1,axis2=2)**0.5)*10; out['bures_fallback']=0
        except:
            out['emd_var']=np.trace(ref_covar)**0.5*10; out['bures_fallback']=1
        ref_distmat=np.linalg.norm(traj[RAND1].xyz[:,None,:]-traj[RAND1].xyz[:,:,None],axis=-1); af_distmat=np.linalg.norm(aftraj.xyz[:,None,:]-aftraj.xyz[:,:,None],axis=-1)
        out['ref_contact_prob']=(ref_distmat<0.8).mean(0); out['af_contact_prob']=(af_distmat<0.8).mean(0); out['crystal_distmat']=crystal_distmat
        out['ref_mean_pairwise_rmsd']=get_rmsds(traj[RAND1].xyz,traj[RAND2].xyz,broadcast=True).mean(); out['af_mean_pairwise_rmsd']=get_rmsds(aftraj.xyz,aftraj.xyz,broadcast=True).mean()
        out['cosine_sim']=(ref_pca.components_[0]*af_pca.components_[0]).sum()
        def emd_ref_af(rc1, afc, K=2):
            rc1=rc1.reshape(rc1.shape[0],-1)[:,:K]; afc=afc.reshape(afc.shape[0],-1)[:,:K]
            dm=(np.square(rc1[:,None]-afc[None]).sum(-1))**0.5/n_atoms**0.5*10
            return {'ref|af':get_wasserstein(dm)}
        out['EMD,ref']=emd_ref_af(ref_coords[RAND1],af_coords_ref_pca); out['EMD,joint']=emd_ref_af(ref_coords_joint_pca[RAND1],af_coords_joint_pca)
        m=metrics(out)
        # 에너지 거리: 후보 실행 = 무작위 반쪽 둘
        Yb=Y[:,bi.ravel()].reshape(len(Y),L,3,3)*10; Ay=A_(Yb)
        Qd={k:Ay[LAB[c]==k] for k in range(3) if (LAB[c]==k).sum()>=10}
        m.update({'E_'+k:v for k,v in edist_rep(Pd,Qd,ne).items()}); m['n_cand']=int(len(Y)); m['n_cand_runs']=[int((LAB[c]==k).sum()) for k in range(3)]
        res['controls'][c]=m
    res['sec']=round(time.time()-t0,1); json.dump(res,open(f,'w'))
    return name,s,res['sec']

if __name__=='__main__':
    names=sys.argv[4].split(',') if len(sys.argv)>4 else sorted(os.listdir(AD))
    tasks=[(n,s) for s in range(NS) for n in names]
    from multiprocessing import Pool
    with Pool(NW) as p:
        for name,s,sec in p.imap_unordered(run,tasks):
            print(time.strftime('%H:%M:%S'),name,'s%d'%s,sec,flush=True)
    print('B7_DONE',flush=True)
