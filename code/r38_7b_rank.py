"""(prereg round36 7b 모델 순위, 서술) 공개 표본 네 세트의 표준 지표 중앙값(원 분석 코드, README 판 환경)과 에너지 거리 Δ² 중앙값을 견준다.
에너지 거리: 공통 계산 코드(모든 실행 쌍, C·B 저장; Δ² = 2C − B_P − B_Q, d2n = Δ²/√변, 비, Δμ²/변, 로그 퍼짐 비).
  맞춘 좌표 = 잔기마다 (CA, C, O) 프레임 좌표(Å), 기준 평균 점 = ATLAS MD 세 반복에서 고르게 고른 프레임(최대 600)의 평균 점.
  MD 실행 = 반복(R1–R3)마다 100 프레임 무작위 — 표적마다 한 번 뽑아 네 세트에 같게 쓴다. 모델 실행 = 250 표본의 무작위 반쪽 둘, 20 나눔 평균(r38_compare 와 같은 구성).
  서열은 block() 으로 맞추고 맞춘 구간만 쓴다(맞춘 길이 < 90% 면 그 표적을 뺀다). 원자는 명시 색인으로 고르므로 mdtraj 판에 무관(결과 3.6b 의 atom_slice 문제 없음).
표준 지표: results/r38_af_bar_readmeenv.txt (AlphaFlow-MD, 표 재현) 와 results/r38_af_bar_readmeenv_<세트>.txt (q_7b_rank.sh) 의 print_analysis 중앙값.
  가까움: 상관·비율·Jaccard·MI ρ 는 클수록, RMWD·PCA W₂ 는 작을수록, 쌍 RMSD·RMSF 는 MD 값과의 차가 작을수록.
순위: 세트마다 에너지 거리 d2n 중앙값(표적 부트스트랩 2000, seed 20261026, 95% 구간). 지표마다 네 세트 순위와 에너지 거리 순위의 Kendall τ.
순위표 자체를 주장하지 않는다(prereg 7b). 순위가 갈리는 세트 쌍은 판정용 대조(결과 3.6c)로 읽는다.
사용: python r38_7b_rank.py <ROOT>   (Colab 기본 python)
출력: derived/r38/b7_rank/<세트>/<표적>.json, results/r38_7b_rank.json, results/r38_7b_rank.txt (사본 results/_backup/<날짜>/)"""
import os, sys, json, glob, time, zlib, itertools, zipfile, shutil, re, difflib
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '2'
import numpy as np, torch, mdtraj as md
from scipy.stats import kendalltau
ROOT = sys.argv[1]; sys.path[:0] = [ROOT + '/framespace_pkg']; import framespace as fs; torch.set_num_threads(2)
src = open(ROOT + '/code/r37_ff.py').read(); exec(src[src.index('AA = dict('):src.index('def stats')])   # at_atoms, block
SETS = ['alphaflow_md_base_202402', 'alphaflow_md_templates_base_202402', 'esmflow_md_base_202402', 'esmflow_md_templates_base_202402']
TXT = {SETS[0]: ROOT + '/results/r38_af_bar_readmeenv.txt', **{s: ROOT + f'/results/r38_af_bar_readmeenv_{s}.txt' for s in SETS[1:]}}
AD = '/content/r38_bar/atlas'; SD = '/content/r38_rank_samples'; OUT = ROOT + '/derived/r38/b7_rank'; RES = ROOT + '/results'
N = 100; NSPLIT = 20


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)
def backup(p):
    b = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + os.path.basename(p))
def frames_t(X): return fs.coordinates_from_backbone(torch.as_tensor(X[:, :, 0]), torch.as_tensor(X[:, :, 1]), torch.as_tensor(X[:, :, 2]))
def cost(A, B): return np.maximum((A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2 * A @ B.T, 0)
def pmean(A, B): return float(np.sqrt(cost(A, B)).mean())
def edist_runs(P, Q, ne):
    C = np.mean([pmean(a, b) for a in P for b in Q]); BP = np.mean([pmean(a, b) for a, b in itertools.combinations(P, 2)]); BQ = np.mean([pmean(a, b) for a, b in itertools.combinations(Q, 2)])
    mP = [a.mean(0) for a in P]; mQ = [b.mean(0) for b in Q]; m2 = lambda x, y: float(((x - y) ** 2).sum())
    cm = np.mean([m2(x, y) for x in mP for y in mQ]); fP = np.mean([m2(x, y) for x, y in itertools.combinations(mP, 2)]); fQ = np.mean([m2(x, y) for x, y in itertools.combinations(mQ, 2)])
    sp = lambda rs: np.mean([((r - r.mean(0)) ** 2).sum(1).mean() for r in rs])
    D2 = 2 * C - BP - BQ; w = (BP + BQ) / 2
    return dict(C=float(C), B_P=float(BP), B_Q=float(BQ), D2=float(D2), d2n=float(D2 / np.sqrt(ne)), wbar=float(w), ratio=float(np.sqrt(max((C / w) ** 2 - 1, 0))),
                dmu=float((cm - (fP + fQ) / 2) / ne), logspread=float(np.log(sp(Q) / sp(P))))


for st in SETS:
    if not os.path.isdir(f'{SD}/{st}') or len(glob.glob(f'{SD}/{st}/*.pdb')) < 82:
        os.makedirs(SD, exist_ok=True); zipfile.ZipFile(f'{ROOT}/raw/models/alphaflow/samples/{st}.zip').extractall(SD); log('풀기', st, len(glob.glob(f'{SD}/{st}/*.pdb')))
    os.makedirs(f'{OUT}/{st}', exist_ok=True)
names = sorted(os.listdir(AD)); log('표적', len(names))
for name in names:
    todo = [st for st in SETS if not os.path.exists(f'{OUT}/{st}/{name}.json')]
    if not todo: continue
    t0 = time.time(); pdb = f'{AD}/{name}/{name}.pdb'; s_md, i_md = at_atoms(md.load(pdb).topology)
    reps = [md.load(f'{AD}/{name}/{name}_prod_R{k}_fit.xtc', top=pdb) for k in (1, 2, 3)]
    rng_md = np.random.default_rng(zlib.crc32(f'7b_rank|md|{name}'.encode()))
    picks = [np.sort(rng_md.choice(r.n_frames, min(N, r.n_frames), replace=False)) for r in reps]
    for st in todo:
        try:
            t = md.load(f'{SD}/{st}/{name}.pdb'); s_m, i_m = at_atoms(t.topology); a0, b0, Lb = block(s_md, s_m)
            if Lb < 0.9 * len(s_md):
                json.dump(dict(name=name, set=st, skip=f'서열 맞춤 {Lb}/{len(s_md)}'), open(f'{OUT}/{st}/{name}.json', 'w')); continue
            imd, im = i_md[a0:a0 + Lb], i_m[b0:b0 + Lb]; L = Lb; ne = L - 1
            MDX = [r.xyz[:, imd.ravel()].reshape(r.n_frames, L, 3, 3) * 10 for r in reps]
            pool = np.concatenate([X[::max(1, len(X) // 200)] for X in MDX])[:600]; ref = fs.mean_point(frames_t(pool))
            A_ = lambda X: fs.align_to_reference(frames_t(X), ref).numpy().reshape(len(X), -1)
            P = [A_(X[p]) for X, p in zip(MDX, picks)]
            Am = A_(t.xyz[:, im.ravel()].reshape(t.n_frames, L, 3, 3) * 10)
            rng = np.random.default_rng(zlib.crc32(f'7b_rank|{st}|{name}'.encode())); outs = []
            for s in range(NSPLIT):
                perm = rng.permutation(len(Am)); h = len(Am) // 2; outs.append(edist_runs(P, [Am[perm[:h]], Am[perm[h:2 * h]]], ne))
            res = {k: float(np.mean([o[k] for o in outs])) for k in outs[0]}
            res.update(name=name, set=st, L=int(L), L_md=len(s_md), L_model=len(s_m), block=[int(a0), int(b0), int(Lb)], n_model=int(len(Am)))
            json.dump(res, open(f'{OUT}/{st}/{name}.json', 'w'))
        except Exception as e:
            log('ERR', st, name, repr(e)[:200])
    log(name, round(time.time() - t0, 1), 's')
log('에너지 거리 끝')

# ---- 순위: 표준 지표 텍스트가 모두 생길 때까지 기다림 (q_7b_rank.sh) ----
for _ in range(24 * 60):
    if all(os.path.exists(p) and 'Per target RMSF r' in open(p).read() for p in TXT.values()): break
    time.sleep(60)
def parse(p):
    out = {}
    for l in open(p):
        m = re.match(r'^(.+?)\s{2,}(-?[0-9.]+)\s*$', l.rstrip())
        if m: out[m.group(1).strip()] = float(m.group(2))
    return out
STD = {st: parse(TXT[st]) for st in SETS if os.path.exists(TXT[st])}
HIGH = ['Pairwise RMSD r', 'Global RMSF r', 'Per target RMSF r', 'PC sim > 0.5 %', 'Weak contacts J', 'Transient contacts J', 'Exposed residue J', 'Exposed MI matrix rho']
LOW = ['RMWD', 'RMWD trans', 'RMWD var', 'MD PCA W2', 'Joint PCA W2']
NEAR = {'Pairwise RMSD': 'MD pairwise RMSD', 'RMSF': 'MD RMSF'}
RNG = np.random.default_rng(20261026); ed = {}
for st in SETS:
    J = [json.load(open(f)) for f in sorted(glob.glob(f'{OUT}/{st}/*.json'))]; J = [j for j in J if 'skip' not in j]
    v = np.array([j['d2n'] for j in J]); bs = [np.median(RNG.choice(v, len(v))) for _ in range(2000)] if len(v) else []
    ed[st] = dict(n=len(v), d2n_median=float(np.median(v)) if len(v) else None, ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if len(v) else None,
                  dmu_median=float(np.median([j['dmu'] for j in J])) if J else None, logspread_median=float(np.median([j['logspread'] for j in J])) if J else None,
                  ratio_median=float(np.median([j['ratio'] for j in J])) if J else None, skipped=sum(1 for f in glob.glob(f'{OUT}/{st}/*.json') if 'skip' in json.load(open(f))))
sets_ok = [st for st in SETS if st in STD and ed[st]['d2n_median'] is not None]
score_ed = [-ed[st]['d2n_median'] for st in sets_ok]; tau = {}
for m in HIGH + LOW + list(NEAR):
    if not all(m in STD[st] for st in sets_ok): continue
    if m in HIGH: sc = [STD[st][m] for st in sets_ok]
    elif m in LOW: sc = [-STD[st][m] for st in sets_ok]
    else: sc = [-abs(STD[st][m] - STD[st][NEAR[m]]) for st in sets_ok]
    t_ = kendalltau(sc, score_ed); tau[m] = dict(tau=None if np.isnan(t_.correlation) else float(t_.correlation), values={st: STD[st][m] for st in sets_ok})
out = dict(written=time.strftime('%F %T'), sets=sets_ok, energy_distance=ed, kendall_tau_vs_energy=tau,
           rule='prereg 7b 모델 순위(서술): 세트마다 표준 지표 중앙값과 에너지 거리 Δ²/√변 중앙값, 지표마다 Kendall τ. 순위표 자체를 주장하지 않는다.')
P_ = RES + '/r38_7b_rank.json'; json.dump(out, open(P_, 'w'), ensure_ascii=False, indent=1); backup(P_)
Lines = [f"7b 모델 순위 ({out['written']}) — 세트 {len(sets_ok)}",
         '에너지 거리 Δ²/√변 중앙값 [95%] (작을수록 MD 에 가까움): ' + '; '.join(f"{st} {ed[st]['d2n_median']:.4f} [{ed[st]['ci'][0]:.4f}, {ed[st]['ci'][1]:.4f}] (n {ed[st]['n']}, 빠짐 {ed[st]['skipped']})" for st in sets_ok),
         '평균 이동 Δμ²/변 · 로그 퍼짐 비 중앙: ' + '; '.join(f"{st} {ed[st]['dmu_median']:.4f} · {ed[st]['logspread_median']:+.2f}" for st in sets_ok),
         '지표별 Kendall τ (에너지 거리 순위와):'] + [f"  {m:24s} τ = {('%.2f' % v['tau']) if v['tau'] is not None else '-'}  " + ', '.join(f"{st.split('_202402')[0]} {x:.2f}" for st, x in v['values'].items()) for m, v in tau.items()]
open(RES + '/r38_7b_rank.txt', 'w').write('\n'.join(Lines) + '\n'); backup(RES + '/r38_7b_rank.txt'); log('\n'.join(Lines)); log('B7_RANK_DONE')
