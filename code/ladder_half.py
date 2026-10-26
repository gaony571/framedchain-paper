"""(prereg round34 ③) 비평형 교란 대조: 같은 온도에서 궤적 앞 반쪽 대 뒤 반쪽을 두 '조건' 으로 둔 Δ².
사다리 온도 행과 같은 통계(ladder_temp_single_v2): 교차 = 앞 반쪽 rep i 대 뒤 반쪽 rep j (i ≠ j), 바닥 = 같은 반쪽의 rep 사이.
Δ²_half = 교차 − (바닥_앞 + 바닥_뒤)/2, √변 수로 나눔. 평형이면 기대값 0. derived/ladder_half/<dom>.json, results/ladder_half.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools
ROOT = sys.argv[1]; NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 400; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
import pandas as pd
OUT = ROOT + '/derived/ladder_half'; os.makedirs(OUT, exist_ok=True)
doms = list(st['domains']); np.random.default_rng(34).shuffle(doms); doms = doms[:NMAX]
for d in doms:
    if os.path.exists(f'{OUT}/{d}.json'): continue
    try:
        rng = np.random.default_rng(17); D, _ = load(d); fold = {T: np.nonzero(D[T]['Q'] >= QT)[0] for T in TEMPS}
        if len(fold['320']) < 100: continue
        ref = fs.mean_point(frames(D['320']['X'][fold['320'][:600]])); res = dict(dom=d, half={})
        for T in TEMPS:
            if len(fold[T]) == 0: continue
            Y = fs.align_to_reference(frames(D[T]['X'][fold[T]]), ref).numpy(); ne = Y.shape[1]; F = Y.reshape(len(Y), -1); rp = D[T]['rep'][fold[T]]
            H1, H2 = {}, {}
            for r in np.unique(rp):
                ix = np.nonzero(rp == r)[0]  # 반복 안의 저장 순서 = 시간 순서; 접힌 프레임만 남긴 뒤 반으로 나눔
                if len(ix) < 80: continue
                h = len(ix) // 2; H1[int(r)], H2[int(r)] = F[ix[:h]], F[ix[h:]]
            R_ = sorted(H1)
            if len(R_) < 2: continue
            N = 50; pk = lambda A: A[sub(np.arange(len(A)), N, rng)]
            prs = [(i, j) for i in R_ for j in R_ if i != j][:10]; cps = list(itertools.combinations(R_, 2))[:10]
            cross = float(np.mean([edist(pk(H1[i]), pk(H2[j])) for i, j in prs]))
            f1 = float(np.mean([edist(pk(H1[i]), pk(H1[j])) for i, j in cps])); f2 = float(np.mean([edist(pk(H2[i]), pk(H2[j])) for i, j in cps]))
            res['half'][T] = dict(ne=int(ne), cross=cross, floor1=f1, floor2=f2, d2n=(cross - (f1 + f2) / 2) / np.sqrt(ne), n_reps=len(R_))
        json.dump(res, open(f'{OUT}/{d}.json', 'w')); print(time.strftime('%H:%M:%S'), d, flush=True)
    except Exception as e: print('ERR', d, repr(e)[:150], flush=True)
cath = {}
for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
    if l.startswith('#') or not l.strip(): continue
    f = l.split(); cath[f[0]] = '.'.join(f[1:5])
J = [json.load(open(q)) for q in glob.glob(OUT + '/*.json')]; rng = np.random.default_rng(3); S = dict(n=len(J))
def sfb(vals, ds):
    s = pd.Series(vals).groupby([cath.get(x, x) for x in ds]).median(); b = [np.median(s.values[rng.integers(0, len(s), len(s))]) for _ in range(2000)]
    return [round(float(s.median()), 4)] + [round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])] + [int(len(s))]
for T in TEMPS:
    sb = [j for j in J if T in j['half']]
    if sb: S[T] = dict(n=len(sb), d2n_half=sfb([j['half'][T]['d2n'] for j in sb], [j['dom'] for j in sb]),
                       frac_pos=round(float(np.mean([j['half'][T]['d2n'] > 0 for j in sb])), 3))
S['compare_rows'] = {'+28K': 0.003, '+59K': 0.020, '+93K': 0.040, '+130K': 0.066}
S['note'] = '[상과 중앙값, 95% 하한, 상한, 상과 수]. 표본 50 프레임(사다리 행은 100 — 표본 수가 다르면 Δ² 의 분산만 다르고 기대값은 같다).'
S['written'] = time.strftime('%F %T'); json.dump(S, open(ROOT + '/results/ladder_half.json', 'w'), indent=1); print(json.dumps(S, indent=1))
open(ROOT + '/results/logs/ladder_half.done', 'w').write(time.strftime('%F %T'))
