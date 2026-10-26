
"""사다리 성분 (prereg round32 ladder_components): (a) 평균 이동 Δμ²/변 수 (쌍선형 바닥 빼기), (b) 로그 총퍼짐 비. results/ladder_components.json, derived/ladder_components_temp.tsv."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
import pandas as pd
rng0 = np.random.default_rng(32)
def mci(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; b = [np.median(v[rng0.integers(0, len(v), len(v))]) for _ in range(2000)]
    return dict(n=int(len(v)), median=round(float(np.median(v)), 4), CI95=[round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])])
# ---- 온도 행 (mdcath_5T, 접힌 프레임, 반복마다 평균 점·퍼짐) ----
rows = []
for d in st['domains']:
    try:
        D, _ = load(d); fold = {T: np.nonzero(D[T]['Q'] >= QT)[0] for T in TEMPS}
        if len(fold['320']) < 100: continue
        ref = fs.mean_point(frames(D['320']['X'][fold['320'][:600]])); M, S = {}, {}
        for T in TEMPS:
            if len(fold[T]) == 0: continue
            Y = fs.align_to_reference(frames(D[T]['X'][fold[T]]), ref).numpy(); ne = Y.shape[1]; F = Y.reshape(len(Y), -1); rp = D[T]['rep'][fold[T]]
            M[T], S[T] = {}, {}
            for r in np.unique(rp):
                Z = F[rp == r]
                if len(Z) >= 40: M[T][int(r)] = Z.mean(0); S[T][int(r)] = float(((Z - Z.mean(0)) ** 2).sum(1).mean())
        if len(M.get('320', {})) < 2: continue
        fl = lambda T: np.mean([((M[T][i] - M[T][j]) ** 2).sum() for i, j in itertools.combinations(sorted(M[T]), 2)]) if len(M[T]) >= 2 else np.nan
        f320 = fl('320')
        for T in TEMPS[1:]:
            if T not in M or len(M[T]) < 2: continue
            cr = np.mean([((M['320'][i] - M[T][j]) ** 2).sum() for i in M['320'] for j in M[T] if i != j])
            d_mu = (cr - (f320 + fl(T)) / 2) / ne
            ls = float(np.log(np.mean(list(S[T].values())) / np.mean(list(S['320'].values()))))
            r3 = sorted(S['320']); h = len(r3) // 2; ls_floor = abs(float(np.log(np.mean([S['320'][k] for k in r3[:h]]) / np.mean([S['320'][k] for k in r3[h:]]))))
            rows.append(dict(dom=d, T=T, ne=ne, d_mu=d_mu, floor_mu=(f320 / 2) / ne, logspread=ls, logspread_floor=ls_floor))
        print(time.strftime('%H:%M:%S'), d, flush=True)
    except Exception as e: print('ERR', d, repr(e)[:150], flush=True)
Dt = pd.DataFrame(rows); Dt.to_csv(ROOT + '/derived/ladder_components_temp.tsv', sep='\t', index=False)
OUT = {'temperature': {f'+{int(T) - 320}K': dict(d_mu=mci(g.d_mu), logspread=mci(g.logspread), floor_mu_median=round(float(g.floor_mu.median()), 4), logspread_floor_median=round(float(g.logspread_floor.median()), 4)) for T, g in Dt.groupby('T')}}
# ---- 모델 행 (ladder_model_error_folded) ----
MODELS = ['alphaflow_md_templates_base_202402', 'esmflow_md_templates_base_202402', 'alphaflow_md_base_202402', 'esmflow_md_base_202402']
mr = {m: [] for m in MODELS}
for p in glob.glob(ROOT + '/derived/ladder_model_error_folded/*.json'):
    j = json.load(open(p))
    if 'floor' not in j: continue
    ne = j['n_common_res'] - 1
    for m in MODELS:
        q = j.get(m)
        if not q or 'skip' in q or 'mean2' not in q: continue
        n_used = max(min(int(round(q['n_samples'] * q['fracfold'])), 250), 2)
        d_mu = (q['mean2'] - j['floor']['mean2'] / 2 - q['spreadA'] / n_used) / ne
        mr[m].append(dict(d_mu=d_mu, floor_mu=(j['floor']['mean2'] / 2) / ne, logspread=float(np.log(q['spreadA'] / q['spreadB'])), logspread_floor=abs(float(np.log(j['floor']['spreadA'] / j['floor']['spreadB'])))))
OUT['model'] = {m: dict(n=len(v), d_mu=mci([x['d_mu'] for x in v]), logspread=mci([x['logspread'] for x in v]), floor_mu_median=round(float(np.median([x['floor_mu'] for x in v])), 4), logspread_floor_median=round(float(np.median([x['logspread_floor'] for x in v])), 4)) for m, v in mr.items() if v}
# ---- 역장 행 (ladder_ff_folded, 반쪽 묶음 대칭) ----
fr = []
for p in glob.glob(ROOT + '/derived/ladder_ff_folded/*.json'):
    j = json.load(open(p))
    if 'cross' not in j: continue
    ne = j['n_res_matched'] - 1; cr = (j['cross']['mean2'] + j['cross2']['mean2']) / 2
    d_mu = (cr - (j['floor_sm']['mean2'] + j['floor_bio']['mean2']) / 2) / ne
    ls = float(np.log(((j['cross']['spreadB'] + j['cross2']['spreadB']) / 2) / ((j['cross']['spreadA'] + j['cross2']['spreadA']) / 2)))
    fr.append(dict(d_mu=d_mu, floor_mu=((j['floor_sm']['mean2'] + j['floor_bio']['mean2']) / 4) / ne, logspread=ls, logspread_floor=abs(float(np.log(j['floor_sm']['spreadA'] / j['floor_sm']['spreadB'])))))
OUT['forcefield'] = dict(n=len(fr), d_mu=mci([x['d_mu'] for x in fr]), logspread=mci([x['logspread'] for x in fr]), floor_mu_median=round(float(np.median([x['floor_mu'] for x in fr])), 4), logspread_floor_median=round(float(np.median([x['logspread_floor'] for x in fr])), 4), note='logspread = log(BioEmu/SM); 반쪽 묶음끼리.')
OUT['note'] = 'd_mu = 평균 이동² / 변 수 (교차 − 바닥 평균; 모델 행은 모델 표본 평균의 표본 분산 tr Σ/n 도 뺌). floor_mu = 실행 하나의 평균 점 흔들림(바닥/2)/변 수 — 비교 기준. logspread = log(총퍼짐 비), 온도는 T/320, 모델은 모델/MD, 역장은 BioEmu/SM.'
OUT['written'] = time.strftime('%F %T'); json.dump(OUT, open(ROOT + '/results/ladder_components.json', 'w'), indent=1); print(json.dumps(OUT, indent=1))
open(ROOT + '/results/logs/ladder_components.done', 'w').write(time.strftime('%F %T'))
