
"""사다리 온도 행 단일 반복 바닥·절대 눈금 (prereg round27 (가)). results/ladder_temp_single.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
OUT = ROOT + '/derived/ladder_temp_single'; os.makedirs(OUT, exist_ok=True)
for d in st['domains']:
    if os.path.exists(f'{OUT}/{d}.json'): continue
    try:
        rng = np.random.default_rng(17); D, _ = load(d); fold = {T: np.nonzero(D[T]['Q'] >= QT)[0] for T in TEMPS}
        if len(fold['320']) < 100: continue
        ref = fs.mean_point(frames(D['320']['X'][fold['320'][:600]]))
        reps = {}
        for T in TEMPS:
            if len(fold[T]) == 0: continue
            Y = fs.align_to_reference(frames(D[T]['X'][fold[T]]), ref).numpy(); ne = Y.shape[1]; F = Y.reshape(len(Y), -1); rp = D[T]['rep'][fold[T]]
            reps[T] = {int(r): F[rp == r] for r in np.unique(rp) if (rp == r).sum() >= 40}
        if len(reps.get('320', {})) < 2: continue
        N = 100; pk = lambda A: A[sub(np.arange(len(A)), N, rng)]
        r320 = sorted(reps['320']); fl = [edist(pk(reps['320'][i]), pk(reps['320'][j])) for i, j in itertools.combinations(r320, 2)][:10]
        res = dict(dom=d, ne=int(ne), floor1=float(np.mean(fl)), n_floor_pairs=len(fl), cross1={})
        for T in TEMPS[1:]:
            if T not in reps or not reps[T]: continue
            prs = [(i, j) for i in r320 for j in sorted(reps[T]) if i != j][:10]
            if prs: res['cross1'][T] = float(np.mean([edist(pk(reps['320'][i]), pk(reps[T][j])) for i, j in prs]))
        json.dump(res, open(f'{OUT}/{d}.json', 'w')); print(time.strftime('%H:%M:%S'), d, flush=True)
    except Exception as e: print('ERR', d, repr(e)[:150], flush=True)
cath = {}
for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
    if l.startswith('#') or not l.strip(): continue
    f = l.split(); cath[f[0]] = '.'.join(f[1:5])
J = [json.load(open(q)) for q in glob.glob(OUT + '/*.json')]; rng = np.random.default_rng(3); S = dict(n=len(J))
C = {json.load(open(q))['dom']: json.load(open(q)).get('ladder', {}) for q in glob.glob(ROOT + '/derived/c1_core/*.json')}
def sfb(vals, doms):
    s = pd.Series(vals).groupby([cath.get(x, x) for x in doms]).median(); b = [np.median(s.values[rng.integers(0, len(s), len(s))]) for _ in range(2000)]
    return [round(float(s.median()), 3)] + [round(float(x), 3) for x in np.percentile(b, [2.5, 97.5])] + [int(len(s))]
import pandas as pd
S['floor1_abs'] = sfb([j['floor1'] / np.sqrt(j['ne']) for j in J], [j['dom'] for j in J])
pooled_floor = [C[j['dom']]['floor_320']['ed'] for j in J if 'floor_320' in C.get(j['dom'], {})]
S['floor_pooled_over_single_median'] = round(float(np.median([C[j['dom']]['floor_320']['ed'] / j['floor1'] for j in J if 'floor_320' in C.get(j['dom'], {})])), 3)
for T in TEMPS[1:]:
    sb = [j for j in J if T in j['cross1']]
    S[T] = dict(n=len(sb), multiple_single=sfb([j['cross1'][T] / j['floor1'] for j in sb], [j['dom'] for j in sb]), cross_abs=sfb([j['cross1'][T] / np.sqrt(j['ne']) for j in sb], [j['dom'] for j in sb]),
                frac_cross_gt_floor=round(float(np.mean([j['cross1'][T] > j['floor1'] for j in sb])), 3))
S['note'] = '배수 = 단일 반복 교차(320 rep i 대 T rep j, i≠j) / 단일 반복 바닥(320 rep i 대 j); 절대 = 에너지 거리/√변 수. [중앙값, 95% 하한, 상한, 상과 수].'
json.dump(S, open(ROOT + '/results/ladder_temp_single.json', 'w'), indent=1); print(json.dumps(S, indent=1))
open(ROOT + '/results/logs/ladder_temp_single.done', 'w').write(time.strftime('%F %T'))
