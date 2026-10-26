
"""사다리 온도 행 단일 반복 바닥·절대 눈금 (prereg round27 (가)). results/ladder_temp_single.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
OUT = ROOT + '/derived/ladder_temp_single_v2'; os.makedirs(OUT, exist_ok=True)
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
        res = dict(dom=d, ne=int(ne), floor1=float(np.mean(fl)), n_floor_pairs=len(fl), cross1={}, floorT={})
        for T in TEMPS[1:]:
            if T not in reps or not reps[T]: continue
            prs = [(i, j) for i in r320 for j in sorted(reps[T]) if i != j][:10]
            if prs: res['cross1'][T] = float(np.mean([edist(pk(reps['320'][i]), pk(reps[T][j])) for i, j in prs]))
            rT = sorted(reps[T])
            if len(rT) >= 2: res['floorT'][T] = float(np.mean([edist(pk(reps[T][i]), pk(reps[T][j])) for i, j in itertools.combinations(rT, 2)][:10]))
        json.dump(res, open(f'{OUT}/{d}.json', 'w')); print(time.strftime('%H:%M:%S'), d, flush=True)
    except Exception as e: print('ERR', d, repr(e)[:150], flush=True)
open(ROOT + '/results/logs/ladder_temp_single_v2.done', 'w').write(time.strftime('%F %T'))
