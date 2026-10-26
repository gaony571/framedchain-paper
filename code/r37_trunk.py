"""(prereg round37) R0 띠 + R1 뒤 절반 판. derived/r37_trunk/<dom>.json, results/r37_trunk.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
import pandas as pd
OUT = ROOT + '/derived/r37_trunk'; os.makedirs(OUT, exist_ok=True)
for d in st['domains']:
    if os.path.exists(f'{OUT}/{d}.json'): continue
    try:
        rng = np.random.default_rng(37); D, _ = load(d); fold = {T: np.nonzero(D[T]['Q'] >= QT)[0] for T in TEMPS}
        if len(fold['320']) < 100: continue
        ref = fs.mean_point(frames(D['320']['X'][fold['320'][:600]])); res = dict(dom=d)
        full, half = {}, {}
        for T in TEMPS:
            if len(fold[T]) == 0: continue
            Y = fs.align_to_reference(frames(D[T]['X'][fold[T]]), ref).numpy(); ne = Y.shape[1]; F = Y.reshape(len(Y), -1); rp = D[T]['rep'][fold[T]]
            full[T], half[T] = {}, {}
            for r in np.unique(rp):
                ix = np.nonzero(rp == r)[0]
                if len(ix) >= 40: full[T][int(r)] = F[ix]
                if len(ix) // 2 >= 40: half[T][int(r)] = F[ix[len(ix) // 2:]]
        res['ne'] = int(ne); N = 100; pk = lambda A: A[sub(np.arange(len(A)), N, rng)]
        A_, B_ = [r for r in (0, 1) if r in full['320']], [r for r in (2, 3, 4) if r in full['320']]
        if len(A_) == 2 and len(B_) >= 2:
            cr = np.mean([edist(pk(full['320'][i]), pk(full['320'][j])) for i in A_ for j in B_]); fa = edist(pk(full['320'][0]), pk(full['320'][1]))
            fb = np.mean([edist(pk(full['320'][i]), pk(full['320'][j])) for i, j in itertools.combinations(B_, 2)])
            res['null320'] = float((cr - (fa + fb) / 2) / np.sqrt(ne))
        H = half
        if len(H.get('320', {})) >= 2:
            r320 = sorted(H['320']); f1 = np.mean([edist(pk(H['320'][i]), pk(H['320'][j])) for i, j in itertools.combinations(r320, 2)][:10]); res['half'] = {}
            for T in TEMPS[1:]:
                if len(H.get(T, {})) < 2: continue
                rT = sorted(H[T]); prs = [(i, j) for i in r320 for j in rT if i != j][:10]
                cr = np.mean([edist(pk(H['320'][i]), pk(H[T][j])) for i, j in prs]); fT = np.mean([edist(pk(H[T][i]), pk(H[T][j])) for i, j in itertools.combinations(rT, 2)][:10])
                res['half'][T] = float((cr - (f1 + fT) / 2) / np.sqrt(ne))
        json.dump(res, open(f'{OUT}/{d}.json', 'w')); print(time.strftime('%H:%M:%S'), d, flush=True)
    except Exception as e: print('ERR', d, repr(e)[:150], flush=True)
cathH = {}
for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
    if l.startswith('#') or not l.strip(): continue
    f = l.split(); cathH[f[0]] = '.'.join(f[1:5])
J = [json.load(open(q)) for q in glob.glob(OUT + '/*.json')]; rng = np.random.default_rng(3); S = dict(n=len(J))
def sfm(vals, ds, n=2000):
    s = pd.Series(vals).groupby([cathH.get(x, x) for x in ds]).median(); b = [np.median(s.values[rng.integers(0, len(s), len(s))]) for _ in range(n)]
    return [round(float(s.median()), 4)] + [round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])] + [int(len(s))]
nul = np.array([j['null320'] for j in J if 'null320' in j]); S['null320'] = dict(n=len(nul), median=round(float(np.median(nul)), 5))
S['band'] = {}
for nm, n in (('+28K', 314), ('+59K', 278), ('+93K', 218), ('+130K', 115), ('n50', 50), ('n150', 150)):
    b = [np.median(nul[rng.integers(0, len(nul), n)]) for _ in range(4000)]; S['band'][nm] = dict(n=n, band95=[round(float(x), 5) for x in np.percentile(b, [2.5, 97.5])])
S['half'] = {}
for T in TEMPS[1:]:
    sb = [j for j in J if T in j.get('half', {})]
    if sb: S['half'][T] = dict(n=len(sb), d2n=sfm([j['half'][T] for j in sb], [j['dom'] for j in sb]))
S['full_for_ratio'] = {'348': 0.0029, '379': 0.0199, '413': 0.0404, '450': 0.0657}
S['written'] = time.strftime('%F %T'); json.dump(S, open(ROOT + '/results/r37_trunk.json', 'w'), indent=1); print(json.dumps(S, indent=1))
open(ROOT + '/results/logs/r37_trunk.done', 'w').write(time.strftime('%F %T'))
