
"""E0 해상도 대조 (prereg round32 E0_resolution). 같은 접힌 프레임에서 세 자 — 프레임 에너지 거리, Cα 에너지 거리, 창 RMSF² (±3) 로그 비 — 를 같은 온도 표지 순열 검정(반복 10 중 5, 252 나눔)으로. 창 RMSF² 는 반복 2 대 3 나눔 귀무 최댓값 판도. results/e0_resolution.json, derived/e0_resolution.tsv."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
import pandas as pd
NF = 40; PERMS = list(itertools.combinations(range(10), 5))
def kabsch_to(P, Q):  # P (f, n, 3) 를 Q (n, 3) 에 맞춤
    Pc = P - P.mean(1, keepdims=True); Qc = Q - Q.mean(0); H = np.einsum('fni,nj->fij', Pc, Qc); U, S, Vt = np.linalg.svd(H); dd = np.sign(np.linalg.det(np.einsum('fij,fjk->fik', U, Vt)))
    Dg = np.ones((len(P), 3)); Dg[:, 2] = dd; Rm = np.einsum('fij,fj,fjk->fik', U, Dg, Vt); return np.einsum('fni,fij->fnj', Pc, Rm)
def ed_from(Dm, ga, gb):
    a, b = np.concatenate(ga), np.concatenate(gb); n, m = len(a), len(b)
    return 2 * Dm[np.ix_(a, b)].mean() - Dm[np.ix_(a, a)].sum() / (n * (n - 1)) - Dm[np.ix_(b, b)].sum() / (m * (m - 1))
rows = []; rng = np.random.default_rng(32)
for d in st['domains']:
    try:
        D, _ = load(d); fold = {T: np.nonzero(D[T]['Q'] >= QT)[0] for T in ('320', '348')}
        reps = {T: [r for r in range(5) if (D[T]['rep'][fold[T]] == r).sum() >= NF] for T in ('320', '348')}
        if len(reps['320']) < 5 or len(reps['348']) < 5: continue
        ref = fs.mean_point(frames(D['320']['X'][fold['320'][:600]]))
        CA0 = D['320']['X'][fold['320'], :, 0]; Q0 = kabsch_to(CA0[:200], CA0[0]).mean(0); Q0 = kabsch_to(CA0[:200], Q0).mean(0)
        blocks, Fr, Ca, rm = [], [], [], []
        for T in ('320', '348'):
            for r in range(5):
                ix = fold[T][D[T]['rep'][fold[T]] == r]; ix = ix[np.linspace(0, len(ix) - 1, NF).round().astype(int)]
                Fr.append(fs.align_to_reference(frames(D[T]['X'][ix]), ref).numpy().reshape(NF, -1)); Ca.append(kabsch_to(D[T]['X'][ix, :, 0], Q0).reshape(NF, -1))
                C = Ca[-1].reshape(NF, -1, 3); r2 = ((C - C.mean(0)) ** 2).sum(-1).mean(0); L = len(r2); w = np.array([r2[max(0, i - 3):i + 4].mean() for i in range(L)]); rm.append(w)
        idx = [np.arange(k * NF, (k + 1) * NF) for k in range(10)]
        out = dict(dom=d)
        for nm, X in (('frame_ed', np.concatenate(Fr)), ('ca_ed', np.concatenate(Ca))):
            Dm = np.sqrt(cost(X, X)); obs = ed_from(Dm, idx[:5], idx[5:]); null = np.array([ed_from(Dm, [idx[k] for k in g], [idx[k] for k in range(10) if k not in g]) for g in PERMS])
            out[nm + '_p'] = float((null >= obs - 1e-12).mean())
        rm = np.array(rm); stat = lambda g, h: float(np.median(np.log(rm[list(h)].mean(0) / rm[list(g)].mean(0))))
        obs = stat(range(5), range(5, 10)); null = np.array([stat(g, [k for k in range(10) if k not in g]) for g in PERMS]); out['wrmsf_p'] = float((null >= obs - 1e-12).mean()); out['wrmsf_obs'] = obs
        sp = [abs(stat(g, [k for k in base if k not in g])) for base in (range(5), range(5, 10)) for g in itertools.combinations(base, 2)]; out['wrmsf_pass_split'] = bool(obs > max(sp))
        rows.append(out); print(time.strftime('%H:%M:%S'), d, flush=True)
    except Exception as e: print('ERR', d, repr(e)[:150], flush=True)
Df = pd.DataFrame(rows); Df['e0'] = Df.dom.isin(E0SET); Df.to_csv(ROOT + '/derived/e0_resolution.tsv', sep='\t', index=False)
def fr(v): v = np.asarray(v, float); b = [v[rng.integers(0, len(v), len(v))].mean() for _ in range(2000)]; return dict(n=int(len(v)), frac=round(float(v.mean()), 3), CI95=[round(float(x), 3) for x in np.percentile(b, [2.5, 97.5])])
OUT = {lab: {k: fr(sb[k] <= 0.05) for k in ('frame_ed_p', 'ca_ed_p', 'wrmsf_p')} | dict(wrmsf_pass_split=fr(sb.wrmsf_pass_split)) for lab, sb in (('all', Df), ('e0_set', Df[Df.e0]))}
OUT['note'] = '온도 표지 순열(반복 10 개 중 5 개를 348 로, 252 나눔), 한쪽 p ≤ 0.05 면 통과. 반복마다 접힌 프레임 40 개. 같은 프레임에서 세 자.'; OUT['written'] = time.strftime('%F %T')
json.dump(OUT, open(ROOT + '/results/e0_resolution.json', 'w'), indent=1); print(json.dumps(OUT, indent=1)); open(ROOT + '/results/logs/e0_resolution.done', 'w').write(time.strftime('%F %T'))
