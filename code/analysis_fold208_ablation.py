
import sys, os, json, time
ROOT = sys.argv[1]; OUTJ = sys.argv[2] if len(sys.argv) > 2 else ROOT + '/results/pu_fold208_ablation.json'
import numpy as np, pandas as pd
st = json.load(open(ROOT + '/results/pu_fold208_setup.json')); blk = st['block']; bi = {n: i for i, n in enumerate(blk)}; N = len(blk)
pf = pd.read_csv(ROOT + '/results/pu_fold208_pairs.tsv', sep='\t')
xb = pd.read_csv(ROOT + '/results/xprot_bench.tsv', sep='\t').set_index('PDB')
topo = np.array([xb.loc[n, 'CATH_ID'].rsplit('.', 1)[0] for n in blk]); cls = np.array([str(xb.loc[n, 'CATH_class']) for n in blk])
ia, ib = pf.a.map(bi).values, pf.b.map(bi).values; n, m = pf.n.values.astype(float), pf.m.values.astype(float); Lb = 0.5 * (n + m)
def full(va, vb=None):
    M = np.full((N, N), np.nan); M[ia, ib] = va; M[ib, ia] = va if vb is None else vb; np.fill_diagonal(M, np.nan); return M
scores = {'TM-align (query-normalized)': full(pf.tm_a.values, pf.tm_b.values),
          'D_TM distance (fixed 4 A)': full(-np.sqrt(np.clip(Lb - pf.S_fix_hi.values, 0, None))),
          'D_TM normalized (d0(Lbar), sqrt(1-S/Lbar))': full(-np.sqrt(np.clip(1 - pf.S_avg_hi.values / Lb, 0, None))),
          'TM-align max(tm_a,tm_b) (short-normalized)': full(np.maximum(pf.tm_a.values, pf.tm_b.values)),
          'TM-align min(tm_a,tm_b) (long-normalized)': full(np.minimum(pf.tm_a.values, pf.tm_b.values)),
          'TM-align approx S/Lbar ((tm_a n + tm_b m)/2 / Lbar)': full(0.5 * (pf.tm_a.values * n + pf.tm_b.values * m) / Lb),
          'ours S_avg_hi / query length': full(pf.S_avg_hi.values / n, pf.S_avg_hi.values / m),
          'ours S_fix_hi / Lbar (fixed 4 A, normalized)': full(pf.S_fix_hi.values / Lb),
          'TM-align mean(tm_a,tm_b) (old 0.285 convention)': full(0.5 * (pf.tm_a.values + pf.tm_b.values))}
done_frac = len(pf) / (N * (N - 1) / 2)
queries = [i for i in range(N) if (topo == topo[i]).sum() > 1]
def ap_per_query(S):
    out = []
    for q in queries:
        s = S[q].copy(); ok = ~np.isnan(s); ok[q] = False
        idx = np.where(ok)[0]
        if len(idx) == 0: out.append(np.nan); continue
        order = idx[np.argsort(-s[idx])]; rel = (topo[order] == topo[q]).astype(float)
        if rel.sum() == 0: out.append(np.nan); continue
        prec = np.cumsum(rel) / np.arange(1, len(rel) + 1); out.append(float((prec * rel).sum() / rel.sum()))
    return np.array(out)
def nn1(S):
    acc = []
    for q in range(N):
        s = S[q].copy(); s[q] = np.nan
        if np.all(np.isnan(s)): continue
        acc.append(cls[int(np.nanargmax(s))] == cls[q])
    return float(np.mean(acc))
res = {'pairs_done_fraction': round(done_frac, 4), 'n_queries': len(queries), 'scores': {}}
aps = {k: ap_per_query(S) for k, S in scores.items()}
rng = np.random.default_rng(0); base = aps['TM-align (query-normalized)']
for k, a in aps.items():
    ok = ~np.isnan(a) & ~np.isnan(base); d = a[ok] - base[ok]
    bs = [np.mean(d[rng.integers(0, len(d), len(d))]) for _ in range(2000)]
    res['scores'][k] = dict(topology_mAP=round(float(np.nanmean(a)), 4), class_1NN=round(nn1(scores[k]), 4), diff_vs_TMalign=round(float(d.mean()), 4), diff_CI95=[round(float(x), 4) for x in np.percentile(bs, [2.5, 97.5])])
ref = aps['TM-align approx S/Lbar ((tm_a n + tm_b m)/2 / Lbar)']; res['vs_TMalign_Lbar'] = {}
for k, x in aps.items():
    ok = ~np.isnan(x) & ~np.isnan(ref); d = x[ok] - ref[ok]
    bs = [np.mean(d[rng.integers(0, len(d), len(d))]) for _ in range(2000)]
    res['vs_TMalign_Lbar'][k] = dict(diff=round(float(d.mean()), 4), CI95=[round(float(v), 4) for v in np.percentile(bs, [2.5, 97.5])])
res['written'] = time.strftime('%Y-%m-%d %H:%M:%S')
json.dump(res, open(OUTJ, 'w'), indent=1, ensure_ascii=False); print(json.dumps(res, indent=1, ensure_ascii=False))
