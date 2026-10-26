
import sys, json, time
ROOT = sys.argv[1]
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
h = pd.read_csv(ROOT + '/results/hu1_rulers_1674.tsv', sep='\t').dropna(subset=['rmsd_tma', 'n_ali'])
n, m = h.n.values.astype(float), h.m.values.astype(float); Lb = 0.5 * (n + m); Lmin = np.minimum(n, m)
r, N = h.rmsd_tma.values, h.n_ali.values.astype(float)
# Kolodny·Koehl·Levitt 2005 식 (1)–(3), 낮을수록 좋음 -> 부호 뒤집어 유사도로
kkl = {'-SI': -(r * Lmin / N), '-SAS': -(100 * r / N), '-MI': (1 + N) / ((1 + r / 1.5) * (1 + Lmin)) - 1}
ours = {'D_TM/sqrtLbar': -np.sqrt(np.clip(Lb - h.S_tm4.values, 0, None) / Lb)}
for l in ['4', '8', '16', 'inf']: ours['l=' + l] = -np.sqrt(np.clip(h['D2_' + ('inf' if l == 'inf' else 'l' + l)].values, 0, None) / (0.75 * (n + m)))
old = {'TM_long': np.minimum(h.tm_a, h.tm_b).values, 'TM_short': np.maximum(h.tm_a, h.tm_b).values, 'GDT': h.gdt_ts_B.values,
       'lDDT': h.lddt_B.values, '-RMSD': -r, '-ESA': -h.esa.values, 'N_ali/Lmin': N / Lmin}
tab = pd.DataFrame({k: {o: spearmanr(v, w)[0] for o, w in ours.items()} for k, v in {**kkl, **old}.items()}).T.round(3)
cross = pd.DataFrame({k: {o: spearmanr(v, w)[0] for o, w in {**old}.items()} for k, v in kkl.items()}).T.round(3)
lr = np.maximum(n, m) / Lmin
def pcorr(a, b):
    Z = np.c_[np.ones(len(a)), rankdata(lr), rankdata(Lb)]
    ra, rb = rankdata(a), rankdata(b)
    ea = ra - Z @ np.linalg.lstsq(Z, ra, rcond=None)[0]; eb = rb - Z @ np.linalg.lstsq(Z, rb, rcond=None)[0]
    return float(np.corrcoef(ea, eb)[0, 1])
part = {k: round(pcorr(v, old['TM_long']), 3) for k, v in {**kkl, **ours}.items()}
strata = {}
for nm, msk in [('lr<=1.25', lr <= 1.25), ('lr>1.5', lr > 1.5)]:
    strata[nm] = {k: round(spearmanr(v[msk], old['TM_long'][msk])[0], 3) for k, v in {**kkl, 'D_TM/sqrtLbar': ours['D_TM/sqrtLbar']}.items()}
    strata[nm]['n'] = int(msk.sum())
out = dict(n_pairs=len(h), vs_ours=tab.to_dict(), kkl_vs_old=cross.to_dict(), partial_vs_TM_long=part, by_length_ratio=strata, written=time.strftime('%F %T'))
json.dump(out, open(ROOT + '/results/hu1_kkl_measures.json', 'w'), indent=1, ensure_ascii=False)
print(len(h)); print(tab.to_string()); print(cross.to_string()); print(part); print(strata)
