
"""접힌 프레임 분석 (2026-10-17; 예측은 results/prereg_2026-10-17.json). 입력 derived/mdcath_5T_folded, derived/mdcath_5T (Q).
(i) C2 Q 점검 (ii) C2 깨끗한 검정 + V2 접힌 판 (iii) H-T 외삽 없는 비교 (iv) C3 스칼라 대조."""
import sys, os, json, time, glob
ROOT = sys.argv[1]
import numpy as np, pandas as pd
from tmtools import tm_align
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import stage0 as s0
FD, RAW5 = ROOT + '/derived/mdcath_5T_folded', ROOT + '/derived/mdcath_5T'; rng = np.random.default_rng(31); OUT = {}
lg = lambda x: np.log(np.maximum(x, 1e-8))
def save(): json.dump(OUT, open(ROOT + '/results/ht_own_alignment.json', 'w'), indent=1, ensure_ascii=False)
cath = {}
for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
    if l.startswith('#') or not l.strip(): continue
    f = l.split(); cath[f[0]] = '.'.join(f[1:5])
D = {}
for p in glob.glob(FD + '/*.npz'):
    z = np.load(p); D[os.path.basename(p)[:-4]] = {k: z[k] for k in z.files}
def sfsum(v, H, n=2000):
    v = pd.Series(np.asarray(v, float)); ok = v.notna().values; s = v[ok].groupby(np.asarray(H)[ok]).median()
    bs = np.array([np.median(s.values[rng.integers(0, len(s), len(s))]) for _ in range(n)])
    return dict(n=int(ok.sum()), n_sf=int(len(s)), sf_median=round(float(s.median()), 4), CI95=[round(float(x), 4) for x in np.percentile(bs, [2.5, 97.5])])
def mean_rep(z, key, R_=None):
    a = z[key] if R_ is None else z[key][R_]; return np.nanmean(a, 0)
HT = pd.read_csv(ROOT + '/results/ht_shared_spread.tsv', sep='\t')[['th', 'me', 'H']]
def shared(a, b):
    r = tm_align(a['ca_ref'].astype(float), b['ca_ref'].astype(float), 'A' * len(a['ca_ref']), 'A' * len(b['ca_ref']))
    i = j = 0; mp = {}
    for x, y in zip(r.seqxA, r.seqyA):
        if x != '-' and y != '-': mp[i] = j
        i += x != '-'; j += y != '-'
    et = np.array([k for k in range(len(a['ca_ref']) - 1) if k in mp and k + 1 in mp and mp[k + 1] == mp[k] + 1])
    return et, (np.array([mp[k] for k in et]) if len(et) else et)
def eiv_resid(y, X, rel):
    Xc = X - X.mean(0); yc = y - y.mean(); S = Xc.T @ Xc / len(y); Se = np.diag(np.diag(S) * (1 - np.clip(rel, 0.05, 1)))
    return yc - Xc @ np.linalg.solve(S - Se + 1e-9 * np.eye(len(S)), Xc.T @ yc / len(y))
z_ = lambda v: (v - np.nanmean(v)) / (np.nanstd(v) + 1e-12)

def shared_own(a, b):
    P, Q = a['ca_ref'].astype(float), b['ca_ref'].astype(float); r = s0.tm_align(P, Q, d0=s0.d0_tm(0.5 * (len(P) + len(Q))))
    M = r['match']; ii, jj = np.nonzero(M); mp = dict(zip(ii.tolist(), jj.tolist()))
    et = np.array([k for k in range(len(P) - 1) if k in mp and k + 1 in mp and mp[k + 1] == mp[k] + 1])
    return et, (np.array([mp[k] for k in et]) if len(et) else et)
shared_tma = shared; shared = shared_own
JAC = []
_orig = shared_own
def shared(a, b):
    et, em = _orig(a, b); e2, m2 = shared_tma(a, b)
    A_ = set(zip(et.tolist(), em.tolist())); B_ = set(zip(e2.tolist(), m2.tolist())); JAC.append(len(A_ & B_) / max(len(A_ | B_), 1)); return et, em

rows = []
for th, me, H in HT.itertuples(index=False):
    if th not in D or me not in D or not all(f'{T}_spread_rep' in D[x] for T in ('320', '348') for x in (th, me)): continue
    a, b = D[th], D[me]; et, em = shared(a, b)
    if len(et) < 20: continue
    rows.append(dict(H=H, th348_me320=float(np.nanmean(lg(mean_rep(a, '348_spread_rep')[et]) - lg(mean_rep(b, '320_spread_rep')[em]))),
                     th320_me320=float(np.nanmean(lg(mean_rep(a, '320_spread_rep')[et]) - lg(mean_rep(b, '320_spread_rep')[em]))),
                     th348_me348=float(np.nanmean(lg(mean_rep(a, '348_spread_rep')[et]) - lg(mean_rep(b, '348_spread_rep')[em])))))
HTd = pd.DataFrame(rows); OUT['iii_HT_direct_own_alignment'] = {c: sfsum(HTd[c], HTd.H) for c in HTd.columns if c != 'H'}
OUT['iii_HT_direct_own_alignment']['note'] = '공유 변의 평균 로그 비(접힌 프레임). 예측: th348_me320 ≈ 0 (약간 양~0), th320_me320 < 0.'; save(); print('iii done', flush=True)

OUT['jaccard_shared_edges_vs_tmalign'] = dict(median=round(float(np.median(JAC)), 3), q25=round(float(np.quantile(JAC, 0.25)), 3), n=len(JAC))
OUT['note'] = 'H-T 외삽 없는 비교를 이 공간의 D_TM 정렬(stage0 tm_align, d0(L̄))로 공유 변을 정해 다시 냄 (prereg round27 (다)). 비교: folded_analyses.json iii_HT_direct (TM-align).'
save(); print(json.dumps(OUT, indent=1)); open(ROOT + '/results/logs/ht_own_alignment.done', 'w').write(time.strftime('%F %T'))
