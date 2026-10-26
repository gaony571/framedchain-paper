"""(prereg round34 ①) θ* 음성 대조와 접힌 판. results/ht_theta_meso.json, derived/ht_theta_meso/*.tsv.
(a) 재현: H-T 쌍의 θ* (전체 프레임)를 이 코드로 다시 내 results/ht_theta_star.tsv 와 견줌.
(b) 음성 대조: 같은 상과의 중온–중온 쌍(생육 < 40 °C). 두 방향 모두 θ* 를 내고, 방향을 무작위로 정한 판과 대칭 판의 상과 중앙값.
(c) 접힌 판: H-T 쌍의 θ* 를 접힌 프레임 퍼짐(mdcath_5T_folded)으로, 그리고 두 쪽 모두 네 온도에서 접힌 비율 ≥ 0.8 인 쌍만."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time
ROOT = sys.argv[1]
import numpy as np, pandas as pd
from tmtools import tm_align
rng = np.random.default_rng(34); lg = lambda x: np.log(np.maximum(x, 1e-8)); OUT = {}
TT = (320.0, 348.0, 379.0, 413.0); TS = ('320', '348', '379', '413')
OD = ROOT + '/derived/ht_theta_meso'; os.makedirs(OD, exist_ok=True)
def save(): OUT['written'] = time.strftime('%F %T'); json.dump(OUT, open(ROOT + '/results/ht_theta_meso.json', 'w'), indent=1, ensure_ascii=False)
cath = {}
for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
    if l.startswith('#') or not l.strip(): continue
    f = l.split(); cath[f[0]] = '.'.join(f[1:5])

def getk(z, T):
    for k in (f'{T}_spread', f'{T}_spread_rep'):
        if k in z.files:
            v = z[k]; return v.mean(0) if v.ndim == 2 else v
    return None
def load_full(d):
    out = {}
    for base in ('/derived/mdcath_spread/', '/derived/mdcath_spread_hot/'):
        p = ROOT + base + d + '.npz'
        if not os.path.exists(p): continue
        z = np.load(p)
        if 'ca0' in z.files and 'ca' not in out: out['ca'] = z['ca0'].astype(float)
        for T in TS:
            v = getk(z, T)
            if v is not None and T not in out: out[T] = v
    return out if 'ca' in out and all(T in out for T in TS) else None
def load_fold(d):
    p = ROOT + '/derived/mdcath_5T_folded/' + d + '.npz'
    if not os.path.exists(p): return None
    z = np.load(p); out = {'ca': z['ca_ref'].astype(float)} if 'ca_ref' in z.files else None
    if out is None: return None
    for T in TS:
        k = f'{T}_spread_rep'
        if k not in z.files: return None
        out[T] = np.nanmean(z[k], 0)
    out['ff'] = {T: (float(np.nanmean(z[f'{T}_fracfold'])) if f'{T}_fracfold' in z.files else np.nan) for T in TS}
    return out
def shared(a, b):
    r = tm_align(a, b, 'A' * len(a), 'A' * len(b)); i = j = 0; mp = {}
    for x, y in zip(r.seqxA, r.seqyA):
        if x != '-' and y != '-': mp[i] = j
        i += x != '-'; j += y != '-'
    et = np.array([k for k in range(len(a) - 1) if k in mp and k + 1 in mp and mp[k + 1] == mp[k] + 1])
    return et, (np.array([mp[k] for k in et]) if len(et) else et), float(max(r.tm_norm_chain1, r.tm_norm_chain2)), float(min(r.tm_norm_chain1, r.tm_norm_chain2))
def theta(hot, ref320):
    """hot: 네 온도의 로그 퍼짐 평균, ref320: 기준 쪽 320 K 값. 조각 선형, 범위 밖은 가장 가까운 두 온도의 기울기."""
    h = np.asarray(hot, float); x = np.asarray(TT)
    if not np.all(np.isfinite(h)): return np.nan, 'nan'
    if h[0] >= ref320:
        s = (h[1] - h[0]) / (x[1] - x[0]); return (x[0] + (ref320 - h[0]) / s, 'below') if s > 0 else (np.nan, 'nan')
    for k in range(3):
        if h[k] < ref320 <= h[k + 1]: return x[k] + (ref320 - h[k]) * (x[k + 1] - x[k]) / (h[k + 1] - h[k]), 'in'
    s = (h[3] - h[2]) / (x[3] - x[2]); return (x[3] + (ref320 - h[3]) / s, 'above') if s > 0 else (np.nan, 'nan')
def pair_theta(A, B, et, em):
    """A 를 '뜨거운 쪽', B 를 기준(320) 으로."""
    hot = [float(np.nanmean(lg(A[T][et]))) for T in TS]; r = float(np.nanmean(lg(B['320'][em])))
    return theta(hot, r)
def sfsum(v, H, n=2000):
    v = pd.Series(np.asarray(v, float)); ok = v.notna().values; s = v[ok].groupby(np.asarray(H)[ok]).median()
    bs = np.array([np.median(s.values[rng.integers(0, len(s), len(s))]) for _ in range(n)])
    return dict(n=int(ok.sum()), n_sf=int(len(s)), sf_median=round(float(s.median()), 2), CI95=[round(float(x), 2) for x in np.percentile(bs, [2.5, 97.5])], frac_pos_sf=round(float((s > 0).mean()), 3))

# ---- (a) 재현 + (c) 접힌 판: H-T 쌍 ----
HT = pd.read_csv(ROOT + '/results/ht_shared_spread.tsv', sep='\t')[['th', 'me', 'H']]
old = pd.read_csv(ROOT + '/results/ht_theta_star.tsv', sep='\t') if os.path.exists(ROOT + '/results/ht_theta_star.tsv') else None
cacheF, cacheD = {}, {}
LF = lambda d: cacheF.setdefault(d, load_full(d)); LD = lambda d: cacheD.setdefault(d, load_fold(d))
rows = []
for th, me, H in HT.itertuples(index=False):
    A, B = LF(th), LF(me)
    r = dict(th=th, me=me, H=H)
    if A is not None and B is not None:
        et, em, tml, tms = shared(A['ca'], B['ca'])
        if len(et) >= 20:
            r['theta_full'], r['kind_full'] = pair_theta(A, B, et, em); r['tm_long'] = tml
    Af, Bf = LD(th), LD(me)
    if Af is not None and Bf is not None:
        et, em, _, _ = shared(Af['ca'], Bf['ca'])
        if len(et) >= 20:
            r['theta_fold'], r['kind_fold'] = pair_theta(Af, Bf, et, em)
            ffs = [Af['ff'][T] for T in TS] + [Bf['ff'][T] for T in TS]; r['minff'] = float(np.nanmin(ffs)) if np.isfinite(ffs).any() else np.nan
    rows.append(r)
R = pd.DataFrame(rows); R.to_csv(OD + '/ht_pairs.tsv', sep='\t', index=False)
OUT['a_reproduce'] = dict(n_pairs=int(R.theta_full.notna().sum()), shift_full=sfsum(R.theta_full - 320, R.H),
    frac_extrap_below=round(float((R.kind_full == 'below').mean()), 3), frac_extrap_above=round(float((R.kind_full == 'above').mean()), 3))
if old is not None:
    c = [x for x in old.columns if 'theta' in x.lower()][:1]
    if c and {'th', 'me'} <= set(old.columns):
        M = R.merge(old[['th', 'me', c[0]]], on=['th', 'me'])
        ok = M.theta_full.notna() & M[c[0]].notna()
        OUT['a_reproduce']['vs_saved'] = dict(column=c[0], n=int(ok.sum()), spearman=round(float(M.theta_full[ok].corr(M[c[0]][ok], method='spearman')), 4),
            median_abs_diff_K=round(float((M.theta_full[ok] - M[c[0]][ok]).abs().median()), 2))
OUT['c_folded'] = dict(all=sfsum(R.get('theta_fold', pd.Series(dtype=float)) - 320, R.H),
    minff_ge_0p8=sfsum((R.theta_fold - 320).where(R.minff >= 0.8), R.H) if 'minff' in R else None,
    n_minff_known=int(R['minff'].notna().sum()) if 'minff' in R else 0,
    frac_extrap_below=round(float((R.get('kind_fold') == 'below').mean()), 3) if 'kind_fold' in R else None)
save(); print('a,c done', json.dumps(OUT, ensure_ascii=False)[:800], flush=True)

# ---- (b) 중온–중온 음성 대조 ----
org = pd.read_csv(ROOT + '/results/pilot_mdcath_organism_temp_v2.tsv', sep='\t')
org = org[org.Topt.notna() & (org.Topt < 40)].copy(); org['H'] = [cath.get(d, h) for d, h in zip(org.dom, org.H)]
htdoms = set(HT.th) | set(HT.me)
tm_min = float(R.tm_long.min()) if 'tm_long' in R and R.tm_long.notna().any() else 0.0
cand = []
for H, g in org.groupby('H'):
    ds = [d for d in g.dom.unique() if os.path.exists(f'{ROOT}/derived/mdcath_spread/{d}.npz')]
    if len(ds) < 2: continue
    prs = [(ds[i], ds[j]) for i in range(len(ds)) for j in range(i + 1, len(ds))]
    rng.shuffle(prs); cand += [(a, b, H) for a, b in prs[:5]]
rows = []
for a, b, H in cand[:1500]:
    A, B = LF(a), LF(b)
    if A is None or B is None: continue
    et, em, tml, _ = shared(A['ca'], B['ca'])
    if len(et) < 20 or tml < tm_min: continue
    tab, kab = pair_theta(A, B, et, em); tba, kba = pair_theta(B, A, em, et)
    o = rng.integers(0, 2)
    rows.append(dict(a=a, b=b, H=H, tm_long=tml, theta_ab=tab, theta_ba=tba, sym=np.nanmean([tab - 320, tba - 320]) if np.isfinite([tab, tba]).any() else np.nan,
                     rand=(tab if o == 0 else tba) - 320, absmean=np.nanmean([abs(tab - 320), abs(tba - 320)]), in_ht=bool(a in htdoms or b in htdoms)))
M = pd.DataFrame(rows); M.to_csv(OD + '/meso_pairs.tsv', sep='\t', index=False)
ht_abs = (R.theta_full - 320).abs()
OUT['b_meso_negative'] = dict(n_pairs=len(M), tm_long_min_matched=round(tm_min, 3), symmetrized=sfsum(M.sym, M.H), random_orientation=sfsum(M.rand, M.H),
    abs_shift_meso=sfsum(M.absmean, M.H), abs_shift_ht=sfsum(ht_abs, R.H),
    frac_pos_pairs_random=round(float((M.rand > 0).mean()), 3), frac_pos_pairs_ht=round(float(((R.theta_full - 320) > 0).mean()), 3),
    note='sym = 두 방향 θ*−320 의 평균(추정량 편향), rand = 방향 무작위 하나, abs = |θ*−320| (상동 쌍 사이의 전형적 크기). H-T 의 +21 K 는 a_reproduce.shift_full.')
save(); print(json.dumps(OUT, indent=1, ensure_ascii=False))
open(ROOT + '/results/logs/ht_theta_meso.done', 'w').write(time.strftime('%F %T'))
