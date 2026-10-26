
"""접힌 프레임 분석 (2026-10-17; 예측은 results/prereg_2026-10-17.json). 입력 derived/mdcath_5T_folded, derived/mdcath_5T (Q).
(i) C2 Q 점검 (ii) C2 깨끗한 검정 + V2 접힌 판 (iii) H-T 외삽 없는 비교 (iv) C3 스칼라 대조."""
import sys, os, json, time, glob
ROOT = sys.argv[1]
import numpy as np, pandas as pd
from tmtools import tm_align
FD, RAW5 = ROOT + '/derived/mdcath_5T_folded', ROOT + '/derived/mdcath_5T'; rng = np.random.default_rng(31); OUT = {}
lg = lambda x: np.log(np.maximum(x, 1e-8))
def save(): json.dump(OUT, open(ROOT + '/results/folded_analyses.json', 'w'), indent=1, ensure_ascii=False)
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
# ---- (i) Q 점검 ----
def qstats(d, T):
    z = np.load(f'{RAW5}/{d}.npz'); qs = np.concatenate([z[k].astype(float) for k in z.files if k.startswith(T + '_') and k.endswith('_Q')])
    return float(qs.mean()), float(np.quantile(qs, 0.1)), float((qs >= 0.75).mean())
rows = []
for th, me, H in HT.itertuples(index=False):
    if th not in D or me not in D: continue
    row = dict(th=th, me=me, H=H)
    for T in ('320', '348', '379'):
        a, b = qstats(th, T), qstats(me, T)
        row[f'dQmean_{T}'], row[f'dQ10_{T}'], row[f'dfold_{T}'] = a[0] - b[0], a[1] - b[1], a[2] - b[2]
        row[f'me_fold_{T}'], row[f'th_fold_{T}'] = b[2], a[2]
    rows.append(row)
Qd = pd.DataFrame(rows)
OUT['i_Q_check'] = {c: sfsum(Qd[c], Qd.H) for c in Qd.columns if c.startswith(('dQ', 'dfold', 'me_fold', 'th_fold'))}
OUT['i_Q_check']['note'] = '호열 − 중온 (쌍, 상과 중앙값). 예측: 379 K 에서 양(중온이 더 풀림), 348 K 에서 작음.'; save(); print('i done', flush=True)
# ---- (ii) C2 ----
rows = []
for th, me, H in HT.itertuples(index=False):
    if th not in D or me not in D: continue
    a, b = D[th], D[me]
    if not all(f'{T}_spread_rep' in a and f'{T}_spread_rep' in b for T in ('320', '348', '379')): continue
    et, em = shared(a, b)
    if len(et) < 20: continue
    g = lambda z, T, idx, R_=None: lg(mean_rep(z, f'{T}_spread_rep', R_)[idx])
    row = dict(H=H)
    # 깨끗한 판: 적응은 320, 온도 반응은 348→379, 반복 5개 모두
    dad = g(a, '320', et) - g(b, '320', em); lme, lth = g(b, '320', em), g(a, '320', et)
    relf = lambda z, idx: (lambda r: 2 * r / (1 + r))(np.corrcoef(lg(mean_rep(z, '320_spread_rep', [0, 1])[idx]), lg(mean_rep(z, '320_spread_rep', [2, 3, 4])[idx]))[0, 1])
    rel = np.array([relf(b, em), relf(a, et)]); X = np.c_[lme, lth]
    for side, dt in (('me', g(b, '379', em) - g(b, '348', em)), ('th', g(a, '379', et) - g(a, '348', et))):
        ok = np.isfinite(dad) & np.isfinite(dt) & np.isfinite(X).all(1)
        if ok.sum() < 15: continue
        row[f'clean_{side}_raw'] = float(pd.Series(dad[ok]).rank().corr(pd.Series(dt[ok]).rank()))
        # (EIV 판은 쓰지 않음: 공변량이 Δ_adapt 의 성분 그 자체라 회귀하면 Δ_adapt 가 지워진다)
    # V2 접힌 판 (앞 판과 같은 짜임: 표준화 후 차감 + 두 쪽 수준 EIV, 반복 분리)
    A_, B_, C_ = [0, 1], [2, 3], [4]
    dad2 = z_(g(a, '348', et, B_)) - z_(g(b, '348', em, B_)); X2 = np.c_[g(b, '348', em, C_), g(a, '348', et, C_)]
    rel2 = np.array([np.corrcoef(g(b, '348', em, C_), g(b, '348', em, A_))[0, 1], np.corrcoef(g(a, '348', et, C_), g(a, '348', et, A_))[0, 1]])
    for side, dt in (('me', g(b, '379', em, A_) - g(b, '348', em, A_)), ('th', g(a, '379', et, A_) - g(a, '348', et, A_))):
        ok = np.isfinite(dad2) & np.isfinite(dt) & np.isfinite(X2).all(1)
        if ok.sum() >= 15 and np.isfinite(rel2).all(): row[f'V2fold_{side}_eiv'] = float(np.corrcoef(eiv_resid(dad2[ok], X2[ok], rel2), eiv_resid(dt[ok], X2[ok], rel2))[0, 1])
    rows.append(row)
C2 = pd.DataFrame(rows)
for pre in ('clean_{}_raw', 'clean_{}_eiv', 'V2fold_{}_eiv'):
    cm, ct = pre.format('me'), pre.format('th')
    if cm in C2 and ct in C2: C2[pre.format('sym')] = (C2[cm] + C2[ct]) / 2
OUT['ii_C2'] = {c: sfsum(C2[c], C2.H) for c in C2.columns if c != 'H'}
OUT['ii_C2']['note'] = '판정은 대칭(sym) 추정량. 예측: clean_sym ≈ 0, V2fold_th 는 앞 판(−0.135)보다 크게 줄어듦.'; save(); print('ii done', flush=True)
# ---- (iii) H-T 외삽 없는 비교 ----
rows = []
for th, me, H in HT.itertuples(index=False):
    if th not in D or me not in D or not all(f'{T}_spread_rep' in D[x] for T in ('320', '348') for x in (th, me)): continue
    a, b = D[th], D[me]; et, em = shared(a, b)
    if len(et) < 20: continue
    rows.append(dict(H=H, th348_me320=float(np.nanmean(lg(mean_rep(a, '348_spread_rep')[et]) - lg(mean_rep(b, '320_spread_rep')[em]))),
                     th320_me320=float(np.nanmean(lg(mean_rep(a, '320_spread_rep')[et]) - lg(mean_rep(b, '320_spread_rep')[em]))),
                     th348_me348=float(np.nanmean(lg(mean_rep(a, '348_spread_rep')[et]) - lg(mean_rep(b, '348_spread_rep')[em])))))
HTd = pd.DataFrame(rows); OUT['iii_HT_direct'] = {c: sfsum(HTd[c], HTd.H) for c in HTd.columns if c != 'H'}
OUT['iii_HT_direct']['note'] = '공유 변의 평균 로그 비(접힌 프레임). 예측: th348_me320 ≈ 0 (약간 양~0), th320_me320 < 0.'; save(); print('iii done', flush=True)
# ---- (iv) C3 스칼라 대조 ----
doms = sorted(d for d in D if all(f'{T}_spread_rep' in D[d] for T in ('320', '348')))
cnts = {}
for d in doms:
    ca = D[d]['ca_ref'].astype(float); M = np.linalg.norm(ca[:, None] - ca[None], axis=-1); c = (M < 10).sum(1) - 1; cnts[d] = (c[:-1] + c[1:]) / 2
Q3 = np.quantile(np.concatenate(list(cnts.values())), [1 / 3, 2 / 3])
CL = [f'{s}{b}{t}' for s in 'HEL' for b in range(3) for t in range(2)]; ci = {c: k for k, c in enumerate(CL)}
E = {}
for d in doms:
    z = D[d]; s = np.array(['H' if c in 'HGI' else 'E' if c in 'EB' else 'L' for c in z['dssp320']]); ess = np.where(s[:-1] == s[1:], s[:-1], 'L'); ne = len(ess)
    term = np.zeros(ne, int); term[:3] = 1; term[-3:] = 1; cls = np.array([ci[f'{a}{b}{t}'] for a, b, t in zip(ess, np.searchsorted(Q3, cnts[d]), term)])
    prof = {}
    for nm, key, per_res in (('frame', 'spread_rep', False), ('rmsf2', 'rmsf2_rep', True), ('wrmsf2', 'wrmsf2_rep', True)):
        def edge(T, R_):
            v = mean_rep(z, f'{T}_{key}', R_); return 0.5 * (v[:-1] + v[1:]) if per_res else v
        dB = lg(edge('348', [2, 3, 4])) - lg(edge('320', [2, 3, 4])); lev = lg(edge('320', [0, 1]))
        ok = np.isfinite(dB) & np.isfinite(lev)
        if ok.sum() < 20: continue
        dB = np.where(ok, dB, np.nan); dB = dB - np.nanmean(dB)
        bins = np.minimum((pd.Series(lev).rank(pct=True).values * 10).astype(int), 9); bm = np.array([np.nanmean(dB[bins == k]) if (bins == k).any() else 0 for k in range(10)])
        res = dB - bm[bins]; prof[nm] = (dB, res - np.nanmean(res))
    if len(prof) == 3: E[d] = dict(cls=cls, **prof)
dd = sorted(E); sfm = {d: cath.get(d, 'solo:' + d) for d in dd}
reps = pd.Series(dd).groupby([sfm[d] for d in dd]).apply(lambda s: s.iloc[rng.integers(0, len(s))]).values
def fpv(lbl, v, K):
    ok = np.isfinite(v); s = np.bincount(lbl[ok], v[ok], K); n = np.bincount(lbl[ok], minlength=K); return np.where(n >= 2, s / np.maximum(n, 1), np.nan)
def stat(F, keep):
    F = np.where(keep[None], F, np.nan); lo, hi = np.nanpercentile(F, [1, 99]); F = np.clip(F, lo, hi); tot, cnt = np.nansum(F, 0), (~np.isnan(F)).sum(0); o = []
    for k in range(len(F)):
        f = F[k]; ok = ~np.isnan(f); cons = (tot - np.where(ok, f, 0)) / np.maximum(cnt - ok, 1); ok = ok & (cnt - ok > 0)
        if ok.sum() >= 6 and np.std(f[ok]) > 1e-12 and np.std(cons[ok]) > 1e-12: o.append(np.corrcoef(f[ok], cons[ok])[0, 1])
    return (float(np.median(o)) if o else np.nan), np.array(o)
c3 = dict(n_domains=len(dd), n_reps=int(len(reps)))
F0 = np.array([fpv(E[d]['cls'], E[d]['frame'][0], len(CL)) for d in reps]); KEEP = (~np.isnan(F0)).mean(0) >= 0.2
for nm in ('frame', 'rmsf2', 'wrmsf2'):
    for j, lab in ((0, 'raw'), (1, 'level_removed_nonpar')):
        F = np.array([fpv(E[d]['cls'], E[d][nm][j], len(CL)) for d in reps]); o, per = stat(F, KEEP)
        nl = [stat(np.array([fpv(np.random.default_rng(k).permutation(E[d]['cls']), E[d][nm][j], len(CL)) for d in reps]), KEEP)[0] for k in range(300)]
        bs = [np.median(per[rng.integers(0, len(per), len(per))]) for _ in range(2000)]
        c3[f'{nm}_{lab}'] = dict(median_corr=round(o, 4), CI95=[round(float(x), 4) for x in np.percentile(bs, [2.5, 97.5])], null_q99=round(float(np.nanquantile(nl, 0.99)), 4))
c3['within_domain_profile_corr'] = {f'frame_vs_{nm}': round(float(np.nanmedian([pd.Series(E[d]['frame'][0]).corr(pd.Series(E[d][nm][0]), method='spearman') for d in dd])), 4) for nm in ('rmsf2', 'wrmsf2')}
c3['note'] = 'Δ = ln 값(348, 반복 2–4) − ln 값(320, 반복 2–4), 접힌 프레임, 도메인 안 평균 뺌. 수준 통제는 반복 0–1 의 320 수준 십분위 칸. Cα 값은 변 = 양 끝 잔기 평균.'
OUT['iv_C3_scalar'] = c3; OUT['written'] = time.strftime('%F %T'); save(); print('iv done', flush=True)
open(ROOT + '/results/logs/folded_analyses.done', 'w').write(time.strftime('%F %T'))
