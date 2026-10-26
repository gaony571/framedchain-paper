
"""(prereg round29) 치환 행 Δ²/√변 수 (c4_sub 에서), 사슬 잇기 v2 (chaining.tsv 에서). results/c4_ladder_row.json, results/chaining_v2.json."""
import sys, json, glob, time
ROOT = sys.argv[1]
import numpy as np, pandas as pd
rng = np.random.default_rng(29)
def med_ci(v, by=None):
    v = np.asarray(v, float)
    if by is None:
        b = [np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(2000)]
    else:
        g = pd.Series(v).groupby(np.asarray(by)); keys = list(g.groups); arr = {k: g.get_group(k).values for k in keys}; b = []
        for _ in range(2000):
            pick = rng.integers(0, len(keys), len(keys)); b.append(np.median(np.concatenate([arr[keys[i]] for i in pick])))
    return dict(n=int(len(v)), median=round(float(np.median(v)), 4), CI95=[round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])])
# ---- 치환 행 ----
rows = []
for p in glob.glob(ROOT + '/derived/c4_sub/*.json'):
    j = json.load(open(p))
    if 'mutants' not in j or 'floor' not in j: continue
    ne = j['nres'] - 1
    for m in j['mutants']:
        if m.get('ed') is None: continue
        rows.append(dict(wt=j['wt'], ne=ne, floor=j['floor'], ed=m['ed'], ddg=m.get('ddg'), single=m.get('single', True), n_floor_runs=len(j.get('floor_all', []))))
D = pd.DataFrame(rows); D['d2n'] = (D.ed - D.floor) / np.sqrt(D['ne']); D['mult'] = D.ed / D.floor
S = dict(n_wt=int(D.wt.nunique()), n_mut=int(len(D)), note='Δ² = ed(변이체 실행 하나, 야생형 G1 묶음) − 바닥(야생형 실행 하나, 같은 G1 묶음의 중앙값). V_mut = V_WT 가정이면 불편. ddG 부호: 음 = 불안정화.')
for lab, sb in (('single_all', D[D.single]), ('double', D[~D.single])):
    if len(sb) == 0: continue
    S[lab] = dict(d2n=med_ci(sb.d2n, sb.wt), mult=med_ci(sb.mult, sb.wt), frac_d2_pos=round(float((sb.d2n > 0).mean()), 3))
s1 = D[D.single & D.ddg.notna()]
for lab, sb in (('destab_lt_-1', s1[s1.ddg < -1]), ('neutral_abs_le_1', s1[s1.ddg.abs() <= 1]), ('stab_gt_1', s1[s1.ddg > 1])):
    if len(sb): S[lab] = dict(d2n=med_ci(sb.d2n, sb.wt), mult=med_ci(sb.mult, sb.wt))
wt_med = D[D.single].groupby('wt').d2n.median(); S['per_wt_median_of_medians'] = med_ci(wt_med.values)
S['written'] = time.strftime('%F %T'); json.dump(S, open(ROOT + '/results/c4_ladder_row.json', 'w'), indent=1, ensure_ascii=False); print(json.dumps(S, indent=1, ensure_ascii=False))
# ---- 사슬 잇기 v2 ----
C = pd.read_csv(ROOT + '/derived/chaining.tsv', sep='\t'); n = dict(X='nX', Y='nY', XY='nXY')
for u, v in (('X', 'Y'), ('X', 'XY'), ('Y', 'XY')):
    Lb = (C[n[u]] + C[n[v]]) / 2; df2 = C[f'{u}{v}_df'] ** 2
    C[f'{u}{v}_dfn'] = np.sqrt(np.clip(df2 / Lb, 0, None))      # 고정 d₀ 평균 정규화 (증명된 거리)
    C[f'{u}{v}_dfn_lin'] = df2 / Lb                              # 고정 d₀ 선형 1 − S/L̄
    C[f'{u}{v}_dn_lin'] = C[f'{u}{v}_dn'] ** 2                   # d₀(L̄) 선형 1 − S/L̄
def bci(v):
    v = np.asarray(v, float); b = [v[rng.integers(0, len(v), len(v))].mean() for _ in range(2000)]; return [round(float(v.mean()), 3)] + [round(float(x), 3) for x in np.percentile(b, [2.5, 97.5])]
TAU = dict(tms=0.5, dn=np.sqrt(0.5), dfn=np.sqrt(0.5), dn_lin=0.5, dfn_lin=0.5)
O = dict(n=len(C), note='dn = D_TM/√L̄ (d₀(L̄), 근사), dfn = 고정 d₀ (세 짝 공통) D_TM/√L̄ (7.4 따름정리로 증명된 거리), *_lin = 1 − S/L̄ (제곱근 없음). τ 는 TM 0.5 에 맞춤.')
for m, tau in TAU.items():
    xy, xx, yy = C[f'XY_{m}'], C[f'XXY_{m}'], C[f'YXY_{m}']
    o = dict(median_X_Y=round(float(xy.median()), 3), median_X_XY=round(float(xx.median()), 3), triangle_violation=bci(xy > xx + yy + 1e-6), triangle_violation_rel1pct=bci(xy > (xx + yy) * 1.01), single_link_merge=bci((np.maximum(xx, yy) < tau) | (xy < tau)))
    mg = []
    for a, b_, c_ in zip(xy, xx, yy):
        d = {('X', 'Y'): a, ('X', 'XY'): b_, ('Y', 'XY'): c_}; (p1, p2), dmin = min(d.items(), key=lambda kv: kv[1])
        if dmin >= tau: mg.append(False); continue
        rest = ({'X', 'Y', 'XY'} - {p1, p2}).pop(); dm = np.mean([d[tuple(sorted((rest, q), key=['X', 'Y', 'XY'].index))] for q in (p1, p2)])
        mg.append(bool({'X', 'Y'} <= {p1, p2} or dm < tau))
    o['average_link_merge'] = bci(mg); O[m] = o
O['written'] = time.strftime('%F %T'); json.dump(O, open(ROOT + '/results/chaining_v2.json', 'w'), indent=1, ensure_ascii=False); print(json.dumps(O, indent=1, ensure_ascii=False))
