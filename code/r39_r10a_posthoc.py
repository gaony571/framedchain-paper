"""(사후 분석, 2026-10-20 UTC; R10a′ 판정 "다르다" 를 본 뒤 — 사용자 "추천대로 하자. 근데 결과 분석을 제대로 해서"; 계획 = Drive JSON round36_r10a_posthoc_2026-10-20, 돌리기 전)
물음: 주 판 잔차 −0.089 (삽입·결실 쌍의 짝 몫이 상동체 회귀선보다 작음)는 삽입·결실의 효과인가, 시험 표본과 기준 표본의 다른 차(끝 겹침·동일성·TM 축 치우침)인가.
자료: derived/r38/r10a/*.json (주 15, 민감도 32), derived/r38/b6/*.json 의 같은 상과 115 (기준). y = 짝 몫/변 (pair_share_per_edge).
공변량: TM (긴 쪽, 등록), TM_s = TM × max(n, m)/min(n, m) (같은 정렬의 짧은 쪽 정규화; n, m = CATH 길이), dfrac = |nA − nB|/max(nA, nB) (D 에 쓴 변 수),
  mfrac = π̄ 짝 수/min(nA, nB), ident (pilot_mdcath_homolog_pairs.tsv; 시험·기준 같은 표), log n̄.
1 표본 비교: 시험(주) 대 기준의 공변량 중앙값.
2 기준 안 진단: 등록 회귀선(y ~ TM) 잔차와 dfrac·ident·mfrac·log n̄ 의 Spearman (상과 부트스트랩 95%).
3 통제 모형(기준에서 최소제곱 → 시험 잔차 상과 중앙값): M0 TM (등록) · M1 TM + dfrac · M2 TM + ident · M3 TM + ident + dfrac · M4 TM_s · M5 TM_s + ident · M6 TM_s + ident + dfrac.
  구간은 둘: (가) 등록과 같은 고정 회귀선 + 시험 상과 부트스트랩, (나) 결합 부트스트랩(기준 상과·시험 상과를 따로 복원 추출하고 회귀선을 다시 맞춤). 판정 규칙은 등록 TOST 와 같음(여유 = 0.5 × 그 모형 잔차 MAD).
4 모형 없는 짝 비교: 시험 쌍마다 표준화한 (TM, ident, dfrac) 공간에서 가장 가까운 기준 5 쌍의 y 평균과의 차; (TM_s, ident) 판도. 시험 상과 중앙값, 시험 상과 부트스트랩 95%.
5 범위·영향: 기준 TM 범위 밖 시험 쌍을 뺀 M0·M3; 시험 상과 하나씩 뺀 M0·M3 추정값의 범위; 민감도 32 쌍의 M3.
6 빈칸 몫(서술): gap_share_per_edge ~ TM (+ dfrac) 의 시험 잔차.
읽는 규칙(돌리기 전에 정함): M3·M6 와 짝 비교에서 결합 구간이 0 을 포함하거나 판정이 "가를 수 없다/동등" 이면 "M0 의 차는 표본 차(끝 겹침·동일성)로 설명되고 삽입·결실 고유의 차는 검출되지 않는다(MDE 함께)".
  통제 뒤에도 0 을 빼고 여유 밖이면 "통제 뒤에도 남는다". 어느 쪽이든 등록 판정("다르다")은 바꾸지 않는다.
사용: python r39_r10a_posthoc.py <ROOT> 출력: results/r39_r10a_posthoc.json"""
import os, sys, json, glob, time, shutil
import numpy as np
from scipy.stats import spearmanr
ROOT = sys.argv[1]; RES = ROOT + '/results'; NB = 2000; SEED = 20261080


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def load():
    import pandas as pd
    hp = pd.read_csv(RES + '/pilot_mdcath_homolog_pairs.tsv', sep='\t'); idn = {}
    for t in hp.itertuples(): idn[(t.a, t.b)] = idn[(t.b, t.a)] = float(t.ident)
    R10 = [json.load(open(f)) for f in sorted(glob.glob(ROOT + '/derived/r38/r10a/*.json'))]; R10 = [r for r in R10 if 'skip' not in r]
    B6 = [json.load(open(f)) for f in glob.glob(ROOT + '/derived/r38/b6/*.json')]; hb = [r for r in B6 if 'skip' not in r and r['level'] == 'H']
    for r in R10 + hb:
        r['ident_hp'] = idn.get((r['a'], r['b']), np.nan); r['TM_s'] = r['TM'] * max(r['n'], r['m']) / min(r['n'], r['m'])
        r['dfrac'] = abs(r['nA'] - r['nB']) / max(r['nA'], r['nB']); r['mfrac'] = r['n_pibar'] / min(r['nA'], r['nB']); r['lognbar'] = float(np.log(r['nbar']))
    return R10, hb


def X(rows, ks):
    return np.column_stack([np.ones(len(rows))] + [np.array([r[k] for r in rows], float) for k in ks])


def Y(rows, yk='pair_share_per_edge'):
    return np.array([r[yk] for r in rows], float)


def fit(ref, ks, yk='pair_share_per_edge'):
    be = np.linalg.lstsq(X(ref, ks), Y(ref, yk), rcond=None)[0]; res = Y(ref, yk) - X(ref, ks) @ be; return be, float(np.median(np.abs(res - np.median(res))))


def sfmed(e, H):
    u = sorted(set(H)); return float(np.median([np.median(e[H == h]) for h in u]))


def groups(rows):
    H = np.array([r['Ha'] for r in rows]); u = sorted(set(H)); return H, u, [np.where(H == h)[0] for h in u]


def fin(rows, ks):
    return [r for r in rows if all(np.isfinite(r[k]) for k in ks)]


def model(test, ref, ks, seed, yk='pair_share_per_edge'):
    test, ref = fin(test, ks), fin(ref, ks); be, mad = fit(ref, ks, yk); margin = 0.5 * mad; e = Y(test, yk) - X(test, ks) @ be; Ht, ut, gt = groups(test); est = sfmed(e, Ht)
    rng = np.random.default_rng(seed); fixed, joint = [], []; Hr, ur, gr = groups(ref)
    for _ in range(NB):
        it = np.concatenate([gt[i] for i in rng.integers(0, len(ut), len(ut))]); fixed.append(sfmed(e[it], Ht[it]))
        ir = np.concatenate([gr[i] for i in rng.integers(0, len(ur), len(ur))]); rr = [ref[i] for i in ir]
        b2 = np.linalg.lstsq(X(rr, ks), Y(rr, yk), rcond=None)[0]; tt = [test[i] for i in it]; e2 = Y(tt, yk) - X(tt, ks) @ b2; joint.append(sfmed(e2, Ht[it]))
    def verdict(bs):
        lo90, hi90 = np.percentile(bs, [5, 95]); lo95, hi95 = np.percentile(bs, [2.5, 97.5])
        v = '동등' if (-margin < lo90 and hi90 < margin) else ('다르다' if ((lo95 > 0 or hi95 < 0) and abs(est) > margin) else '가를 수 없다')
        return dict(ci90=[float(lo90), float(hi90)], ci95=[float(lo95), float(hi95)], se=float(np.std(bs)), mde80=float(2.8 * np.std(bs)), verdict=v)
    return dict(covariates=list(ks), n_test=len(test), n_sf_test=len(ut), n_ref=len(ref), coef=[float(x) for x in be], mad=mad, margin=margin, est=est, fixed_line=verdict(fixed), joint=verdict(joint))


def sboot_spearman(ref, xk, e, seed):
    Hr, ur, gr = groups(ref); x = np.array([r[xk] for r in ref], float); ok = np.isfinite(x); rng = np.random.default_rng(seed); bs = []
    for _ in range(NB):
        i = np.concatenate([gr[j] for j in rng.integers(0, len(ur), len(ur))]); i = i[ok[i]]; bs.append(spearmanr(x[i], e[i]).correlation)
    return dict(rho=float(spearmanr(x[ok], e[ok]).correlation), ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])


def nn_match(test, ref, ks, k, seed):
    test, ref = fin(test, ks), fin(ref, ks); Z = np.array([[r[c] for c in ks] for r in ref], float); mu, sd = Z.mean(0), Z.std(0); Zr = (Z - mu) / sd; yr = Y(ref)
    d = []
    for r in test:
        z = (np.array([r[c] for c in ks], float) - mu) / sd; nn = np.argsort(((Zr - z) ** 2).sum(1))[:k]; d.append(r['pair_share_per_edge'] - yr[nn].mean())
    d = np.array(d); Ht, ut, gt = groups(test); est = sfmed(d, Ht); rng = np.random.default_rng(seed); bs = []
    for _ in range(NB):
        it = np.concatenate([gt[i] for i in rng.integers(0, len(ut), len(ut))]); bs.append(sfmed(d[it], Ht[it]))
    return dict(covariates=list(ks), k=k, est=est, ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], per_pair=[[r['a'], r['b'], round(float(x), 4)] for r, x in zip(test, d)])


def main():
    R10, hb = load(); main_ = [r for r in R10 if r.get('main')]; out = dict(written=time.strftime('%F %T'), n_main=len(main_), n_sens=len(R10), n_ref=len(hb), posthoc=True)
    ks_all = ('TM', 'TM_s', 'ident_hp', 'dfrac', 'mfrac', 'nbar')
    out['covariates'] = {k: dict(test_main=float(np.nanmedian([r[k] for r in main_])), test_sens=float(np.nanmedian([r[k] for r in R10])), ref=float(np.nanmedian([r[k] for r in hb])),
                                 ref_range=[float(np.nanmin([r[k] for r in hb])), float(np.nanmax([r[k] for r in hb]))], test_range=[float(np.nanmin([r[k] for r in main_])), float(np.nanmax([r[k] for r in main_]))]) for k in ks_all}
    be0, _ = fit(hb, ('TM',)); e0 = Y(hb) - X(hb, ('TM',)) @ be0
    out['ref_residual_vs'] = {k: sboot_spearman(hb, k, e0, SEED + i) for i, k in enumerate(('dfrac', 'ident_hp', 'mfrac', 'lognbar', 'TM_s'))}
    M = {'M0 TM (등록)': ('TM',), 'M1 TM + dfrac': ('TM', 'dfrac'), 'M2 TM + ident': ('TM', 'ident_hp'), 'M3 TM + ident + dfrac': ('TM', 'ident_hp', 'dfrac'),
         'M4 TM_s': ('TM_s',), 'M5 TM_s + ident': ('TM_s', 'ident_hp'), 'M6 TM_s + ident + dfrac': ('TM_s', 'ident_hp', 'dfrac')}
    out['models_main'] = {nm: model(main_, hb, ks, SEED + 10 + i) for i, (nm, ks) in enumerate(M.items())}
    out['models_sens32'] = {nm: model(R10, hb, ks, SEED + 30 + i) for i, (nm, ks) in enumerate(M.items()) if nm.startswith(('M0', 'M3', 'M6'))}
    out['nn_match'] = dict(TM_ident_dfrac=nn_match(main_, hb, ('TM', 'ident_hp', 'dfrac'), 5, SEED + 50), TMs_ident=nn_match(main_, hb, ('TM_s', 'ident_hp'), 5, SEED + 51),
                           TM_only=nn_match(main_, hb, ('TM',), 5, SEED + 52))
    lo, hi = min(r['TM'] for r in hb), max(r['TM'] for r in hb); inr = [r for r in main_ if lo <= r['TM'] <= hi]
    out['in_range'] = dict(n=len(inr), dropped=[[r['a'], r['b'], r['TM']] for r in main_ if not (lo <= r['TM'] <= hi)], M0=model(inr, hb, ('TM',), SEED + 60), M3=model(inr, hb, ('TM', 'ident_hp', 'dfrac'), SEED + 61))
    loso = {}
    for nm, ks in (('M0', ('TM',)), ('M3', ('TM', 'ident_hp', 'dfrac'))):
        mt = fin(main_, ks); be, _ = fit(fin(hb, ks), ks); e = Y(mt) - X(mt, ks) @ be; H = np.array([r['Ha'] for r in mt]); vals = {h: sfmed(e[H != h], H[H != h]) for h in sorted(set(H))}
        loso[nm] = dict(min=float(min(vals.values())), max=float(max(vals.values())), by_sf={h: round(v, 4) for h, v in vals.items()})
    out['leave_one_sf_out'] = loso
    out['per_pair'] = [dict(ab=[r['a'], r['b']], sf=r['Ha'], TM=round(r['TM'], 3), TM_s=round(r['TM_s'], 3), ident=round(r['ident_hp'], 3), dfrac=round(r['dfrac'], 3), mfrac=round(r['mfrac'], 3),
                            y=round(r['pair_share_per_edge'], 4), res_M0=round(float(r['pair_share_per_edge'] - X([r], ('TM',)) @ be0), 4)) for r in main_]
    out['gap_share'] = dict(M0=model(main_, hb, ('TM',), SEED + 70, 'gap_share_per_edge'), M1=model(main_, hb, ('TM', 'dfrac'), SEED + 71, 'gap_share_per_edge'))
    out['rule'] = '사후 분석(판정 바꾸지 않음). 읽는 규칙은 docstring 과 Drive JSON round36_r10a_posthoc_2026-10-20 그대로.'
    p = RES + '/r39_r10a_posthoc.json'; json.dump(out, open(p + '.tmp', 'w'), ensure_ascii=False, indent=1); os.replace(p + '.tmp', p)
    b = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + 'r39_r10a_posthoc.json')
    log('R10A_POSTHOC_DONE')


if __name__ == '__main__':
    main()
