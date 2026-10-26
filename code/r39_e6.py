"""(round39 E6 위치 판정; 2026-10-21 UTC, Drive JSON round39_amend1_2026-10-21 — r39_unfold_D 값을 보기 전에 짬) 위 닻 띠와 그림 5 칸의 위치.
규칙(동결 문서 E6 그대로): 기본값 = 위 닻 없음. 아래를 모두 통과할 때만 "문턱으로 정의한 띠" 로 그린다.
  띠 = mdCATH 413·450 K 접힘 대 풀림(문턱 (0.7, 0.3)·(0.75, 0.4)·(0.8, 0.5)) Δ²_D/√n̄ 의 상과 중앙값(상과 안 중앙값의 중앙값), 상과 복원 2,000 의 95% (seed 20261079). 띠 여섯.
  칸 = 그림 5 의 D 칸, 값 = 단위(도메인·사슬·쌍) Δ²_D/√n̄ 의 상과 중앙값(점):
    온도 걸음 348·379·413·450 (줄기 336 D, 뒤 절반), R3′ (양쪽 50–100 ns, r41_dff50 = 그림 5 주 판; 전체 창 r38_dff_D 는 함께 적되 판정 밖),
    7c aSAMt320·AlphaFlow (P2 67, compare_D2, 전체 창 Q 0.75), 상동체 연속 축(B6, H 수준) TM 구간별, 자연 삽입·결실(R10a′ 주 판).
  (1) 위치(아래 < 띠 하한 ≤ 안 ≤ 띠 상한 < 위)가 여섯 띠에서 모두 같다(세 문턱 × 413·450).
  (2) 정규화 대안: 칸 값 = 단위 Δ²_D/n̄^α̂ 의 상과 중앙값, α̂ = 그 띠의 보조 집합에서 E1 과 같은 IRLS (r39_core.irls). 띠마다 칸 순위의 Kendall τ(Δ²/√n̄ 대 Δ²/n̄^α̂) ≥ 0.9.
  셋 다 통과 → "문턱으로 정의한 띠", 아니면 "위 닻 없음(Δ²_D 자체)". 배수는 쓰지 않는다. Majewski 는 이 스크립트 밖(외부 확인).
사용: python r39_e6.py <ROOT>   (r39_D_unfold.py run 뒤)   출력: results/r39_E6.json"""
import os, sys, json, glob, time, shutil, warnings
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import kendalltau
warnings.filterwarnings('ignore', category=RuntimeWarning)
ROOT = sys.argv[1]; RES = ROOT + '/results'; DER = ROOT + '/derived'; NB = 2000; SEED = 20261079
THR = ('0.7-0.3', '0.75-0.4', '0.8-0.5'); TU = ('413', '450')
_s = open(ROOT + '/code/r39_core.py').read()
exec(_s[_s.index('def log(*a):'):_s.index('# ---------- 자료 ----------')])     # log, jdump, cath_sf, sf_groups, sf_meds, boot_joint, ci, summ (그대로)
exec(_s[_s.index('def _c(n, y, w, a):'):_s.index('def wmedian(')])               # _c, _rss, _rss_fe, nls, irls (그대로)


def units_unfold(T, thr):
    out = []
    for f in sorted(glob.glob(DER + '/r39/unfold_D/*.json')):
        j = json.load(open(f)); c = j.get('cells', {}).get(f'{T}|{thr}', {})
        if 'D' in c: out.append(dict(u=j['dom'], dom=j['dom'], D2=c['D']['D2'], n=float(j['ne'])))
    return out


def cells():
    C = {}
    J = [json.load(open(f)) for f in sorted(glob.glob(DER + '/r39/D/*.json'))]
    for T in ('348', '379', '413', '450'):
        C[f'온도 +{dict(zip(("348", "379", "413", "450"), (28, 59, 93, 130)))[T]} K (뒤 절반)'] = [dict(u=j['dom'], dom=j['dom'], D2=j['cells'][f'MD{T}-MD320|half']['D']['D2'], n=float(j['ne']))
                                                                                         for j in J if f'MD{T}-MD320|half' in j.get('cells', {})]
    for tag, sub in (('R3′ 역장+규약+20 K (50–100 ns, 주)', 'r41/dff50'), ('R3′ (전체 창, 참고)', 'r38/dff_D')):
        C[tag] = [dict(u=j['chain'], dom=j['dom'], D2=j['ff']['D2'], n=float(j['ne'] if 'ne' in j else j['L'] - 1)) for j in (json.load(open(f)) for f in sorted(glob.glob(f'{DER}/{sub}/*.json'))) if 'ff' in j]
    doms67 = open(DER + '/r38/p2_67.txt').read().split()
    for tag, key in (('7c aSAMt 320 K', 'aSAMt320-MD320|full|0.75'), ('7c AlphaFlow', 'AlphaFlow-MD320|full|0.75')):
        rows = []
        for d in doms67:
            f = f'{DER}/r38/compare_D2/{d}.json'
            if not os.path.exists(f): continue
            c = json.load(open(f)).get('cells', {}).get(key)
            if c and c.get('d2n'): rows.append(dict(u=d, dom=d, D2=c['D2'], n=float((c['D2'] / c['d2n']) ** 2)))
        C[tag] = rows
    B6 = [json.load(open(f)) for f in glob.glob(DER + '/r38/b6/*.json')]; B6 = [j for j in B6 if 'skip' not in j and j.get('D') and j.get('level') == 'H']
    for b in sorted({j['bin'] for j in B6}):
        C[f'상동체 TM {b}'] = [dict(u=f"{j['a']}|{j['b']}", dom=j['a'], sf=j['Ha'], D2=j['D']['D2'], n=float(j['nbar'])) for j in B6 if j['bin'] == b]
    R10 = [json.load(open(f)) for f in glob.glob(DER + '/r38/r10a/*.json')]
    C['자연 삽입·결실 (R10a′ 주)'] = [dict(u=f"{j['a']}|{j['b']}", dom=j['a'], sf=j['Ha'], D2=j['D']['D2'], n=float(j['nbar'])) for j in R10 if 'skip' not in j and j.get('main') and j.get('D')]
    return {k: v for k, v in C.items() if v}


def sfmed(rows, sfmap, alpha=None):
    v = np.array([r['D2'] / (r['n'] ** (0.5 if alpha is None else alpha)) for r in rows]); u, g = sf_groups([r.get('sf') or sfmap.get(r['dom'], r['dom']) for r in rows])
    return float(np.nanmedian(sf_meds(v, g, len(u)))), len(u)


def main():
    sfmap = cath_sf(); rng = np.random.default_rng(SEED); CL = cells(); out = dict(written=None, bands={}, cells={}, rule=__doc__.split('\n사용')[0])
    for k, rows in CL.items():
        est, ns = sfmed(rows, sfmap); out['cells'][k] = dict(n_units=len(rows), n_sf=ns, d2n_sf_median=est, nbar_median=float(np.median([r['n'] for r in rows])))
    pos = {}; taus = {}
    for T in TU:
        for thr in THR:
            key = f'{T}|{thr}'; U = units_unfold(T, thr)
            if len(U) < 5: out['bands'][key] = dict(n=len(U)); continue
            v = np.array([r['D2'] / np.sqrt(r['n']) for r in U]); u, g = sf_groups([sfmap.get(r['dom'], r['dom']) for r in U]); pt, bs = boot_joint([v], g, len(u), rng)
            lo, hi = ci(bs[0]); a, gam, it = irls(np.array([r['n'] for r in U]), np.array([r['D2'] for r in U]))
            out['bands'][key] = dict(n=len(U), n_sf=len(u), d2n_sf_median=float(pt[0]), ci95=[lo, hi], alpha_irls=float(a), gamma=float(gam))
            pos[key] = {k: ('아래' if c['d2n_sf_median'] < lo else '위' if c['d2n_sf_median'] > hi else '안') for k, c in out['cells'].items()}
            main_k = [k for k in out['cells'] if '참고' not in k]; x1 = [out['cells'][k]['d2n_sf_median'] for k in main_k]; x2 = [sfmed(CL[k], sfmap, a)[0] for k in main_k]
            taus[key] = float(kendalltau(x1, x2).correlation)
    main_k = [k for k in out['cells'] if '참고' not in k]; keys = list(pos)
    stable = {k: len({pos[b][k] for b in keys}) == 1 for k in main_k} if keys else {}
    c1 = bool(keys) and len(keys) == 6 and all(stable.values()); c3 = bool(taus) and min(taus.values()) >= 0.9
    out.update(positions=pos, position_stable_by_cell=stable, kendall_tau_by_band=taus, crit_position_all_six=c1, crit_tau_min_ge_0_9=c3,
               verdict='문턱으로 정의한 띠' if (c1 and c3) else '위 닻 없음 (Δ²_D 자체)', n_bands=len(keys), written=time.strftime('%F %T'))
    jdump(out, RES + '/r39_E6.json'); log('R39_E6', out['verdict'], '위치 안정', c1, 'τ 최소', min(taus.values()) if taus else None); log('R39_E6_DONE')


if __name__ == '__main__':
    main()
