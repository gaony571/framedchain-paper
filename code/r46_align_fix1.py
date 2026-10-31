"""(고침 = Drive results/prereg_2026-10-30_addendum9_fix1.json; 원 등록 addendum9 "r46_align") 순위를 평균 순위로 바꿔 r46_align 의 통계를 다시 낸다.
r46_align.py (a2c0156b) 의 함수를 그대로 불러 쓰고, 순위 함수만 scipy.stats.rankdata(method='average') 로 바꾼다(analyse_sf·boot 가 쓰는 전역 ranks).
1단계(derived/r46/align_pairs.jsonl, tm_check.json pass = true)는 그대로 쓴다. 판정 규칙·문장 갈래·seed·재표집 수는 원 등록 그대로.
점검 2a: argsort 판(원래 ranks)으로 전체 C 기준 판의 n·n_sf 가 r38_c.json 과 같고 점추정 |차| < 1e-3 (C-1, C-2).
점검 2b: 평균 순위 판을 등록 seed(20261071, 20261074)로 내어 est·lo·hi 가 고침 기록에 적은 값과 |차| < 1e-9.
서술: argsort 판의 기준·A·B0·B2 (같은 seed 흐름; 판정에 쓰지 않음).
사용: python r46_align_fix1.py <ROOT>   출력: results/r46_align_fix1.json"""
import sys, json, time, os
import numpy as np
ROOT = sys.argv[1]
_src = open(ROOT + '/code/r46_align.py').read(); _src = _src[:_src.index("if __name__ == '__main__':")]
sys.argv = ['r46_align', ROOT, 'none']; exec(_src)   # r46_align 의 함수(boot, ci, verdict, load_fp, corr2_rows, ORDER, …)와 그 안의 r38_c·r38_p3 함수
from scipy.stats import rankdata
ranks_argsort = ranks
def ranks_avg(v): return rankdata(v, method='average').astype(float)
REC = dict(C1=dict(est=0.22460683061158487, lo=0.11919817035584895, hi=0.3326286464040772),
           C2=dict(est=0.18094574136027425, lo=0.0786982706674338, hi=0.28307810014941687))
FO = RES + '/r46_align_fix1.json'


def with_ranks(rk, fn):
    g = globals(); old = g['ranks']; g['ranks'] = rk
    try: return fn()
    finally: g['ranks'] = old


def main():
    tmc = json.load(open(FC)) if os.path.exists(FC) else {}
    if not tmc.get('pass'): log('R46F STOP: tm_check 통과 파일 없음'); return
    base = merged(OUT, OUT + '/anm_pairs.jsonl'); reg = json.load(open(RES + '/r38_c.json'))
    q1, q2 = with_ranks(ranks_argsort, lambda: (analyse_sf(base, 'main', z=('control',), seed=20261071), analyse_sf(base, 'main', z=('control', 'r_anm'), seed=20261074)))
    ok_a = all(q[k] == reg[c]['main'][k] for q, c in ((q1, 'C1'), (q2, 'C2')) for k in ('n', 'n_sf')) and abs(q1['est'] - reg['C1']['main']['est']) < 1e-3 and abs(q2['est'] - reg['C2']['main']['est']) < 1e-3
    a1, a2 = with_ranks(ranks_avg, lambda: (analyse_sf(base, 'main', z=('control',), seed=20261071), analyse_sf(base, 'main', z=('control', 'r_anm'), seed=20261074)))
    dev_b = max(abs(a[k] - REC[c][k]) for a, c in ((a1, 'C1'), (a2, 'C2')) for k in ('est', 'lo', 'hi'))
    out = dict(written=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()),
               check=dict(tm_check=tmc, check_2a=dict(ok=bool(ok_a), argsort_C1=q1, argsort_C2=q2, registered_C1=reg['C1']['main'], registered_C2=reg['C2']['main']),
                          check_2b=dict(ok=bool(dev_b < 1e-9), dev=float(dev_b), avg_C1=a1, avg_C2=a2)))
    if not (ok_a and dev_b < 1e-9):
        out['stopped'] = '점검 2 실패 — 판정을 내지 않음'; json.dump(out, open(FO, 'x'), ensure_ascii=False, indent=1); log('R46F STOP', ok_a, dev_b); return
    AL, c2r = load_fp(), corr2_rows(); rows = []
    for r in base:
        m = r['main']; key = tuple(r['ab']); a = AL.get(key)
        if a is None or key not in c2r or not all(m.get(c) is not None for c in ('r_dyn', 'ident', 'control', 'r_anm')): continue
        if a['n_res'] < MIN_RES or a['prec0'] is None: continue
        cv = c2r[key].get('cover') or {}
        rows.append(dict(ab=r['ab'], H=r['H'], main=dict(m, prec0=a['prec0'], prec2=a['prec2'], rec0=cv.get('tm_med_by_pibar|t0'), rec2=cv.get('tm_med_by_pibar|t2')),
                         alt=dict(r_dyn=m['r_dyn'], control=m['control'], r_anm=m['r_anm'], ident=a['ident_tm_res']),
                         ap=dict(r_dyn=m['r_dyn'], control=m['control'], r_anm=m['r_anm'], ident=a['ident_aprime'], n=a['n_aprime'])))
    Z1, Z2 = ('control',), ('control', 'r_anm')
    specs = [('base_C1', 'main', 'r_dyn', 'ident', Z1), ('base_C2', 'main', 'r_dyn', 'ident', Z2),
             ('A_C1', 'alt', 'r_dyn', 'ident', Z1), ('A_C2', 'alt', 'r_dyn', 'ident', Z2),
             ('B0_C1', 'main', 'r_dyn', 'ident', Z1 + ('prec0',)), ('B0_C2', 'main', 'r_dyn', 'ident', Z2 + ('prec0',)),
             ('B2_C1', 'main', 'r_dyn', 'ident', Z1 + ('prec2',)), ('B2_C2', 'main', 'r_dyn', 'ident', Z2 + ('prec2',)),
             ('cpl_rdyn_prec0', 'main', 'r_dyn', 'prec0', Z2), ('cpl_ident_prec0', 'main', 'ident', 'prec0', Z2),
             ('cpl_rdyn_prec2', 'main', 'r_dyn', 'prec2', Z2), ('cpl_ident_prec2', 'main', 'ident', 'prec2', Z2)]
    est, bs, nsf = with_ranks(ranks_avg, lambda: boot(rows, specs, SEED, NBOOT))
    res = {nm: dict(est=est[nm], lo=ci(bs[nm])[0], hi=ci(bs[nm])[1]) for nm in est}
    diff = {f'base-{v}_{c}': dict(est=est[f'base_{c}'] - est[f'{v}_{c}'], lo=ci(bs[f'base_{c}'] - bs[f'{v}_{c}'])[0], hi=ci(bs[f'base_{c}'] - bs[f'{v}_{c}'])[1])
            for v in ('A', 'B0', 'B2') for c in ('C1', 'C2')}
    est2, bs2, _ = with_ranks(ranks_avg, lambda: boot(rows, specs[:8], SEED2, NBOOT))
    seed2 = {nm: dict(lo=ci(bs2[nm])[0], hi=ci(bs2[nm])[1]) for nm in est2}
    estq, bsq, _ = with_ranks(ranks_argsort, lambda: boot(rows, specs[:8], SEED, NBOOT))
    argsort_desc = {nm: dict(est=estq[nm], lo=ci(bsq[nm])[0], hi=ci(bsq[nm])[1]) for nm in estq}
    rows_ap = [r for r in rows if r['ap']['n'] >= MIN_RES and r['ap']['ident'] is not None]
    rows_rc = [r for r in rows if r['main']['rec0'] is not None and r['main']['rec2'] is not None]
    desc_an = with_ranks(ranks_avg, lambda: dict(
        Aprime_C1=analyse_sf(rows_ap, 'ap', z=Z1, seed=SEED + 10, B=NBOOT), Aprime_C2=analyse_sf(rows_ap, 'ap', z=Z2, seed=SEED + 11, B=NBOOT),
        recall0_C1=analyse_sf(rows_rc, 'main', z=Z1 + ('rec0',), seed=SEED + 12, B=NBOOT), recall0_C2=analyse_sf(rows_rc, 'main', z=Z2 + ('rec0',), seed=SEED + 13, B=NBOOT),
        recall2_C1=analyse_sf(rows_rc, 'main', z=Z1 + ('rec2',), seed=SEED + 14, B=NBOOT), recall2_C2=analyse_sf(rows_rc, 'main', z=Z2 + ('rec2',), seed=SEED + 15, B=NBOOT)))
    med = lambda k, key='main': float(np.median([r[key][k] for r in rows]))
    desc = dict(n_pairs=len(rows), n_sf=nsf, n_aprime=len(rows_ap), n_recall=len(rows_rc),
                ident_pibar_median=med('ident'), ident_tm_res_median=med('ident', 'alt'), prec0_median=med('prec0'), prec2_median=med('prec2'),
                spearman_ident_pibar_tm=float(spearman(np.array([r['main']['ident'] for r in rows]), np.array([r['alt']['ident'] for r in rows]))))
    vA = verdict(res['A_C1'], res['A_C2']); vB0 = verdict(res['B0_C1'], res['B0_C2']); vB2 = verdict(res['B2_C1'], res['B2_C2'])
    vB = min((vB0, vB2), key=ORDER.get)
    reading = 'both_hold' if vA == vB == 'holds' else 'both_fail' if vA == vB == 'fails' else 'mixed'
    not_all = any(d['lo'] > 0 for d in diff.values())   # 짝지은 차 여섯 모두가 대상
    if reading == 'both_hold': case = 'hold_not_all' if not_all else 'hold'
    elif reading == 'both_fail': case = 'fail'
    elif vB == 'fails': case = 'A_ok_B_fail'
    elif vA == 'fails': case = 'A_fail_B_ok'
    else: case = 'mixed_partly'
    out.update(desc=desc, est=res, paired_diff=diff, seed2=seed2, argsort_descriptive=argsort_desc, descriptive=desc_an,
               verdict=dict(A=vA, B0=vB0, B2=vB2, B=vB, reading=reading, not_all=bool(not_all), sentence_case=case, abstract_change=bool(vA == 'fails')),
               rule='고침 addendum9_fix1: 모든 순위를 평균 순위로. 판정 규칙은 원 등록 r46_align 그대로(판마다 C-1·C-2 하한 > 0 → holds, C-1 만 → partly, 아니면 fails; B = B0·B2 중 약한 쪽; '
                    '문장 갈래 hold / hold_not_all (짝지은 차 여섯 중 하나라도 하한 > 0) / fail / A_ok_B_fail / A_fail_B_ok / mixed_partly; 초록은 A 가 fails 일 때만 단서). seed 20261121 10,000 회 한 흐름. argsort 판은 서술.')
    json.dump(out, open(FO, 'x'), ensure_ascii=False, indent=1)
    log('R46F', out['verdict']); log('R46F_DONE')


if __name__ == '__main__':
    try: main()
    except Exception as e: traceback.print_exc(); log('ERR r46f', repr(e)[:300])
