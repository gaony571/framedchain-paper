"""(결과 뒤 진단 — 1차 P1 판정(2026-10-20 12:57 UTC, 서지 않음) 뒤에 짬; 판정을 바꾸지 않음. 사용자 "1단계 하고") r38_p1.py 의 함수를 그대로 불러 쓴다.
부분:
  het : 쌍 사이 실제 흩어짐을 대칭 판(MM)의 LR(320) 에서 잰다: Var(LR320) − 잡음 분산(데워진 쪽 320 K 반복 SE² + 기준 쪽 SE²) → sd_θ = sd_실제 / b.
         그 흩어짐으로 조건 1 (θ 335·350·380) 과 조건 2 (차 15·21·30 K) 의 검정력을 다시 낸다(sim_once·cond2 그대로, 반복 100, 부트스트랩 1,000, seed 20261046).
  win : 뒤 절반 창과 50–100 ns 창의 쌍별 견줌(두 창 모두 θ* 가 난 쌍): 검열 비율, 수치 θ* 차, LR(320) 상관, 반복 잡음 SE.
         같은 쌍 집합에서 두 창의 조건 2 를 다시 내어 쌍 구성 차를 뗀다.
  unc : 검열 없는 짝 통계(계산 전에 Drive JSON 에 정의를 찍음): 쌍마다 LR(320) = log(S_데워진 쪽(320)/S_기준 쪽(320)) (주 칸, 공유 변 GPA).
         상과 요약 = 쌍 LR(320) 의 중앙값. U1 = 호열–중온 상과 요약의 상과 중앙값(음수면 같은 320 K 에서 호열 쪽이 덜 퍼짐),
         U2 = 짝지은 상과 차(호열–중온 − 대칭 판)의 상과 중앙값. 상과 부트스트랩 2,000 (seed 20261045) 95% 구간. K 환산 = −LR / b (b = 0.171/28 K).
사용: python r38_p1_diag.py <ROOT> het,win,unc 출력: results/r38_p1_diag.json"""
import sys, os, json
ROOT = sys.argv[1]; PARTSD = sys.argv[2].split(',')
_s = open(ROOT + '/code/r38_p1.py').read(); _s = _s[:_s.index("ALL1 = [")]
sys.argv = ['r38_p1', ROOT, 'none', '1']; exec(_s)
PJ1 = json.load(open(OUT1 + '/pairs.json')); SP = load_spread(PJ1['domains']); FO = RES + '/r38_p1_diag.json'
OUTD = json.load(open(FO)) if os.path.exists(FO) else {}


def lr320(r, cell=MAIN):
    key = f"{r['a']}__{r['b']}"; H_, R_ = SP.get(r['hot']), SP.get(r['ref'])
    if H_ is None or R_ is None or key not in H_['pairs'] or key not in R_['pairs'] or H_['pairs'][key]['n_shared'] < MINSH1: return None
    k = ckey('320', *cell); hc, rc = H_['pairs'][key]['cells'], R_['pairs'][key]['cells']
    if k not in hc or k not in rc or hc[k]['gpa'] is None or rc[k]['gpa'] is None: return None
    se = lambda runs: float(np.std(np.log(runs), ddof=1) / np.sqrt(len(runs))) if len(runs) > 1 else 0.0
    return dict(H=r['H'], lr=float(np.log(hc[k]['gpa'] / rc[k]['gpa'])), se2=se(hc[k]['runs']) ** 2 + se(rc[k]['runs']) ** 2)


def part_het():
    mm = [x for x in (lr320(r) for r in PJ1['MM']) if x is not None]; ht = [x for x in (lr320(r) for r in PJ1['HT']) if x is not None]
    v_mm = float(np.var([x['lr'] for x in mm], ddof=1)); n_mm = float(np.mean([x['se2'] for x in mm])); sd_real = float(np.sqrt(max(v_mm - n_mm, 0))); sd_th = sd_real / B_SLOPE
    v_ht = float(np.var([x['lr'] for x in ht], ddof=1)); n_ht = float(np.mean([x['se2'] for x in ht]))
    nz = [z for z in (sim_noise(SP, r) for r in PJ1['HT']) if z is not None]; nzm = [z for z in (sim_noise(SP, r) for r in PJ1['MM']) if z is not None]
    rng = np.random.default_rng(20261046); rows = []
    def simit(zs, th0, sdv):
        items = {}
        for z in zs:
            th = th0 + rng.normal(0, sdv); er = rng.normal(0, z['se_r'])
            items.setdefault(z['H'], []).append(theta([B_SLOPE * (T - th) + rng.normal(0, s) - er for T, s in zip(z['temps'], z['se_h'])], z['temps'])[:2])
        return items
    for th0 in (335.0, 350.0, 380.0):
        ok = []
        for _ in range(100):
            it = {h: sf_summary(v, 'ht') for h, v in simit(nz, th0, sd_th).items()}; ok.append(cond1(it, rng, B=1000)['passed'])
        rows.append(dict(test='cond1', theta=th0, sd_theta=sd_th, power=float(np.mean(ok)))); log('p1 het', rows[-1])
    for dlt in (15.0, 21.0, 30.0):
        ok = []
        for _ in range(100):
            T_s = {h: sf_summary(v, 'ht') for h, v in simit(nz, 320.0 + dlt, sd_th).items()}; C_s = {h: sf_summary(v, 'mm') for h, v in simit(nzm, 320.0, sd_th).items()}
            ok.append(cond2(T_s, C_s, rng, B=1000)['passed'])
        rows.append(dict(test='cond2', delta=dlt, sd_theta=sd_th, power=float(np.mean(ok)))); log('p1 het', rows[-1])
    OUTD['het'] = dict(n_mm=len(mm), n_ht=len(ht), var_lr320_mm=v_mm, noise_var_mm=n_mm, sd_real_lr=sd_real, sd_theta_K=sd_th, var_lr320_ht=v_ht, noise_var_ht=n_ht,
                       sd_theta_assumed_in_registered_sim=15.0, power=rows, note='결과 뒤 진단: 판정을 바꾸지 않음')


def part_win():
    A, H = ('abs', 0.75), MAIN; out = {}
    for side, prs in (('HT', PJ1['HT']), ('MM', PJ1['MM'])):
        th_h = [pair_theta(SP, r, H) for r in prs]; th_a = [pair_theta(SP, r, A) for r in prs]; both = [(x, y) for x, y in zip(th_h, th_a) if x[0] != 'skip' and y[0] != 'skip']
        cen = lambda ts, k: float(np.mean([t[0] == k for t in ts])) if ts else None
        nn = [(y[1] - x[1]) for x, y in both if x[0] == 'num' and y[0] == 'num']
        l_h = [lr320(r, H) for r in prs]; l_a = [lr320(r, A) for r in prs]; ll = [(x['lr'], y['lr']) for x, y in zip(l_h, l_a) if x and y]
        se_h = [x['se2'] for x in l_h if x]; se_a = [x['se2'] for x in l_a if x]
        out[side] = dict(n_both=len(both), left_half=cen([x for x, _ in both], 'lt'), left_abs=cen([y for _, y in both], 'lt'), right_half=cen([x for x, _ in both], 'gt'), right_abs=cen([y for _, y in both], 'gt'),
                         n_num_both=len(nn), median_theta_abs_minus_half=float(np.median(nn)) if nn else None, lr320_corr=float(np.corrcoef(*zip(*ll))[0, 1]) if len(ll) > 2 else None,
                         median_lr320_half=float(np.median([a for a, _ in ll])) if ll else None, median_lr320_abs=float(np.median([b for _, b in ll])) if ll else None,
                         noise_sd_lr320_half=float(np.sqrt(np.median(se_h))) if se_h else None, noise_sd_lr320_abs=float(np.sqrt(np.median(se_a))) if se_a else None)
    keepH = [r for r in PJ1['HT'] if pair_theta(SP, r, H)[0] != 'skip' and pair_theta(SP, r, A)[0] != 'skip']; keepM = [r for r in PJ1['MM'] if pair_theta(SP, r, H)[0] != 'skip' and pair_theta(SP, r, A)[0] != 'skip']
    for nm, cell in (('half_same_pairs', H), ('abs_same_pairs', A)):
        v = verdict_set(keepH, keepM, SP, cell=cell, seed=20261047)[0]; out[nm] = dict(cond1=dict(median=v['cond1']['median'], lo=v['cond1']['ci95'][0], left=v['cond1']['frac_left_cens']),
                                                                               cond2=dict(d=v['cond2'].get('d'), lo=v['cond2'].get('ci95', [None])[0], version=v['cond2'].get('version'), n_sf_both=v['cond2'].get('n_sf_both')))
    OUTD['win'] = dict(out, note='결과 뒤 진단: 같은 쌍 집합에서 두 창을 견줌(쌍 구성 차를 뗌). 50–100 ns 창은 N = 50, 뒤 절반은 N = 100')
    log('p1 win', json.dumps(OUTD['win'], ensure_ascii=False, default=str)[:1500])


def part_unc():
    rng = np.random.default_rng(20261045); out = {}
    for cell_nm, cell in (('main', MAIN), ('q070', ('half', 0.7)), ('q080', ('half', 0.8)), ('abs_50_100ns', ('abs', 0.75))):
        def sfs_c(prs):
            by = {}
            for x in (lr320(r, cell) for r in prs):
                if x is not None: by.setdefault(x['H'], []).append(x['lr'])
            return {h: float(np.median(v)) for h, v in by.items()}
        T_s, C_s = sfs_c(PJ1['HT']), sfs_c(PJ1['MM']); tv = np.array(list(T_s.values())); both = sorted(set(T_s) & set(C_s)); dv = np.array([T_s[h] - C_s[h] for h in both])
        bt = [np.median(tv[rng.integers(0, len(tv), len(tv))]) for _ in range(2000)]; bd = [np.median(dv[rng.integers(0, len(dv), len(dv))]) for _ in range(2000)] if len(dv) else []
        u1 = dict(n_sf=len(tv), median=float(np.median(tv)), ci95=[float(x) for x in np.percentile(bt, [2.5, 97.5])])
        u2 = dict(n_sf=len(dv), median=float(np.median(dv)) if len(dv) else None, ci95=[float(x) for x in np.percentile(bd, [2.5, 97.5])] if len(dv) else None)
        for u in (u1, u2):
            if u['median'] is not None: u['K'] = -u['median'] / B_SLOPE; u['K_ci95'] = sorted([-x / B_SLOPE for x in u['ci95']])
        out[cell_nm] = dict(U1=u1, U2=u2, mm_sf_median=float(np.median(list(C_s.values()))) if C_s else None)
        log('p1 unc', cell_nm, out[cell_nm])
    OUTD['unc'] = dict(out, rule='결과 뒤 분석(정의는 계산 전 Drive JSON round36_p1_diag_2026-10-20): U1 = 호열–중온 상과 요약(쌍 LR(320) 중앙값)의 상과 중앙값, U2 = 짝지은 상과 차의 상과 중앙값; 상과 부트스트랩 2,000; 음수 = 같은 320 K 에서 호열 쪽이 덜 퍼짐; K = −LR/b')


if __name__ == '__main__':
    for part in PARTSD:
        log('부분 시작', part)
        try: globals()['part_' + part](); jdump(OUTD, FO, bk=True); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('P1_DIAG_DONE', ','.join(PARTSD))
