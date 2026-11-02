
"""(prereg round29 E0_and_HT) E0 변 단위 퍼짐 기준 통과율, 에너지 거리 Δ²>0 비율, H-T 세균만 + 문턱. results/e0_ht_extra.json."""
import sys, os, json, glob, time, itertools
ROOT = sys.argv[1]
import numpy as np, pandas as pd
rng = np.random.default_rng(29); lg = lambda x: np.log(np.maximum(x, 1e-12)); OUT = {}
st = json.load(open(ROOT + '/results/mdcath_5T_setup.json')); E0SET, C1SET = set(st['e0_domains']), set(st['c1_domains'])
# (1) 퍼짐 기준 (mdcath_spread, 전체 프레임)
splits = [c for c in itertools.combinations(range(5), 2)]
def stat(a, b): return float(np.median(lg(a) - lg(b)))
rows = []
for p in sorted(glob.glob(ROOT + '/derived/mdcath_spread/*.npz')):
    z = np.load(p)
    if not all(k in z.files for k in ('320_spread_rep', '348_spread_rep')) or z['320_spread_rep'].shape[0] != 5 or z['348_spread_rep'].shape[0] != 5: continue
    s3, s4 = z['320_spread_rep'], z['348_spread_rep']; obs = stat(s4.mean(0), s3.mean(0)); null = []
    for S_ in (s3, s4):
        for g in splits:
            h = [k for k in range(5) if k not in g]; null.append(abs(stat(S_[list(g)].mean(0), S_[h].mean(0))))
    d = os.path.basename(p)[:-4]; rows.append(dict(dom=d, obs=obs, null_max=max(null), null_q95=float(np.quantile(null, 0.95)), e0=d in E0SET, c1=d in C1SET))
Dp = pd.DataFrame(rows); Dp['pass'] = Dp.obs > Dp.null_max
def frac_ci(v):
    v = np.asarray(v, float); b = [v[rng.integers(0, len(v), len(v))].mean() for _ in range(2000)]; return dict(n=int(len(v)), frac=round(float(v.mean()), 3), CI95=[round(float(x), 3) for x in np.percentile(b, [2.5, 97.5])])
OUT['E0_spread'] = dict(all=frac_ci(Dp['pass']), e0_set=frac_ci(Dp[Dp.e0]['pass']), c1_set=frac_ci(Dp[Dp.c1]['pass']), obs_median=round(float(Dp.obs.median()), 4), harmonic=round(float(np.log(348 / 320)), 4),
                        note='통계 = 변 중앙값 log(평균 퍼짐 348 / 320) (반복 5씩), 귀무 = 같은 온도 반복 2 대 3 나눔의 같은 통계 절댓값 20개, 통과 = 통계 > 귀무 최댓값. 전체 프레임(Q 거르기 없음).')
# (2) 에너지 거리 Δ² (ladder_temp_single_v2, 접힌 프레임)
L = [json.load(open(q)) for q in glob.glob(ROOT + '/derived/ladder_temp_single_v2/*.json')]
d2 = [(j['cross1']['348'] - (j['floor1'] + j['floorT']['348']) / 2) / np.sqrt(j['ne']) for j in L if '348' in j.get('cross1', {}) and '348' in j.get('floorT', {})]
mult = [j['cross1']['348'] / j['floor1'] for j in L if '348' in j.get('cross1', {})]
OUT['E0_energy'] = dict(frac_d2_pos=frac_ci(np.array(d2) > 0), d2n_median=round(float(np.median(d2)), 4), single_multiple_median=round(float(np.median(mult)), 3), frac_cross_gt_floor320=frac_ci(np.array(mult) > 1), n=len(d2), note='접힌 프레임, 단일 반복 짝 평균(쌍마다 값은 저장되지 않아 95 백분위 판정은 못 함 — 서술).')
# (3) H-T 세균만 + 문턱
src = open(ROOT + '/code/c3_saturation_ht_variants.py').read(); i0 = src.index('# (2) H-T 변형'); i1 = src.index('P = D[D.prok == True]')
SPD = ROOT + '/derived/mdcath_spread'
cath = {}
for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
    if l.startswith('#') or not l.strip(): continue
    f = l.split(); cath[f[0]] = '.'.join(f[1:5])
exec(src[i0:i1])
P = D[D.prok == True]; B = P[(P['T'] >= 60) | (P['T'] < 40)].assign(hot=lambda d: (d['T'] >= 60).astype(float))
OUT['HT_bacteria_threshold'] = dict(bacteria_threshold_no_lineage=run(B[B.sk == 'Bacteria'], ['hot'], lineage=False, scale=1.0), counts=dict(bacteria_hot=int(((B.sk == 'Bacteria') & (B.hot == 1)).sum()), bacteria_cold=int(((B.sk == 'Bacteria') & (B.hot == 0)).sum())),
                                    note='호열(≥60) 대 중온(<40) 차(ln 퍼짐 348), 세균만, 상과 고정효과, 상과 부트스트랩 500. 쌍 분석 세균–세균 −0.233 과 견줌.')
OUT['written'] = time.strftime('%F %T'); json.dump(OUT, open(ROOT + '/results/e0_ht_extra.json', 'w'), indent=1, ensure_ascii=False); print(json.dumps(OUT, indent=1, ensure_ascii=False))
