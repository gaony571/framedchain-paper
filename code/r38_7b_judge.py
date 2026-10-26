"""(prereg round36 7b) 판정용 나눔 출력의 요약과 판정. 규칙은 2026-10-19 14:56:51 UTC 부분 동결 원문
(Drive results/prereg_2026-10-17.json 의 round36_partial_freeze_7b7c → sections.sec_7b)을 따른다.
입력: derived/r38/b7/<표적>__s<k>.json (r38_7b_split.py; 82 표적 × 20 나눔 × 대조 8)
출력: results/r38_7b_summary.json, results/r38_7b_summary.txt (사본 results/_backup/<날짜>/)
가까움의 방향:
  거리꼴(작을수록 가깝다): E_D2(에너지 거리 Δ²), rmwd, emd_mean(RMWD 평균 항), emd_var(RMWD Bures 항), refemd(MD PCA W₂), jointemd(결합 PCA W₂),
                         pairwise_absdiff = |모델 쌍 RMSD − MD 쌍 RMSD| (쌍 RMSD 를 표적별로 읽는 꼴; 이 스크립트에서 정의)
  유사꼴(클수록 가깝다): rmsf_pearson(표적별 RMSF r), cosine_sim(PC 유사도), weak_contacts_iou, transient_contacts_iou (접촉 Jaccard)
표적마다 나눔마다 "대조가 진짜보다 가깝다"(동점은 0.5) → 20 나눔 평균 p_t (prereg: 표적마다 판정 비율을 평균한다).
지표·대조의 "가깝다고 하는 표적 비율" = p_t 의 표적 평균, 표적 부트스트랩 2000 (seed 20261025) 95% 구간. "멀다고 한다" = 1 − 그것.
볼 수 없음: prereg 가 식으로 정한 동점 쌍(표적별 RMSF r·PC 유사도 × 줄임 0.8·0.6·부풂 1.2)은 이름으로 "볼 수 없음 (식으로 정해짐)" 이고,
  실측 동점 비율(|대조 − 진짜| ≤ 1e-6 · max(1, |진짜|))과 상대 차 중앙값은 확인으로만 낸다(원 코드가 줄인 구조를 결정에 다시 중첩해 정확히 같지 않을 수 있음).
  그 밖의 쌍은 실측 동점이 99% 이상일 때만 "볼 수 없음 (실측 동점)".
판정: 오판 = 가깝다고 하는 표적 비율의 95% 하한 > 0.5. 시험 예측 1 (E_D2, 모든 대조, 멀다 ≥ 0.95), 2 (rmwd, 줄임·잘라내기, ≥ 0.70),
      3 (접촉 Jaccard 둘, 섞음, ≥ 0.90), 4 (PCA W₂ 둘, 줄임에서 가깝다고 하는 비율, 미리 정하지 않음).
식으로 정해짐: 진짜 후보 Δ²/√변 — 표적마다 20 나눔 평균 → 82 표적 평균, 표준오차 = 표적 사이 sd/√82, |평균| < 2 표준오차면 선다.
표적 사이 지표(나눔마다 82 표적으로): 쌍 RMSD r (Pearson), PC sim > 0.5 비율. 전체 RMSF r 는 원자별 RMSF 를 저장하지 않아 낼 수 없다(줄임·부풂의 동점은 정의로 섬).
사용: python r38_7b_judge.py <ROOT>"""
import os, sys, json, glob, time, shutil
import numpy as np
ROOT = sys.argv[1]; IN = ROOT + '/derived/r38/b7'; RES = ROOT + '/results'
CT = ['shrink0.8', 'shrink0.6', 'inflate1.2', 'trim20', 'trim40', 'scramble', 'xtalnoise']
DIST = ['E_D2', 'rmwd', 'emd_mean', 'emd_var', 'refemd', 'jointemd', 'pairwise_absdiff']
SIM = ['rmsf_pearson', 'cosine_sim', 'weak_contacts_iou', 'transient_contacts_iou']
NB = 2000; RNG = np.random.default_rng(20261025)
FORMULA = {(m, c) for m in ('rmsf_pearson', 'cosine_sim') for c in ('shrink0.8', 'shrink0.6', 'inflate1.2')}

J = {}
for f in sorted(glob.glob(IN + '/*__s*.json')):
    j = json.load(open(f)); J.setdefault(j['name'], {})[j['split']] = j['controls']
names = sorted(J); nsplit = {n: len(J[n]) for n in names}
def val(c, m):
    if m == 'pairwise_absdiff': return abs(c['af_pairwise'] - c['md_pairwise'])
    v = c.get(m); return float('nan') if v is None else float(v)

def boot(pt):
    pt = np.asarray([x for x in pt if np.isfinite(x)])
    if len(pt) == 0: return dict(mean=None, lo=None, hi=None, n_targets=0)
    bs = np.array([RNG.choice(pt, len(pt)).mean() for _ in range(NB)])
    return dict(mean=float(pt.mean()), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), n_targets=int(len(pt)))

table = {}
for m in DIST + SIM:
    sign = -1 if m in DIST else 1                       # 가깝다 = sign·(대조 − 진짜) > 0
    for c in CT:
        pt, ties, cells, rel = [], 0, 0, []
        for n in names:
            v = []
            for s, ctl in J[n].items():
                a, t = val(ctl[c], m), val(ctl['true'], m)
                if not (np.isfinite(a) and np.isfinite(t)): continue
                cells += 1; rel.append(abs(a - t) / max(abs(t), 1e-12))
                if abs(a - t) <= 1e-6 * max(1.0, abs(t)): ties += 1; v.append(0.5)
                else: v.append(1.0 if sign * (a - t) > 0 else 0.0)
            pt.append(np.mean(v) if v else float('nan'))
        b = boot(pt); tf = ties / cells if cells else None
        verdict = '볼 수 없음 (식으로 정해짐)' if (m, c) in FORMULA else ('볼 수 없음 (실측 동점)' if tf is not None and tf >= 0.99 else ('오판' if b['lo'] is not None and b['lo'] > 0.5 else '아님'))
        table[f'{m}|{c}'] = dict(metric=m, control=c, closer=b, farther=None if b['mean'] is None else dict(mean=1 - b['mean'], lo=1 - b['hi'], hi=1 - b['lo']),
                                tie_frac=tf, median_rel_diff=float(np.median(rel)) if rel else None, n_cells=cells, verdict=verdict)

def far(m, c): return table[f'{m}|{c}']['farther']
tests = {}
t1 = {c: far('E_D2', c)['mean'] for c in CT}
tests['1_energy_distance_all_controls_farther_ge_0.95'] = dict(by_control=t1, passed={c: bool(v >= 0.95) for c, v in t1.items()}, all_passed=bool(all(v >= 0.95 for v in t1.values())))
t2 = {c: far('rmwd', c)['mean'] for c in ('shrink0.8', 'shrink0.6', 'trim20', 'trim40')}
tests['2_rmwd_shrink_trim_farther_ge_0.70'] = dict(by_control=t2, passed={c: bool(v >= 0.70) for c, v in t2.items()}, all_passed=bool(all(v >= 0.70 for v in t2.values())))
t3 = {m: far(m, 'scramble')['mean'] for m in ('weak_contacts_iou', 'transient_contacts_iou')}
tests['3_contact_jaccard_scramble_farther_ge_0.90'] = dict(by_metric=t3, passed={m: bool(v >= 0.90) for m, v in t3.items()}, all_passed=bool(all(v >= 0.90 for v in t3.values())))
tests['4_pca_w2_shrink_closer_fraction (미리 정하지 않음)'] = {f'{m}|{c}': table[f'{m}|{c}']['closer'] for m in ('refemd', 'jointemd') for c in ('shrink0.8', 'shrink0.6')}

# 시험 1 의 여유: 표적마다 (대조 − 진짜) Δ²/√변 의 20 나눔 평균 ÷ 표준오차
margin = {}
for c in CT:
    z = []
    for n in names:
        d = np.array([J[n][s][c]['E_d2n'] - J[n][s]['true']['E_d2n'] for s in J[n]])
        if len(d) > 1 and d.std(ddof=1) > 0: z.append(d.mean() / (d.std(ddof=1) / np.sqrt(len(d))))
    z = np.array(z); margin[c] = dict(median=float(np.median(z)), frac_gt_2=float((z > 2).mean()), frac_lt_0=float((z < 0).mean()), n=int(len(z)))

# 식으로 정해짐: 진짜 후보 Δ²/√변 ≈ 0
tm = np.array([np.mean([J[n][s]['true']['E_d2n'] for s in J[n]]) for n in names])
se = float(tm.std(ddof=1) / np.sqrt(len(tm)))
sanity = dict(mean=float(tm.mean()), se=se, z=float(tm.mean() / se), n_targets=int(len(tm)), stands=bool(abs(tm.mean()) < 2 * se),
              pilot=[0.0036, -0.0016, 0.0006], rule='|평균| < 2 표준오차면 선다')

# 표적 사이 지표(나눔마다 82 표적)
splits = sorted({s for n in names for s in J[n]})
def cross(c):
    r, pc = [], []
    for s in splits:
        ns_ = [n for n in names if s in J[n]]
        a = np.array([J[n][s][c]['af_pairwise'] for n in ns_]); md = np.array([J[n][s][c]['md_pairwise'] for n in ns_])
        r.append(float(np.corrcoef(a, md)[0, 1])); pc.append(float(np.mean([J[n][s][c]['cosine_sim'] > 0.5 for n in ns_])))
    return np.array(r), np.array(pc)
rt, pct = cross('true'); acr = {}
for c in CT:
    r, pc = cross(c)
    acr[c] = dict(pairwise_r_mean=float(r.mean()), pairwise_r_true=float(rt.mean()), pairwise_r_tie=bool(np.allclose(r, rt, rtol=1e-6, atol=1e-9)),
                  pc_sim_gt05_mean=float(pc.mean()), pc_sim_gt05_true=float(pct.mean()), pc_sim_tie=bool(np.allclose(pc, pct, atol=1e-12)))

bures = int(sum(J[n][s][c]['bures_fallback'] for n in names for s in J[n] for c in ['true'] + CT))
mis = sorted(k for k, v in table.items() if v['verdict'] == '오판'); cannot = sorted(k for k, v in table.items() if v['verdict'].startswith('볼 수 없음'))
out = dict(written=time.strftime('%F %T'), rule_source='prereg_2026-10-17.json → round36_partial_freeze_7b7c.sections.sec_7b (2026-10-19 14:56:51 UTC)',
           n_targets=len(names), splits_per_target=sorted(set(nsplit.values())), n_files=sum(nsplit.values()), bures_fallback_total=bures,
           sanity_true_candidate=sanity, tests=tests, energy_margin=margin, table=table, across_targets=acr,
           misjudged=mis, cannot_see=cannot, box_title='상자 유지' if (mis or cannot) else '"고차원 경험 W₂ 계열의 주의" 로 바꿈 (오판·볼 수 없음 둘 다 없음)')
P = RES + '/r38_7b_summary.json'; json.dump(out, open(P, 'w'), ensure_ascii=False, indent=1)
L = [f"7b 판정 요약 ({out['written']}) — 표적 {len(names)}, 나눔 {out['splits_per_target']}, 파일 {out['n_files']}, Bures 대체 {bures}",
     f"식으로 정해짐: 진짜 후보 Δ²/√변 평균 {sanity['mean']:.5f} ± {se:.5f} (z {sanity['z']:.2f}) → {'선다' if sanity['stands'] else '서지 않음'}"]
for k, v in tests.items():
    if k.startswith('4'): L.append(f"시험 {k}: " + ', '.join(f"{kk} {vv['mean']:.2f} [{vv['lo']:.2f}, {vv['hi']:.2f}]" for kk, vv in v.items())); continue
    L.append(f"시험 {k}: " + ', '.join(f"{kk} {vv:.3f}" for kk, vv in (v.get('by_control') or v.get('by_metric')).items()) + f" → {'선다' if v['all_passed'] else '서지 않음'}")
L.append('에너지 거리 여유(표적별 z 중앙값, z>2 비율): ' + ', '.join(f"{c} {margin[c]['median']:.1f}/{margin[c]['frac_gt_2']:.2f}" for c in CT))
L.append('지표 × 대조: 가깝다고 하는 표적 비율 [95%] (판정)')
for m in DIST + SIM:
    L.append(f"  {m:24s} " + ' | '.join(f"{c} {table[f'{m}|{c}']['closer']['mean']:.2f}[{table[f'{m}|{c}']['closer']['lo']:.2f},{table[f'{m}|{c}']['closer']['hi']:.2f}]{'*오판' if table[f'{m}|{c}']['verdict']=='오판' else ('=동점' if table[f'{m}|{c}']['verdict'].startswith('볼') else '')}" if table[f'{m}|{c}']['closer']['mean'] is not None else f"{c} -" for c in CT))
L.append('표적 사이: ' + ', '.join(f"{c} 쌍RMSD r {acr[c]['pairwise_r_mean']:.3f}{'(동점)' if acr[c]['pairwise_r_tie'] else ''} PC>0.5 {acr[c]['pc_sim_gt05_mean']:.2f}{'(동점)' if acr[c]['pc_sim_tie'] else ''}" for c in CT) + f" | 진짜 쌍RMSD r {acr[CT[0]]['pairwise_r_true']:.3f} PC>0.5 {acr[CT[0]]['pc_sim_gt05_true']:.2f}")
L.append(f"오판: {mis or '없음'} | 볼 수 없음: {cannot or '없음'} | 상자: {out['box_title']}")
open(RES + '/r38_7b_summary.txt', 'w').write('\n'.join(L) + '\n')
bk = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(bk, exist_ok=True)
for f in (P, RES + '/r38_7b_summary.txt'): shutil.copy(f, bk + '/' + time.strftime('%H%M%S_') + os.path.basename(f))
print('\n'.join(L))
