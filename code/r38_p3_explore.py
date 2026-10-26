"""(탐색 — 등록 판정 아님) 1차 P3 가 '판정 불가' 로 끝난 뒤, 검정력이 어디서 새는지 보는 진단. results/r38_p3_explore.json.
r38_p3.py 의 analyse 와 같은 점추정(순위의 가중 편 Pearson, 쌍 무게 1/상과 쌍 수)에 부트스트랩만 바꿔 본다:
  reg   : 등록 판 — 상과 v ~ Exp(1) × 두 도메인 w ~ Exp(1) (재현 확인)
  sf    : 상과 v ~ Exp(1) 만 (H-풀 쌍은 모두 같은 상과 안이라 도메인이 상과에 포함됨 — 내포 구조의 한 단계 군집 부트스트랩)
  dom   : 두 도메인 w 만
  sfres : 상과를 다항으로 다시 뽑기(고전 군집 부트스트랩)
그리고 Kish 실효 n(무게만), 상과 크기별 추정, 같은 무게에서 상과 수를 k 배로 늘렸을 때의 어림 반폭.
사용: CUDA_VISIBLE_DEVICES= python r38_p3_explore.py <ROOT>"""
import sys, json, time
ROOT = sys.argv[1]
P3 = {'__name__': 'p3x'}; _s = open(ROOT + '/code/r38_p3.py').read(); _s = _s[:_s.index("ALL = ['pairs'")]
sys.argv = ['r38_p3', ROOT, 'none', '1']; exec(_s, P3)
np = P3['np']; ranks, pcorr = P3['ranks'], P3['pcorr']
rows = [json.loads(l) for l in open(P3['OUT'] + '/pairs_out.jsonl')]; rows = [r for r in rows if r['main'].get('r_dyn') is not None]
Y = ranks(np.array([r['main']['r_dyn'] for r in rows])); X = ranks(np.array([r['main']['ident'] for r in rows])); Z = [ranks(np.array([r['main']['control'] for r in rows]))]
H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); cnt = np.bincount(si); base = 1.0 / cnt[si]
doms = sorted({r['ab'][0] for r in rows} | {r['ab'][1] for r in rows}); di = {d: i for i, d in enumerate(doms)}
ia = np.array([di[r['ab'][0]] for r in rows]); ib = np.array([di[r['ab'][1]] for r in rows])
est = pcorr(Y, X, Z, base); out = dict(written=time.strftime('%F %T'), n=len(rows), n_sf=len(sfl), n_dom=len(doms), est=est,
                                       kish_n_eff_weights_only=float(base.sum() ** 2 / (base ** 2).sum()), note='탐색. 등록 판정(r38_p3.json)을 바꾸지 않는다')
rng = np.random.default_rng(20261050); B = 2000; bs = {k: [] for k in ('reg', 'sf', 'dom', 'sfres')}
for b in range(B):
    v = rng.exponential(1, len(sfl)); w = rng.exponential(1, len(doms))
    bs['reg'].append(pcorr(Y, X, Z, base * v[si] * w[ia] * w[ib])); bs['sf'].append(pcorr(Y, X, Z, base * v[si])); bs['dom'].append(pcorr(Y, X, Z, base * w[ia] * w[ib]))
    c = np.bincount(rng.integers(0, len(sfl), len(sfl)), minlength=len(sfl)); bs['sfres'].append(pcorr(Y, X, Z, base * c[si]) if c[si].sum() > 0 else np.nan)
for k, v in bs.items():
    v = np.array(v, float); v = v[np.isfinite(v)]; lo, hi = np.percentile(v, [2.5, 97.5]); se = float(v.std())
    out[k] = dict(lo=float(lo), hi=float(hi), se=se, n_eff_fisher=float(1.06 / se ** 2 + 3), lo_one_sided_95=float(np.percentile(v, 5)))
# 상과 크기별 (쌍 수 1–2 / 3–10 / > 10)
for name, m in (('sf_pairs_1_2', cnt[si] <= 2), ('sf_pairs_3_10', (cnt[si] >= 3) & (cnt[si] <= 10)), ('sf_pairs_gt10', cnt[si] > 10)):
    out[name] = dict(n=int(m.sum()), n_sf=int(len(set(si[m]))), est=float(pcorr(Y[m], X[m], [Z[0][m]], base[m])) if m.sum() > 20 else None)
# 상과 수를 k 배로: 상과 군집 SE 는 1/√(상과 수) 에 비례한다고 보고 반폭 어림
se_sf = out['sf']['se']; se_reg = out['reg']['se']
out['scale_n_sf'] = {f'x{k}': dict(n_sf=len(sfl) * k, half_width_sf=float(1.96 * se_sf / np.sqrt(k)), half_width_reg=float(1.96 * se_reg / np.sqrt(k))) for k in (1, 2, 3, 4)}
out['sf_needed_for_lo_gt0_at_est'] = dict(sf=int(np.ceil(len(sfl) * (1.96 * se_sf / est) ** 2)), reg=int(np.ceil(len(sfl) * (1.96 * se_reg / est) ** 2)))
json.dump(out, open(ROOT + '/results/r38_p3_explore.json', 'w'), ensure_ascii=False, indent=1, default=float)
print(json.dumps(out, ensure_ascii=False, indent=1, default=float)); print('P3X_DONE')
