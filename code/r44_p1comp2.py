"""(round44, 등록 = Drive results/prereg_2026-10-28_addendum5.json "r44_p1comp" + 고침 _addendum5_fix1.json 08:59:08 UTC) P1′ 재현의 쌍 구성 점검.
r44_w2p1.py 의 부분 p1comp 는 pairs.json 의 tm·ident 칸이 모두 비어 있어 값을 내기 전에 멈췄다. 고침대로: 구조 유사도 = 평균 점 사이 정규화 편집 거리 dhat
(= √(max(0, 1 − 2S/(‖m_a‖² + ‖m_b‖²))), S = align.jsonl 의 평균 점 정렬 점수, m = mean/<d>.npz), ident 는 뺌, 보정 회귀 변수 = [log n_shared, 공유 비율, dhat].
쓰는 쌍·통계·seed 는 등록 그대로: 쌍 = r38_p1n.lr320_of 가 값을 내는 쌍(320 K, 뒤 절반, Q 0.7); (가) 상과 짝지은 차(HT − MM)의 상과 중앙, 재추출 2,000 seed 20261101;
(나) HT 안·MM 안 Spearman(로그 비, 구성); (다) MM 에서 맞춘 선형식을 뺀 U2 (상과 중앙 차의 상과 중앙, 재추출 2,000 seed 20261049, 한쪽 97.5%).
사용: python r44_p1comp2.py <ROOT>   출력: results/r44_p1comp.json"""
import os, sys, json, time
import numpy as np
ROOT = sys.argv[1]; RES44 = ROOT + '/results'


def rk(v): return np.argsort(np.argsort(np.asarray(v, float))).astype(float)
def wcor(u, v, w):
    w = w / w.sum(); u = u - (w * u).sum(); v = v - (w * v).sum(); return float((w * u * v).sum() / np.sqrt((w * u * u).sum() * (w * v * v).sum()))


_s = open(ROOT + '/code/r38_p1n.py').read(); _s = _s[:_s.index("ALLN = [")]
sys.argv = ['r38_p1n', ROOT, 'none', '1']; g = {'__name__': 'r38_p1n_lib'}; exec(_s, g)   # lr320_of, load_spread, sf_med, OUT1 (= derived/r38/p1n)
lr320_of, sf_med, OUT1 = g['lr320_of'], g['sf_med'], g['OUT1']
pj = json.load(open(OUT1 + '/pairs.json')); SP = g['load_spread'](pj['domains']); cell = ('half', 0.7)
AL = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(OUT1 + '/align.jsonl')}; NM = {}
def norm2(d):
    if d not in NM: NM[d] = float((np.load(f'{OUT1}/mean/{d}.npz')['m'].astype(np.float64) ** 2).sum())
    return NM[d]


def rows_of(prs):
    out = []
    for r in prs:
        v = lr320_of(SP, r, cell)
        if v is None or not np.isfinite(v): continue
        ns = SP[r['hot']]['pairs'][f"{r['a']}__{r['b']}"]['n_shared']; S = float(AL[(r['a'], r['b'])]['S'])
        dh = float(np.sqrt(max(0.0, 1 - 2 * S / (norm2(r['a']) + norm2(r['b'])))))
        out.append(dict(r=r, lr=v, logn=float(np.log(ns)), frac=float(ns / min(r['La'] - 1, r['Lb'] - 1)), dhat=dh, logL=float(np.log(max(r['La'], r['Lb'])))))
    return out


HT, MM = rows_of(pj['HT']), rows_of(pj['MM']); res = dict(written=time.strftime('%F %T UTC', time.gmtime()), n_HT=len(HT), n_MM=len(MM))
comp = {}; rng = np.random.default_rng(20261101)
for k in ('logn', 'frac', 'dhat', 'logL'):
    th = sf_med([dict(H=x['r']['H'], v=x[k]) for x in HT], lambda q: q['v']); tc = sf_med([dict(H=x['r']['H'], v=x[k]) for x in MM], lambda q: q['v'])
    both = sorted(set(th) & set(tc)); dv = np.array([th[h] - tc[h] for h in both]); bd = [np.median(dv[rng.integers(0, len(dv), len(dv))]) for _ in range(2000)]
    comp[k] = dict(n_sf=len(both), HT_median=float(np.median([x[k] for x in HT])), MM_median=float(np.median([x[k] for x in MM])),
                   sf_paired_diff=float(np.median(dv)), lo=float(np.percentile(bd, 2.5)), hi=float(np.percentile(bd, 97.5)))
res['composition'] = comp
sp = lambda xs, k: wcor(rk([x['lr'] for x in xs]), rk([x[k] for x in xs]), np.ones(len(xs)))
res['pair_level_spearman'] = {f'{grp}|{k}': sp(xs, k) for grp, xs in (('HT', HT), ('MM', MM)) for k in ('logn', 'frac', 'dhat')}
A = lambda xs: np.c_[np.ones(len(xs)), [x['logn'] for x in xs], [x['frac'] for x in xs], [x['dhat'] for x in xs]]
beta = np.linalg.lstsq(A(MM), np.array([x['lr'] for x in MM]), rcond=None)[0]


def u2(adj):
    f = (lambda x: x['lr'] - float(A([x])[0] @ beta)) if adj else (lambda x: x['lr'])
    T = sf_med([dict(H=x['r']['H'], x=x) for x in HT], lambda q: f(q['x'])); C = sf_med([dict(H=x['r']['H'], x=x) for x in MM], lambda q: f(q['x']))
    both = sorted(set(T) & set(C)); dv = np.array([T[h] - C[h] for h in both]); rg = np.random.default_rng(20261049)
    bd = [np.median(dv[rg.integers(0, len(dv), len(dv))]) for _ in range(2000)]
    return dict(n_sf=len(both), U2=float(np.median(dv)), lo025=float(np.percentile(bd, 2.5)), hi975=float(np.percentile(bd, 97.5)))


res['beta_MM'] = dict(intercept=float(beta[0]), logn=float(beta[1]), frac=float(beta[2]), dhat=float(beta[3]))
res['U2_main_recomputed'] = u2(False); res['U2_adjusted'] = u2(True)
a, m = res['U2_adjusted'], res['U2_main_recomputed']
res['predictions'] = dict(pc1_adjusted_upper_negative=bool(a['hi975'] < 0), pc2_change_lt_005=bool(abs(a['U2'] - m['U2']) < 0.05))
res['rule'] = ('등록 r44_p1comp (prereg_2026-10-28_addendum5.json) + 고침 _addendum5_fix1.json: PC1 보정 판 U2 의 97.5% 상한 < 0; '
               'PC2 |보정 판 U2 − 주 판 U2| < 0.05. 구조 유사도는 tm 대신 dhat(평균 점 사이 정규화 편집 거리).')
json.dump(res, open(RES44 + '/r44_p1comp.json', 'w'), ensure_ascii=False, indent=1); print(json.dumps(res, ensure_ascii=False, indent=1), flush=True)
