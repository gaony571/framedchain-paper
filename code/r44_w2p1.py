"""(round44, 등록 = Drive results/prereg_2026-10-28_addendum5.json "r44_w2gap"·"r44_p1comp", 2026-10-28 08:53:17 UTC) 외부 검토 열한째의 두 점검.
w2    : 연속 축(B6, W2 가 있는 122 쌍)에서 Spearman(W₂, Δ²_D), 편 Spearman(W₂, Δ²_D | w̄) (보조 | w̄, log n̄), ρ(W₂, w̄)·ρ(Δ²_D, w̄)·ρ(W₂, C).
        구간 = 도메인 Exp(1) 무게의 곱 2,000 (seed 20261100; r44_b6size 와 같은 방법).
p1comp: P1′ 재현(320 K, 뒤 절반, Q 0.7)의 호열–중온(HT)·중온–중온(MM) 쌍 구성(log n_shared, 공유 비율, tm, ident, log 긴 쪽 길이)의 상과 짝지은 차,
        쌍 수준 Spearman(로그 비, 구성), MM 에서 맞춘 선형 보정(log n_shared, 공유 비율, tm)을 뺀 U2 (r38_p1n.ustats 와 같은 상과 중앙 차, 한쪽 97.5%).
사용: python r44_w2p1.py <ROOT> w2,p1comp   출력: results/r44_w2gap.json, results/r44_p1comp.json"""
import os, sys, json, time, glob
import numpy as np
ROOT = sys.argv[1]; PARTS44 = sys.argv[2].split(','); RES44 = ROOT + '/results'


def rk(v): return np.argsort(np.argsort(np.asarray(v, float))).astype(float)
def wcor(u, v, w):
    w = w / w.sum(); u = u - (w * u).sum(); v = v - (w * v).sum(); return float((w * u * v).sum() / np.sqrt((w * u * u).sum() * (w * v * v).sum()))
def wres(u, Z, w):
    A = np.c_[np.ones(len(u)), np.asarray(Z).reshape(len(u), -1)]; W = np.sqrt(w); b = np.linalg.lstsq(A * W[:, None], u * W, rcond=None)[0]; return u - A @ b
def dump(o, f):
    json.dump(o, open(f, 'w'), ensure_ascii=False, indent=1); print(json.dumps(o, ensure_ascii=False, indent=1), flush=True)


def part_w2():
    rows = []
    for f in sorted(glob.glob(f'{ROOT}/derived/r38/b6/*.json')):
        z = json.load(open(f)); Dz = z.get('D') or {}
        if z.get('W2') is not None and Dz.get('D2') is not None and Dz.get('wbar') is not None:
            rows.append((z['a'], z['b'], float(z['W2']), float(Dz['D2']), float(Dz['wbar']), float(Dz['C']), float(z['nbar'])))
    Y, X, Zw, Cc, Zn = rk([r[2] for r in rows]), rk([r[3] for r in rows]), rk([r[4] for r in rows]), rk([r[5] for r in rows]), rk(np.log([r[6] for r in rows]))
    doms = sorted({x for r in rows for x in r[:2]}); ia = np.array([doms.index(r[0]) for r in rows]); ib = np.array([doms.index(r[1]) for r in rows])
    stat = dict(rho_W_D=lambda w: wcor(Y, X, w), rho_W_wbar=lambda w: wcor(Y, Zw, w), rho_D_wbar=lambda w: wcor(X, Zw, w), rho_W_C=lambda w: wcor(Y, Cc, w),
                partial_wbar=lambda w: wcor(wres(Y, Zw, w), wres(X, Zw, w), w), partial_wbar_logn=lambda w: wcor(wres(Y, np.c_[Zw, Zn], w), wres(X, np.c_[Zw, Zn], w), w))
    w1 = np.ones(len(rows)); est = {k: f(w1) for k, f in stat.items()}; rng = np.random.default_rng(20261100); bs = {k: [] for k in stat}
    for _ in range(2000):
        v = rng.exponential(1, len(doms)); w = v[ia] * v[ib]
        for k, f in stat.items(): bs[k].append(f(w))
    out = dict(written=time.strftime('%F %T UTC', time.gmtime()), n_pairs=len(rows), n_dom=len(doms),
               **{k: dict(est=est[k], lo=float(np.percentile(bs[k], 2.5)), hi=float(np.percentile(bs[k], 97.5))) for k in stat})
    out['predictions'] = dict(pw1_partial_ge_08=bool(est['partial_wbar'] >= 0.8), pw2_W_more_with_wbar=bool(est['rho_W_wbar'] > est['rho_D_wbar']))
    out['rule'] = '등록 r44_w2gap (prereg_2026-10-28_addendum5.json): PW1 편 Spearman(W₂, Δ²_D | w̄) 점추정 ≥ 0.8; PW2 ρ(W₂, w̄) > ρ(Δ²_D, w̄).'
    dump(out, RES44 + '/r44_w2gap.json')


def part_p1comp():
    _s = open(ROOT + '/code/r38_p1n.py').read(); _s = _s[:_s.index("ALLN = [")]
    sys.argv = ['r38_p1n', ROOT, 'none', '1']; g = {'__name__': 'r38_p1n_lib'}; exec(_s, g)   # lr320_of, load_spread, sf_med, OUT1 (= derived/r38/p1n)
    lr320_of, sf_med = g['lr320_of'], g['sf_med']; pj = json.load(open(g['OUT1'] + '/pairs.json')); SP = g['load_spread'](pj['domains']); cell = ('half', 0.7)
    def rows_of(prs):
        out = []
        for r in prs:
            v = lr320_of(SP, r, cell)
            if v is None or not np.isfinite(v): continue
            ns = SP[r['hot']]['pairs'][f"{r['a']}__{r['b']}"]['n_shared']
            out.append(dict(r=r, lr=v, logn=float(np.log(ns)), frac=float(ns / min(r['La'] - 1, r['Lb'] - 1)), tm=float(r['tm']), ident=float(r['ident']), logL=float(np.log(max(r['La'], r['Lb'])))))
        return out
    HT, MM = rows_of(pj['HT']), rows_of(pj['MM']); res = dict(written=time.strftime('%F %T UTC', time.gmtime()), n_HT=len(HT), n_MM=len(MM))
    comp = {}; rng = np.random.default_rng(20261101)
    for k in ('logn', 'frac', 'tm', 'ident', 'logL'):
        th = sf_med([dict(H=x['r']['H'], v=x[k]) for x in HT], lambda q: q['v']); tc = sf_med([dict(H=x['r']['H'], v=x[k]) for x in MM], lambda q: q['v'])
        both = sorted(set(th) & set(tc)); dv = np.array([th[h] - tc[h] for h in both]); bd = [np.median(dv[rng.integers(0, len(dv), len(dv))]) for _ in range(2000)]
        comp[k] = dict(n_sf=len(both), HT_median=float(np.median([x[k] for x in HT])), MM_median=float(np.median([x[k] for x in MM])),
                       sf_paired_diff=float(np.median(dv)), lo=float(np.percentile(bd, 2.5)), hi=float(np.percentile(bd, 97.5)))
    res['composition'] = comp
    sp = lambda xs, k: wcor(rk([x['lr'] for x in xs]), rk([x[k] for x in xs]), np.ones(len(xs)))
    res['pair_level_spearman'] = {f'{grp}|{k}': sp(xs, k) for grp, xs in (('HT', HT), ('MM', MM)) for k in ('logn', 'frac', 'tm')}
    A = lambda xs: np.c_[np.ones(len(xs)), [x['logn'] for x in xs], [x['frac'] for x in xs], [x['tm'] for x in xs]]
    beta = np.linalg.lstsq(A(MM), np.array([x['lr'] for x in MM]), rcond=None)[0]
    def u2(adj):
        f = (lambda x: x['lr'] - float(A([x])[0] @ beta)) if adj else (lambda x: x['lr'])
        T = sf_med([dict(H=x['r']['H'], x=x) for x in HT], lambda q: f(q['x'])); C = sf_med([dict(H=x['r']['H'], x=x) for x in MM], lambda q: f(q['x']))
        both = sorted(set(T) & set(C)); dv = np.array([T[h] - C[h] for h in both]); rg = np.random.default_rng(20261049)
        bd = [np.median(dv[rg.integers(0, len(dv), len(dv))]) for _ in range(2000)]
        return dict(n_sf=len(both), U2=float(np.median(dv)), lo025=float(np.percentile(bd, 2.5)), hi975=float(np.percentile(bd, 97.5)))
    res['beta_MM'] = dict(intercept=float(beta[0]), logn=float(beta[1]), frac=float(beta[2]), tm=float(beta[3]))
    res['U2_main_recomputed'] = u2(False); res['U2_adjusted'] = u2(True)
    a, m = res['U2_adjusted'], res['U2_main_recomputed']
    res['predictions'] = dict(pc1_adjusted_upper_negative=bool(a['hi975'] < 0), pc2_change_lt_005=bool(abs(a['U2'] - m['U2']) < 0.05))
    res['rule'] = '등록 r44_p1comp (prereg_2026-10-28_addendum5.json): PC1 보정 판 U2 의 97.5% 상한 < 0; PC2 |보정 판 U2 − 주 판 U2| < 0.05.'
    dump(res, RES44 + '/r44_p1comp.json')


if __name__ == '__main__':
    for part in PARTS44:
        print(time.strftime('%H:%M:%S'), '부분 시작', part, flush=True); globals()['part_' + part](); print(time.strftime('%H:%M:%S'), '부분 끝', part, flush=True)
