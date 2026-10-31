"""(등록 = Drive results/prereg_2026-10-30_addendum10.json "r47_ties") 동점 순위 점검: 동점이 가장 많고 구간 끝이 0 에 붙은 곳을 평균 순위로 다시 낸다.
부분 (서로 독립; 부분 × 순위 방식마다 프로세스 하나로 병렬):
  ties      : 동일성·순위 변수의 동점 구조(서로 다른 값 수, 동점 몫 T = Σ(t³ − t)/(n³ − n)) — C, H-풀(첫 연구), 비선형 스크립트의 H 행; ANM 예측(13·15 Å, GNM)의 도메인 안 동점 수. 결과 변수와 무관한 자료 기술.
  p3 MODE   : 첫 연구 1b = r38_p3.analyse(rows, 'main') (seed 20261028, 상과 × 두 도메인 Exp(1) 곱 2,000) 와
              상과 무게 구간 = r38_p3_explore2 의 real['main'] (seed 20261051 한 흐름, 반복마다 상과 v 다음 도메인 w, 1,000 회, 상과만 무게 판).
  nl MODE   : r38_c_nl 의 H 행(load(P3, anm_pairs_hpool, 'main'))에서 선형 두 변수 판(control + r_anm)과 선형 8 변수 판 (seed 20261053, 2,000). GBM(비선형 판)은 다시 적합하지 않음.
  anm MODE  : r43_sens s3 의 H 행(pair_rows(P3); 쌍별 r_ANM 판은 r43_sens 그대로 — 연속 예측의 Spearman)에서 C-2 = analyse_sf(rows, 'main', z = (control, r_판), seed 20261077),
              판 = ANM 10·13·15 Å, GNM 7.3 Å.
  combine   : 저장값(results/r38_p3.json, r38_p3_explore2.json, r38_c_nl.json, r43_sens.json, r46_align_fix1.json 의 C 점검)과 두 판을 맞대어 results/r47_ties.json.
MODE = argsort (등록 판의 순위: np.argsort 두 번) | average (scipy.stats.rankdata method='average'). 순위 함수는 원래 함수가 부르는 전역 ranks 를 바꾼다.
같은 런타임에서 두 판을 낸다. 잼: 추정치마다 점추정·구간 끝의 |평균 순위 − argsort|, |평균 순위 − 저장값|, |argsort − 저장값|의 최댓값, 그리고 0 에 대한 구간 위치
(위·아래·걸침)가 저장값·argsort·평균 순위 사이에서 바뀌는지.
사용: python r47_ties.py <ROOT> <부분> [MODE]   출력: derived/r47/<부분>_<MODE>.json, results/r47_ties.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
import sys, json, glob, time, traceback
import numpy as np
ROOT = sys.argv[1]; PART = sys.argv[2]; MODE = sys.argv[3] if len(sys.argv) > 3 else 'none'
RES, DER = ROOT + '/results', ROOT + '/derived'; OUT47 = DER + '/r47'; os.makedirs(OUT47, exist_ok=True)
P3D, CD = DER + '/r38/p3', DER + '/r38/c'


def log(*a): print(time.strftime('%H:%M:%S'), PART, MODE, *a, flush=True)
def rk_argsort(v): return np.argsort(np.argsort(v)).astype(float)
def rk_average(v):
    from scipy.stats import rankdata
    return rankdata(v, method='average').astype(float)
RK = dict(argsort=rk_argsort, average=rk_average)
def triple(d, est='est', lo='lo', hi='hi'): return None if d is None else dict(est=float(d[est]), lo=float(d[lo]), hi=float(d[hi]))
def save(o): json.dump(dict(o, written=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), part=PART, mode=MODE), open(f'{OUT47}/{PART}_{MODE}.json', 'x'), ensure_ascii=False, indent=1)


def load_ns(path, cut, argv):
    """원래 스크립트의 앞부분을 새 이름 공간에 실행(원래 결과 파일에 쓰는 부분은 실행하지 않음)."""
    src = open(path).read(); src = src[:src.index(cut)] if cut else src
    g = {'__name__': 'r47_' + os.path.basename(path)}; old = sys.argv; sys.argv = argv
    try: exec(src, g)
    finally: sys.argv = old
    return g


def part_p3():
    g = load_ns(ROOT + '/code/r38_p3.py', "ALL = ['pairs'", ['r38_p3', ROOT, 'none', '1']); g['ranks'] = RK[MODE]
    rows = [json.loads(l) for l in open(g['OUT'] + '/pairs_out.jsonl')]
    m = g['analyse'](rows, 'main')   # 등록 1b: seed 20261028, NB 2,000, 상과 × 두 도메인
    # r38_p3_explore2 real['main']: seed 20261051 의 첫 판, 반복마다 v(상과) 다음 w(도메인), 1,000 회
    rr = [r for r in rows if r['main'].get('r_dyn') is not None]; rkf = g['ranks']
    Y = rkf(np.array([r['main']['r_dyn'] for r in rr])); X = rkf(np.array([r['main']['ident'] for r in rr])); Z = [rkf(np.array([r['main']['control'] for r in rr]))]
    H = np.array([r['H'] for r in rr]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); base = 1.0 / np.bincount(si)[si]
    doms = sorted({r['ab'][0] for r in rr} | {r['ab'][1] for r in rr}); di = {d: i for i, d in enumerate(doms)}
    ia = np.array([di[r['ab'][0]] for r in rr]); ib = np.array([di[r['ab'][1]] for r in rr])
    rng = np.random.default_rng(20261051); sf = []
    for _ in range(1000):
        v = rng.exponential(1, len(sfl)); w = rng.exponential(1, len(doms)); sf.append(g['pcorr'](Y, X, Z, base * v[si]))
    est = g['pcorr'](Y, X, Z, base)
    save(dict(first_1b=dict(n=m['n'], n_sf=m['n_sf'], **triple(m)), first_sf_only=dict(n=len(rr), est=float(est), lo=float(np.percentile(sf, 2.5)), hi=float(np.percentile(sf, 97.5)))))
    log('p3 끝', round(m['est'], 6))


def part_nl():
    g = load_ns(ROOT + '/code/r38_c_nl.py', "def nl_resid", ['r38_c_nl', ROOT, 'none']); rkf = RK[MODE]; pcorr = g['pcorr']; ZK = g['ZK']
    rows = g['load'](P3D, CD + '/anm_pairs_hpool.jsonl', 'main')
    H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); base = 1.0 / np.bincount(si)[si]
    Y = rkf(np.array([r['y'] for r in rows])); X = rkf(np.array([r['x'] for r in rows])); Zm = np.array([[r[k] for k in ZK] for r in rows], float); Zr = [rkf(Zm[:, j]) for j in range(Zm.shape[1])]
    Z2 = [Zr[ZK.index('control')], Zr[ZK.index('r_anm')]]
    est = dict(lin2=pcorr(Y, X, Z2, base), linall=pcorr(Y, X, Zr, base))
    rng = np.random.default_rng(20261053); bs = dict(lin2=[], linall=[])
    for _ in range(2000):   # r38_c_nl.analyse 와 같은 무게 흐름(반복마다 상과 Exp(1) 한 번; 비선형 판은 난수를 쓰지 않음)
        w = base * rng.exponential(1, len(sfl))[si]; bs['lin2'].append(pcorr(Y, X, Z2, w)); bs['linall'].append(pcorr(Y, X, Zr, w))
    save({f'H_{k}': dict(n=len(rows), n_sf=len(sfl), est=float(est[k]), lo=float(np.percentile(bs[k], 2.5)), hi=float(np.percentile(bs[k], 97.5))) for k in est})
    log('nl 끝', {k: round(v, 6) for k, v in est.items()})


def part_anm():
    g = load_ns(ROOT + '/code/r43_sens.py', "if __name__ == '__main__':", ['r43_sens', ROOT, 'none', '1']); g['ranks'] = RK[MODE]
    rows = g['pair_rows'](P3D); o = {}
    for v in g['VARS']:
        m = g['analyse_sf'](rows, 'main', z=('control', 'r_' + v), seed=20261077); o[f'H_C2_{v}'] = dict(n=m['n'], n_sf=m['n_sf'], **triple(m))
    save(o); log('anm 끝', {k: round(x['est'], 6) for k, x in o.items()})


def tie_stats(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; n = len(v); _, t = np.unique(v, return_counts=True)
    return dict(n=int(n), distinct=int(len(t)), max_tie=int(t.max()) if n else 0, T=float(((t ** 3 - t).sum()) / (n ** 3 - n)) if n > 1 else None)


def part_ties():
    o = {}
    for nm, d in (('C', CD), ('H_first', P3D)):
        rr = [json.loads(l)['main'] for l in open(d + '/pairs_out.jsonl')]; rr = [m for m in rr if m.get('r_dyn') is not None]
        o[nm] = {k: tie_stats([m.get(k, np.nan) if m.get(k) is not None else np.nan for m in rr]) for k in ('ident', 'r_dyn', 'control')}
    g = load_ns(ROOT + '/code/r38_c_nl.py', "def nl_resid", ['r38_c_nl', ROOT, 'none']); hr = g['load'](P3D, CD + '/anm_pairs_hpool.jsonl', 'main')
    o['H_nl_rows'] = dict(x=tie_stats([r['x'] for r in hr]), y=tie_stats([r['y'] for r in hr]), **{k: tie_stats([r[k] for r in hr]) for k in g['ZK']})
    doms = json.load(open(P3D + '/pairs.json'))['domains']; cnt = {v: [0, 0] for v in ('anm13', 'anm15', 'gnm')}
    for dd in doms:
        f = f'{DER}/r43/anm/{dd}.npz'
        if not os.path.exists(f): continue
        z = np.load(f)
        for v in cnt:
            if v in z.files: a = z[v]; cnt[v][0] += 1; cnt[v][1] += int(len(a) - len(np.unique(a)) > 0)
    o['anm_pred_domains_with_ties'] = {v: dict(n_domains=c[0], with_ties=c[1]) for v, c in cnt.items()}
    import platform, numpy, scipy
    flags = sorted(set(open('/proc/cpuinfo').read().split('flags')[1].split('\n')[0].replace(':', ' ').split())) if os.path.exists('/proc/cpuinfo') else None
    o['env'] = dict(python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__, avx512f=bool(flags and 'avx512f' in flags), avx2=bool(flags and 'avx2' in flags),
                    cpu_model=next((l.split(':', 1)[1].strip() for l in open('/proc/cpuinfo') if l.startswith('model name')), None) if os.path.exists('/proc/cpuinfo') else None)
    save(o); log('ties 끝', o['H_first']['ident'])


def zside(t): return None if t is None else ('above' if t['lo'] > 0 else 'below' if t['hi'] < 0 else 'straddle')
def dmax(a, b): return None if a is None or b is None else float(max(abs(a[k] - b[k]) for k in ('est', 'lo', 'hi')))


def part_combine():
    J = lambda f: json.load(open(f'{OUT47}/{f}.json'))
    av = {**J('p3_average'), **J('nl_average'), **J('anm_average')}; ag = {**J('p3_argsort'), **J('nl_argsort'), **J('anm_argsort')}
    st = {}
    try: st['first_1b'] = triple(json.load(open(RES + '/r38_p3.json'))['main'])
    except Exception: st['first_1b'] = None
    try:
        e = json.load(open(RES + '/r38_p3_explore2.json'))['real_sf_only']['main']; st['first_sf_only'] = dict(est=e['est'], lo=e['ci_sf'][0], hi=e['ci_sf'][1])
    except Exception: st['first_sf_only'] = None
    try:
        e = json.load(open(RES + '/r38_c_nl.json'))['H_main']
        for k in ('lin2', 'linall'): st[f'H_{k}'] = dict(est=e['est'][k], lo=e['ci95'][k][0], hi=e['ci95'][k][1])
    except Exception: st['H_lin2'] = st['H_linall'] = None
    try:
        e = json.load(open(RES + '/r43_sens.json'))['s3']['H']['judge']
        for v in ('anm10', 'anm13', 'anm15', 'gnm'): st[f'H_C2_{v}'] = triple(e[v]['main'])
    except Exception:
        for v in ('anm10', 'anm13', 'anm15', 'gnm'): st.setdefault(f'H_C2_{v}', None)
    rows = {}
    for k in sorted(av):
        a, q, s = triple(av[k]), triple(ag[k]), st.get(k)
        sides = dict(stored=zside(s), argsort=zside(q), average=zside(a))
        rows[k] = dict(n=av[k].get('n'), stored=s, argsort=q, average=a, d_avg_argsort=dmax(a, q), d_avg_stored=dmax(a, s), d_argsort_stored=dmax(q, s), sides=sides,
                       crossing=bool(len({x for x in sides.values() if x is not None}) > 1))
    fx = json.load(open(RES + '/r46_align_fix1.json'))['check']; c2a, c2b = fx['check_2a'], fx['check_2b']
    for c in ('C1', 'C2'):
        reg, q, a = triple(c2a[f'registered_{c}']), triple(c2a[f'argsort_{c}']), triple(c2b[f'avg_{c}'])
        rows[f'C_{c}'] = dict(n=c2a[f'argsort_{c}']['n'], stored=reg, argsort=q, average=a, d_avg_argsort=dmax(a, q), d_avg_stored=dmax(a, reg), d_argsort_stored=dmax(q, reg),
                              sides=dict(stored=zside(reg), argsort=zside(q), average=zside(a)), crossing=bool(len({zside(reg), zside(q), zside(a)}) > 1), source='r46_align_fix1 check')
    mx = lambda key: max((r[key] for r in rows.values() if r[key] is not None), default=None)
    out = dict(written=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), estimates=rows, ties=J('ties_none'),
               max_d_avg_stored=mx('d_avg_stored'), max_d_avg_argsort=mx('d_avg_argsort'), max_d_argsort_stored=mx('d_argsort_stored'),
               crossings=[k for k, r in rows.items() if r['crossing']],
               rule='등록 r47_ties: 원고 본문은 등록된 argsort 값(저장값)을 그대로 둔다. 방법 절에 평균 순위와의 최대 |차|(max_d_avg_stored)와 동점 순서에 따른 최대 |차|를 "at most …" 로 적는다. '
                    '0 에 대한 구간 위치가 저장값·argsort·평균 순위 사이에서 바뀌는 추정치(crossings)가 있으면 원고의 그 자리에 두 판(저장값과 평균 순위)을 함께 보고한다.')
    json.dump(out, open(RES + '/r47_ties.json', 'x'), ensure_ascii=False, indent=1); log('combine 끝', out['max_d_avg_stored'], out['crossings'])


if __name__ == '__main__':
    try: globals()['part_' + PART](); log('R47_PART_DONE')
    except Exception as e: traceback.print_exc(); log('ERR', repr(e)[:300])
