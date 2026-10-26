"""(round39 D 판; 2026-10-20 UTC 사용자 결정 "(가) 줄기 336 전체 D, 맨 마지막에") 줄기 derived/mdcath_5T 336 도메인의 같은 서열 칸을 D 로.
계산은 r38_e7_30.one() 과 같은 식이다(실행마다 N = 100, 실행 쌍마다 프레임 쌍 2,000, 창 = 전체·뒤 절반, 칸 = MD T 대 320 K + R0 2:3 나눔 10, 같은 cell_of).
더해서 실행 쌍 표(rp)를 저장한다: 창 × 온도 쌍마다 [실행 a, 실행 b, d, d², a(2,000), a(전수), a²(2,000), a²(전수), 실행 평균 점 D², 실행 평균 ²].
round39 의 부분집합 통계(E3 의 n 별 R0 띠, E4 의 공유 · 표적 3 나눔, E0′ 의 C·B)는 이 표에서 CPU 로 낸다.
seed 꼬리표 = 'r39_D'. check 부분은 꼬리표를 'r38_e7_30' 으로 바꿔 E7 30 출력(derived/r38/s0/e7_30/<dom>.json)과 칸이 같은지 본다.
사용: python r39_D_stem.py <ROOT> run,sum | check <도메인> 출력: derived/r39/D/<dom>.json, results/r39_D_stem.json"""
import sys, os, json
ROOT = sys.argv[1]; PART = sys.argv[2].split(','); ARG = sys.argv[3] if len(sys.argv) > 3 else None
_e = open(ROOT + '/code/r38_e7_30.py').read(); _e = _e[:_e.index('\nfor d in SEL:')]
sys.argv = ['r38_e7_30', ROOT]; exec(_e)
OUT4 = ROOT + '/derived/r39/D'; os.makedirs(OUT4, exist_ok=True)
STEM = sorted(os.path.basename(f)[:-4] for f in glob.glob(ROOT + '/derived/mdcath_5T/*.npz'))
KEYS = ('d', 'd2', 'a2000', 'aall', 'a2_2000', 'a2all', 'muD', 'muA')


def one_rp(d, tag='r39_D'):
    """r38_e7_30.one() 과 같은 계산(같은 호출 순서) + 실행 쌍 표."""
    t0 = time.time(); rng = seeded(tag, d); _runs.clear(); _md.clear()
    ref = ref_of(d); res = dict(dom=d, cells={}, nruns={}, diag={}, rp={}, rp_keys=list(KEYS))
    for W_ in ('full', 'half'):
        R = {T: dom_runs(d, T, W_, TT)[0] for T in TT}; L = dom_runs(d, '320', W_, TT)[1]; ne = L - 1; res.update(L=L, ne=ne)
        for T in TT: res['nruns'][f'{T}|{W_}'] = len(R[T])
        res['rp'][W_] = dict(runs={T: [r for r, _ in R[T]] for T in TT})
        Ts = [T for T in TT if len(R[T]) >= 2]
        if '320' not in Ts: continue
        E = {T: coords(np.concatenate([X for _, X in R[T]])) for T in Ts}
        A = {T: fs.align_to_reference(frames(np.concatenate([X for _, X in R[T]])), ref).numpy().reshape(len(R[T]) * N, -1) for T in Ts}
        MP = {T: [fs.mean_point(E[T][k * N:(k + 1) * N]) for k in range(len(R[T]))] for T in Ts}
        AM = {T: [A[T][k * N:(k + 1) * N].mean(0) for k in range(len(R[T]))] for T in Ts}
        pm = {}
        def go(Ta, Tb, rr):
            ms = MeanStart(E[Ta], E[Tb]); o = runpair_means(ms, rr, rng)
            ea = np.sqrt(((A[Ta][o['IA']] - A[Tb][o['IB']]) ** 2).sum(1))
            o['a2000'] = [float(ea[o['GR'] == g].mean()) for g in range(len(rr))]
            o['aall'] = [float(np.sqrt(cost(A[Ta][a * N:(a + 1) * N], A[Tb][b * N:(b + 1) * N])).mean()) for a, b in rr]
            res['diag'][f'{Ta}-{Tb}|{W_}'] = dict(widened=o['widened'], kmax=o['kmax'], n_pairs=len(o['dv']), S_mean=ms.S_mean)
            mD = {rp: float(fs.sq_distance(MP[Ta][rp[0]], MP[Tb][rp[1]])[0]) for rp in rr}
            mA = {rp: float(((AM[Ta][rp[0]] - AM[Tb][rp[1]]) ** 2).sum()) for rp in rr}
            o['d2'] = [float((o['dv'][o['GR'] == g] ** 2).mean()) for g in range(len(rr))]
            o['a2_2000'] = [float((ea[o['GR'] == g] ** 2).mean()) for g in range(len(rr))]
            o['a2all'] = [float(cost(A[Ta][a * N:(a + 1) * N], A[Tb][b * N:(b + 1) * N]).mean()) for a, b in rr]
            out = {**{k: {rp: o[k][g] for g, rp in enumerate(rr)} for k in ('d', 'a2000', 'aall', 'd2', 'a2_2000', 'a2all')}, 'muD': mD, 'muA': mA}
            res['rp'][W_][f'{Ta}-{Tb}'] = [[int(rp[0]), int(rp[1])] + [out[k][rp] for k in KEYS] for rp in rr]
            return out
        def cell_of(cr, fp, fq):
            c = {VN[v]: estats(np.mean(cr[v]), np.mean(fp[v]), np.mean(fq[v]), ne) for v in ('d', 'a2000', 'aall')}
            for v, nm in (('muD', 'dmu_D'), ('muA', 'dmu_A'), ('d2', 'dmu2_D'), ('a2_2000', 'dmu2_A2000'), ('a2all', 'dmu2_A')): c[nm] = float((np.mean(cr[v]) - (np.mean(fp[v]) + np.mean(fq[v])) / 2) / ne)
            return c
        n320 = len(R['320']); pm[('320', '320')] = go('320', '320', list(itertools.combinations(range(n320), 2)))
        for T in Ts[1:]:
            nT = len(R[T]); pm[(T, T)] = go(T, T, list(itertools.combinations(range(nT), 2)))
            pm[('320', T)] = go('320', T, list(itertools.product(range(n320), range(nT))))
            sv = lambda key: {v: list(pm[key][v].values()) for v in ('d', 'a2000', 'aall', 'muD', 'muA', 'd2', 'a2_2000', 'a2all')}
            cell = cell_of(sv(('320', T)), sv(('320', '320')), sv((T, T))); cell.update(nP=n320, nQ=nT); res['cells'][f'MD{T}-MD320|{W_}'] = cell
        if n320 == 5:
            f = lambda v, a, b: pm[('320', '320')][v][(min(a, b), max(a, b))]
            for i, g2 in enumerate(itertools.combinations(range(5), 2)):
                g3 = [k for k in range(5) if k not in g2]
                sv = lambda pairs: {v: [f(v, a, b) for a, b in pairs] for v in ('d', 'a2000', 'aall', 'muD', 'muA', 'd2', 'a2_2000', 'a2all')}
                res['cells'][f'null320|{W_}|{i}'] = cell_of(sv([(a, b) for a in g2 for b in g3]), sv([g2]), sv(list(itertools.combinations(g3, 2))))
    res['sec'] = time.time() - t0
    return res


def part_check():
    """E7 30 의 한 도메인을 seed 꼬리표 'r38_e7_30' 으로 다시 내어 칸이 같은지(상대 차) 본다. 출력은 쓰지 않는다."""
    d = ARG; r = one_rp(d, tag='r38_e7_30'); old = json.load(open(f'{OUT2}/{d}.json'))['cells']; worst = 0.0; n = 0
    for k, c in old.items():
        for v, x in c.items():
            ys = x if isinstance(x, dict) else {'_': x}; zs = r['cells'][k][v] if isinstance(x, dict) else {'_': r['cells'][k][v]}
            for q, y in ys.items():
                if isinstance(y, (int, float)) and np.isfinite(y): worst = max(worst, abs(zs[q] - y) / max(abs(y), 1e-9)); n += 1
    print('R39_CHECK', d, '칸', len(old), '값', n, '최대 상대 차', worst, '실행 쌍 표', {W: len(v) - 1 for W, v in r['rp'].items()}, round(r['sec']), 's', flush=True)


def part_run():
    for d in STEM:
        fo = f'{OUT4}/{d}.json'
        if os.path.exists(fo): continue
        try:
            r = one_rp(d); jdump(r, fo); log('r39_D', d, len(r['cells']), '칸', round(r['sec']), 's')
        except Exception as e:
            traceback.print_exc(); log('ERR r39_D', d, repr(e)[:200])
    log('R39_D_RUN_DONE')


def part_sum():
    J = [json.load(open(f)) for f in sorted(glob.glob(OUT4 + '/*.json'))]; cnt = {}
    for j in J:
        for k in j['cells']:
            if k.startswith('MD'): cnt[k] = cnt.get(k, 0) + 1
    jdump(dict(written=time.strftime('%F %T'), n_stem=len(STEM), n_done=len(J), n_with_cells=sum(1 for j in J if any(k.startswith('MD') for k in j['cells'])), cells=cnt,
               rule='round39 D 판(사용자 결정 2026-10-20): r38_e7_30.one() 과 같은 계산 + 실행 쌍 표. 판정·분석은 round39 CPU 스크립트에서.'), RES + '/r39_D_stem.json', bk=True)
    log('R39_D_DONE', len(J), cnt)


if __name__ == '__main__':
    for part in PART:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
