"""(prereg round36 3 절) E7 30 도메인: 같은 서열 칸의 D 판. 판정 규칙 = 2026-10-19 14:56:51 UTC 부분 동결 원문 sec_3
(Δ², w̄, 비 셋 모두 Spearman ≥ 0.9, D 판 / 맞춘 좌표 판 비의 중앙값 ∈ [0.8, 1.2]).
도메인 = derived/r38/s0/e7_30.json (줄기 336, seed 20261024, 첫 세션 도메인 뺌). r38_s0.py 의 함수를 그대로 불러 쓴다(그 파일은 바꾸지 않음).
계산은 r38_s0.py e7 부분과 같다: 같은 표본(실행마다 N = 100, 실행 쌍마다 프레임 쌍 2,000)에서 D(줄인 계산)와 맞춘 좌표;
칸 = MD T 대 320 K (348·379·413·450; 전체·뒤 절반) + R0 (320 K 반복 2 대 3, 10 가지). 요약의 주 판도 e7 과 같다(온도 칸 + R0 는 도메인·창마다 10 나눔 중앙값 한 칸).
기준 평균 점 = 320 K 접힌 프레임(반복 순서로 이어 앞 600)의 평균 점(r38_md_ref 와 같은 정의; 줄기 도메인에는 r38 MD 기준이 없어 여기서 만든다).
더한 서술 (2026-10-19, E7 5 실패 뒤·이 30 을 돌리기 전에 정함; 판정에 쓰지 않음): 평균 이동 Δμ²/변 을 두 판으로 —
  맞춘 좌표(실행 평균 사이 유클리드², r38_compare 와 같은 추정량: 교차 − (바닥_P + 바닥_Q)/2) 와
  D 판(실행마다 fs.mean_point, 실행 평균 점 사이 D², 128 출발; 같은 추정량) — 그리고 R0 나눔에서 두 판 Δμ² 의 중심.
  또 지수 2 에너지 통계 Δμ²₂ = (C₂ − (B₂_P + B₂_Q)/2)/변 (C₂·B₂ = 실행 쌍마다 프레임 쌍 거리² 의 평균): 맞춘 좌표(전수)에서는 등록된 Δμ² 와 정확히 같고
  (쌍선형 항등식, 퍼짐 항이 지워짐), D 판은 같은 프레임 쌍 2,000 의 D² 로 낸다(평균 점 없이 실제 구조 쌍만).
사용: python r38_e7_30.py <ROOT>
출력: derived/r38/s0/e7_30/<dom>.json, results/r38_e7_30.json (사본 results/_backup/<날짜>/)"""
import sys, os, json, glob, time, itertools
ROOT = sys.argv[1]
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)
OUT2 = OUT + '/e7_30'; os.makedirs(OUT2, exist_ok=True)
SEL = json.load(open(OUT + '/e7_30.json'))['e7_30']
VN = {'d': 'D', 'a2000': 'A2000', 'aall': 'Aall'}


def ref_of(d):
    MD, Qf, L = md_load(d, TT); f = [X[Qf(X) >= QT] for X in MD['320']]
    allf = np.concatenate([x for x in f if len(x)])[:600]
    return torch.as_tensor(fs.mean_point(frames(allf)))


def one(d):
    t0 = time.time(); rng = seeded('r38_e7_30', d); _runs.clear(); _md.clear()
    ref = ref_of(d); res = dict(dom=d, cells={}, nruns={}, diag={})
    for W_ in ('full', 'half'):
        R = {T: dom_runs(d, T, W_, TT)[0] for T in TT}; L = dom_runs(d, '320', W_, TT)[1]; ne = L - 1; res.update(L=L, ne=ne)
        for T in TT: res['nruns'][f'{T}|{W_}'] = len(R[T])
        Ts = [T for T in TT if len(R[T]) >= 2]
        if '320' not in Ts: continue
        E = {T: coords(np.concatenate([X for _, X in R[T]])) for T in Ts}
        A = {T: fs.align_to_reference(frames(np.concatenate([X for _, X in R[T]])), ref).numpy().reshape(len(R[T]) * N, -1) for T in Ts}
        MP = {T: [fs.mean_point(E[T][k * N:(k + 1) * N]) for k in range(len(R[T]))] for T in Ts}      # 실행 평균 점 (D 판)
        AM = {T: [A[T][k * N:(k + 1) * N].mean(0) for k in range(len(R[T]))] for T in Ts}            # 실행 평균 (맞춘 좌표)
        pm, mu = {}, {}
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
            return {**{k: {rp: o[k][g] for g, rp in enumerate(rr)} for k in ('d', 'a2000', 'aall', 'd2', 'a2_2000', 'a2all')}, 'muD': mD, 'muA': mA}
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


def summarize():
    rows = []
    for fn in sorted(glob.glob(OUT2 + '/*.json')):
        j = json.load(open(fn))
        for k, c in j.get('cells', {}).items(): rows.append(dict(dom=j['dom'], cell=k, kind='R0' if k.startswith('null') else 'T', win=k.split('|')[1], c=c))
    def cmp(sub_, ref, qs=('D2', 'wbar', 'ratio')):
        o = {}
        for q in qs:
            x = np.array([r['c']['D'][q] for r in sub_]); y = np.array([r['c'][ref][q] for r in sub_]); ok = np.isfinite(x) & np.isfinite(y) & (np.abs(y) > 1e-12)
            o[q] = dict(spearman=float(spearmanr(x[ok], y[ok]).correlation) if ok.sum() > 2 else None, median_ratio=float(np.median(x[ok] / y[ok])) if ok.any() else None, n=int(ok.sum()))
        return o
    r0 = {}
    for r in rows:
        if r['kind'] == 'R0': r0.setdefault((r['dom'], r['win']), []).append(r)
    med = lambda xs, v, q: float(np.median([x['c'][v][q] for x in xs]))
    r0rows = [dict(dom=k[0], cell=f'R0|{k[1]}', kind='R0med', win=k[1], c={v: {q: med(xs, v, q) for q in ('D2', 'wbar', 'ratio')} for v in ('D', 'A2000', 'Aall')}) for k, xs in r0.items()]
    main = [r for r in rows if r['kind'] == 'T'] + r0rows
    T = [r for r in rows if r['kind'] == 'T']; N0 = [r for r in rows if r['kind'] == 'R0']
    xd = np.array([r['c']['dmu_D'] for r in T]); xa = np.array([r['c']['dmu_A'] for r in T]); okm = np.abs(xa) > 1e-12
    x2 = np.array([r['c']['dmu2_D'] for r in T]); x2a = np.array([r['c']['dmu2_A'] for r in T])
    out = dict(written=time.strftime('%F %T'), n_dom=len({r['dom'] for r in rows}), n_T=len(T), n_R0_split=len(N0),
               main_D_vs_Aall=cmp(main, 'Aall'), main_D_vs_A2000=cmp(main, 'A2000'), T_only_D_vs_Aall=cmp(T, 'Aall'), all_splits_D_vs_Aall=cmp(rows, 'Aall'),
               dmu_T_cells=dict(spearman=float(spearmanr(xd, xa).correlation) if len(T) > 2 else None, median_ratio_D_over_A=float(np.median(xd[okm] / xa[okm])) if okm.any() else None, n=len(T)),
               dmu_R0=dict(D_mean=float(np.mean([r['c']['dmu_D'] for r in N0])) if N0 else None, D_sd=float(np.std([r['c']['dmu_D'] for r in N0])) if N0 else None,
                           A_mean=float(np.mean([r['c']['dmu_A'] for r in N0])) if N0 else None, A_sd=float(np.std([r['c']['dmu_A'] for r in N0])) if N0 else None, n=len(N0)),
               dmu2_T_cells=dict(spearman_vs_dmu_A=float(spearmanr(x2, xa).correlation) if len(T) > 2 else None, median_ratio_vs_dmu_A=float(np.median(x2[okm] / xa[okm])) if okm.any() else None,
                                 identity_aligned_max_abs=float(np.max(np.abs(x2a - xa))) if len(T) else None, n=len(T)),
               dmu2_R0=dict(D_mean=float(np.mean([r['c']['dmu2_D'] for r in N0])) if N0 else None, D_sd=float(np.std([r['c']['dmu2_D'] for r in N0])) if N0 else None, n=len(N0)),
               rule='prereg 3 절: Δ², w̄, 비 셋 모두 Spearman ≥ 0.9, D/맞춘 좌표 비의 중앙값 ∈ [0.8, 1.2]. 주 판 = 온도 칸 + R0(도메인·창마다 10 나눔 중앙값 한 칸), 기준 = 맞춘 좌표 전수. Δμ² 두 판은 서술(판정에 쓰지 않음)')
    m_ = out['main_D_vs_Aall']
    out['verdict'] = '통과' if all(m_[q]['spearman'] is not None and m_[q]['spearman'] >= 0.9 and 0.8 <= m_[q]['median_ratio'] <= 1.2 for q in m_) else '실패'
    jdump(out, RES + '/r38_e7_30.json', bk=True)
    log('E7_30 요약', out['verdict'], {q: (round(m_[q]['spearman'] or 0, 3), round(m_[q]['median_ratio'] or 0, 3)) for q in m_}, 'Δμ²', out['dmu_T_cells'], out['dmu_R0'], 'Δμ²₂', out['dmu2_T_cells'], out['dmu2_R0'])


for d in SEL:
    fo = f'{OUT2}/{d}.json'
    if os.path.exists(fo): continue
    try:
        r = one(d); jdump(r, fo); log('e7_30', d, len(r['cells']), '칸', round(r['sec']), 's')
    except Exception as e:
        traceback.print_exc(); log('ERR e7_30', d, repr(e)[:200])
summarize(); log('E7_30_DONE')
