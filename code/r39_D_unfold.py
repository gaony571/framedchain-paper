"""(round39 E6 D 판; 2026-10-20 UTC 사용자 결정 "(가)") 같은 서열 접힘 대 풀림을 D 로: 줄기 336 의 mdCATH 413·450 K.
반복마다 접힘(Q ≥ qf)·풀림(Q < qu) 무리, 문턱 셋 (0.7, 0.3)·(0.75, 0.4)·(0.8, 0.5) (시연 설계 E6 의 문턱 민감도). 창 = 전체(round37 R7 과 같음).
실행 = (반복, 상태) 무리에서 N = 100 프레임(r38 공통 코드 규칙; round37 보조는 무리 ≥ 40 이었으므로 도메인 수가 줄 수 있다). 상태마다 무리 ≥ 2 인 칸만.
교차 = 접힘 반복 i 대 풀림 반복 j (i ≠ j, R7 과 같음), 바닥 = 같은 상태의 서로 다른 반복 쌍. 실행 쌍마다 프레임 쌍 2,000 (줄인 계산, r38_s0).
Q 는 320 K 첫 프레임의 천연 접촉(r38 공통 md_load). 같은 표본의 맞춘 좌표 값(접힘 프레임 앞 600 의 평균 점에 맞춤)도 함께 낸다.
사용: python r39_D_unfold.py <ROOT> [run,sum] 출력: derived/r39/unfold_D/<dom>.json, results/r39_unfold_D.json"""
import sys, os, json
ROOT = sys.argv[1]; PART = sys.argv[2].split(',') if len(sys.argv) > 2 else ['run', 'sum']
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)
OUTU = ROOT + '/derived/r39/unfold_D'; os.makedirs(OUTU, exist_ok=True)
STEM = sorted(os.path.basename(f)[:-4] for f in glob.glob(ROOT + '/derived/mdcath_5T/*.npz'))
THR = ((0.7, 0.3), (0.75, 0.4), (0.8, 0.5))


def one_u(d):
    t0 = time.time(); rng = seeded('r39_unfold_D', d); _runs.clear(); _md.clear()
    MD, Qf, L = md_load(d, ('413', '450')); ne = L - 1; res = dict(dom=d, L=L, ne=ne, cells={})
    for T in ('413', '450'):
        if T not in MD: continue
        QX = [Qf(X) for X in MD[T]]
        for qf, qu in THR:
            key = f'{T}|{qf}-{qu}'; G = {'f': [], 'u': []}
            for r, X in enumerate(MD[T]):
                for s, ix in (('f', np.where(QX[r] >= qf)[0]), ('u', np.where(QX[r] < qu)[0])):
                    if len(ix) >= N: G[s].append((r, X[np.sort(rng.choice(ix, N, replace=False))]))
            rf = [r for r, _ in G['f']]; ru = [r for r, _ in G['u']]; res['cells'][key] = dict(reps_f=rf, reps_u=ru)
            crp = [(i, j) for i in range(len(rf)) for j in range(len(ru)) if rf[i] != ru[j]]
            if len(rf) < 2 or len(ru) < 2 or not crp: continue
            E = {s: coords(np.concatenate([X for _, X in G[s]])) for s in G}
            fall = np.concatenate([X for _, X in G['f']]); ref = torch.as_tensor(fs.mean_point(frames(fall[:600])))
            A = {s: fs.align_to_reference(frames(np.concatenate([X for _, X in G[s]])), ref).numpy().reshape(len(G[s]) * N, -1) for s in G}
            def go(a, b, rr):
                o = runpair_means(MeanStart(E[a], E[b]), rr, rng)
                return [dict(d=o['d'][g], a=float(np.sqrt(np.maximum(cost(A[a][x * N:(x + 1) * N], A[b][y * N:(y + 1) * N]), 0)).mean())) for g, (x, y) in enumerate(rr)]
            cr = go('f', 'u', crp); ff = go('f', 'f', list(itertools.combinations(range(len(rf)), 2))); uu = go('u', 'u', list(itertools.combinations(range(len(ru)), 2)))
            m = lambda xs, k: float(np.mean([x[k] for x in xs]))
            res['cells'][key].update(D=estats(m(cr, 'd'), m(ff, 'd'), m(uu, 'd'), ne), A=estats(m(cr, 'a'), m(ff, 'a'), m(uu, 'a'), ne), n_cross=len(crp))
    res['sec'] = time.time() - t0
    return res


def part_run():
    for d in STEM:
        fo = f'{OUTU}/{d}.json'
        if os.path.exists(fo): continue
        try:
            r = one_u(d); jdump(r, fo); n_ = sum('D' in c for c in r['cells'].values())
            if n_: log('r39_unfold_D', d, n_, '칸', round(r['sec']), 's')
        except Exception as e:
            traceback.print_exc(); log('ERR r39_unfold_D', d, repr(e)[:200])
    log('R39_UNFOLD_RUN_DONE')


def part_sum():
    J = [json.load(open(f)) for f in sorted(glob.glob(OUTU + '/*.json'))]; out = {}
    for T in ('413', '450'):
        for qf, qu in THR:
            key = f'{T}|{qf}-{qu}'; v = [j['cells'][key]['D']['d2n'] for j in J if 'D' in j['cells'].get(key, {})]; a = [j['cells'][key]['A']['d2n'] for j in J if 'D' in j['cells'].get(key, {})]
            if v: out[key] = dict(n=len(v), D_median=float(np.median(v)), D_iqr=[float(x) for x in np.percentile(v, [25, 75])], A_median=float(np.median(a)), ratio_D_over_A=float(np.median(np.array(v) / np.array(a))))
    jdump(dict(written=time.strftime('%F %T'), n_done=len(J), cells=out, round37_aligned={'413': 0.359, '450': 0.403},
               rule='round39 E6 D 판: 접힘 대 풀림 Δ²_D/√변, 문턱 셋, 무리 N = 100. 띠(상과 중앙값 95% 구간)와 위치 안정 판정은 round39 CPU 분석에서.'), RES + '/r39_unfold_D.json', bk=True)
    log('R39_UNFOLD_DONE', {k: (v['n'], round(v['D_median'], 3)) for k, v in out.items()})


if __name__ == '__main__':
    for part in PART:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
