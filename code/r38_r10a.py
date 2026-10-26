"""(prereg round36 2 절, R10a′) 자연 삽입·결실 쌍(results/r37_trunk_pairs.json 의 r10a_prime: 주 19 = 동일성 ≥ 0.25, 민감도 40 ⊇ 주)의 앙상블 D.
계산은 r38_b6.one() 그대로(320 K 전체 창, 접힌 프레임 Q ≥ 0.75, 반복마다 N = 100, 프레임 쌍 2,000, 줄인 계산; Δ²_D, 짝·빈칸 몫, D_matched, π̄·TM-align 고정 대응,
정적 BnB). 더하는 것: 같은 서열 자기 비교의 빈칸 몫(바닥 실행 쌍의 빈칸 평균; one() 안의 runpair_means 호출을 그대로 받아 씀).
stats (prereg 2 절 예측, 등록 2026-10-20 보충 = Drive JSON round36_r10a_2026-10-20):
  1 (서술) 관행 열(TM-align 짝 잔기 고정 대응) Δ²/√n̄ 이 R0 띠(r38_s0.R0_BAND) 안이고 Δ²_D/√n̄ 이 띠 위인 주 판 쌍의 비율 (≥ 70% 예측).
  2 (시험) 빈칸 초과 = 교차 빈칸 평균 − (A 자기 + B 자기)/2 (D² 단위) 가 0.75 × 삽입 잔기 수와 맞는가: 쌍별 상대 오차, 주 판 중앙 상대 오차 < 10%.
  3 (판정, 2차) 짝 몫(변당)의 상동체 회귀선(B6 같은 상과 쌍, 짝 몫 ~ TM 긴 쪽, 최소제곱) 잔차의 상과 중앙값: TOST(여유 = 회귀 잔차 MAD 의 0.5 배,
    상과 부트스트랩 2,000 의 90% 구간이 ±여유 안 → 동등). 구간이 0 을 빼고 여유 밖 → 다르다. 그 밖 → 가를 수 없다(검출 가능한 최소 차 [x]).
    민감도: 동일성 공변량 판(B6 상동체 쌍의 동일성 = pilot_mdcath_homolog_pairs.tsv), TM 보정 판(TM × 긴 쪽 길이 / π̄ 짝 수), 민감도 40 쌍 판.
고침 1 (2026-10-20 UTC, 앙상블 계산 전; Drive JSON round36_r10a_2026-10-20.amend1): 등록의 "틀린 쌍은 정렬이 빈칸을 어디에 두었는지 보고한다" 를 구현 —
  red_pairs 를 같은 계산 그대로 다시 써서 변별 빈칸 질량을 남기고(gap_regions: π̄ 기준 끝 겹침·건너뜀 창·중심), 사후 서술 판(끝 겹침을 뺀 초과)과
  기준 회귀선에서 시험 쌍과 겹치는 B6 쌍을 뺀 민감도를 더함. 등록 판정 규칙(예측 1·2·3)은 그대로.
사용: python r38_r10a.py <ROOT> ens|stats 출력: derived/r38/r10a/<a>__<b>.json, results/r38_r10a.json"""
import sys, os, json
ROOT = sys.argv[1]; PARTS_R = sys.argv[2].split(',')
_b = open(ROOT + '/code/r38_b6.py').read(); _b = _b[:_b.index("\ndef part_ens():")]
sys.argv = ['r38_b6', ROOT, 'none']; exec(_b)
OUTR = ROOT + '/derived/r38/r10a'; os.makedirs(OUTR, exist_ok=True)
_rpm0 = runpair_means; _CAP = []; GL = dict(on=False, rec=None)
def runpair_means(ms, rr, rng, want_gap=False):
    GL['on'] = want_gap; o = _rpm0(ms, rr, rng, want_gap); GL['on'] = False
    if want_gap: o['loc'] = GL['rec']; o['pi'] = [x.tolist() for x in np.nonzero(ms.match.cpu().numpy())]
    _CAP.append(o); return o


def red_pairs(ms, ia, ib, k0=None, iters=10, want_gap=False, membudget=3e9):
    """(2026-10-20 고침 1, 빈칸 자리) r38_b6.red_pairs 와 같은 계산·같은 연산 순서. GL['on'] 이면 쌍마다 마지막 풀이의 짝 없는 변 표시를 남겨
    변별 빈칸 질량 평균(MA: A 변, MB: B 변; 합 = gv 평균)과 짝 없는 빈도(UA, UB)를 GL['rec'] 에 둔다. 값(d, 빈칸 몫, S)은 바꾸지 않는다."""
    nl = RED['n_local']; k0 = RED['k0'] if k0 is None else k0
    EA, EB = ms.EA, ms.EB; n, m = EA.shape[1], EB.shape[1]
    ia = torch.as_tensor(np.asarray(ia), device=EA.device); ib = torch.as_tensor(np.asarray(ib), device=EA.device)
    nA, nB = (EA ** 2).sum((-1, -2)), (EB ** 2).sum((-1, -2))
    S = EA.new_zeros(len(ia)); G = EA.new_zeros(len(ia)) if want_gap else None; kused = np.zeros(len(ia), int); widened = 0
    rec = GL['on'] and want_gap
    if rec: UA = torch.zeros(len(ia), n, dtype=torch.bool, device=EA.device); UB = torch.zeros(len(ia), m, dtype=torch.bool, device=EA.device)
    todo = torch.arange(len(ia), device=EA.device); kk = k0
    while len(todo):
        lo, w = band_windows(ms.center, m, kk); step = max(1, int(membudget / (n * w * 9 * 8 * 4 * nl)))
        edge_all = torch.zeros(len(todo), dtype=torch.bool, device=EA.device)
        for c in range(0, len(todo), step):
            sel = todo[c:c + step]; a, b = EA[ia[sel]], EB[ib[sel]]
            sc, _, M, edge = align_band(a, b, lo, w, ms.starts(ia[sel], ib[sel], n_local=nl), iters=iters)
            S[sel] = sc; edge_all[c:c + step] = edge
            if want_gap: G[sel] = (nA[ia[sel]] * (~M.any(2))).sum(1) + (nB[ib[sel]] * (~M.any(1))).sum(1)
            if rec: UA[sel] = ~M.any(2); UB[sel] = ~M.any(1)
        kused[todo.cpu().numpy()] = kk
        if w >= m: break
        todo = todo[edge_all]
        if kk == k0: widened = int(len(todo))
        kk *= 2
    d2 = (sq_norm(EA)[ia] + sq_norm(EB)[ib] - 2 * S).clamp_min(0)
    if rec: GL['rec'] = dict(MA=(UA * nA[ia]).mean(0).cpu().numpy().tolist(), MB=(UB * nB[ib]).mean(0).cpu().numpy().tolist(), UA=UA.double().mean(0).cpu().numpy().tolist(), UB=UB.double().mean(0).cpu().numpy().tolist())
    return d2.sqrt().cpu().numpy(), (G.cpu().numpy() if want_gap else None), S.cpu().numpy(), kused, widened


def gap_regions(cr, fa, fb, w=3):
    """빈칸 초과를 변 자리로 나눈다. 자리는 π̄ (두 평균 점 정렬)의 짝에서: 끝 겹침 = 첫·마지막 짝 바깥 변, 건너뜀 창 = 이웃한 두 짝 (i1, j1)·(i2, j2) 의 걸음이 (1, 1) 이 아닌 곳의 [i1 − w, i2 + w]·[j1 − w, j2 + w], 나머지 = 짝지어진 중심.
    변별 초과 = 교차 빈칸 질량 − 같은 쪽 자기 비교 빈칸 질량(자기 비교의 두 쪽 평균). 합은 등록된 빈칸 초과와 같다."""
    pa, pb = np.array(cr['pi'][0]), np.array(cr['pi'][1]); o = np.argsort(pa); pa, pb = pa[o], pb[o]
    sA = (np.array(fa['loc']['MA']) + np.array(fa['loc']['MB'])) / 2; sB = (np.array(fb['loc']['MA']) + np.array(fb['loc']['MB'])) / 2
    eA = np.array(cr['loc']['MA']) - sA; eB = np.array(cr['loc']['MB']) - sB; n, m = len(eA), len(eB)
    ovA = (np.arange(n) < pa.min()) | (np.arange(n) > pa.max()); ovB = (np.arange(m) < pb.min()) | (np.arange(m) > pb.max())
    jw_A = np.zeros(n, bool); jw_B = np.zeros(m, bool); jumps = []
    for k in range(len(pa) - 1):
        if (pa[k + 1] - pa[k], pb[k + 1] - pb[k]) != (1, 1):
            jumps.append([int(pa[k]), int(pa[k + 1]), int(pb[k]), int(pb[k + 1])]); jw_A[max(0, pa[k] - w):min(n, pa[k + 1] + w + 1)] = True; jw_B[max(0, pb[k] - w):min(m, pb[k + 1] + w + 1)] = True
    jw_A &= ~ovA; jw_B &= ~ovB; core_A = ~(ovA | jw_A); core_B = ~(ovB | jw_B)
    return dict(n_pibar=int(len(pa)), jumps=jumps, overhang=float(eA[ovA].sum() + eB[ovB].sum()), jump_window=float(eA[jw_A].sum() + eB[jw_B].sum()), core=float(eA[core_A].sum() + eB[core_B].sum()),
                n_overhang_edges=[int(ovA.sum()), int(ovB.sum())], forced_min=0.75 * abs(n - m), total=float(eA.sum() + eB.sum()),
                top_edges_A=[[int(i), round(float(eA[i]), 3)] for i in np.argsort(-eA)[:5]], top_edges_B=[[int(i), round(float(eB[i]), 3)] for i in np.argsort(-eB)[:5]],
                excess_A=[round(float(x), 4) for x in eA], excess_B=[round(float(x), 4) for x in eB])


def cath_H():
    h = {}
    for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
        if l.startswith('#') or not l.strip(): continue
        f = l.split(); h[f[0]] = ('.'.join(f[1:5]), int(f[10]))
    return h


def part_ens():
    TJ = json.load(open(RES + '/r37_trunk_pairs.json'))['r10a_prime']; H = cath_H(); main = {(x['a'], x['b']) for x in TJ['main']}
    for x in TJ['sens']:
        fo = f"{OUTR}/{x['a']}__{x['b']}.json"
        if os.path.exists(fo): continue
        pr = dict(a=x['a'], b=x['b'], n=H[x['a']][1], m=H[x['b']][1], TM=float(x['tm_long']), level='H', Ha=H[x['a']][0], Hb=H[x['b']][0], bin='r10a',
                  ident=float(x['ident']), ln=int(x['ln']), main=(x['a'], x['b']) in main)
        try:
            _CAP.clear(); r = one(pr)
            if 'skip' not in r:
                cr, fa, fb = _CAP[0], _CAP[1], _CAP[2]
                r.update(gap_cross=float(cr['gv'].mean()), gap_self_A=float(fa['gv'].mean()), gap_self_B=float(fb['gv'].mean()))
                try: r['gap_regions'] = gap_regions(cr, fa, fb)
                except Exception as e: r['gap_regions'] = dict(err=repr(e)[:200])
            jdump(r, fo); log('r10a', x['a'], x['b'], r.get('skip', ''), round(r.get('sec', 0)), 's')
        except Exception as e:
            traceback.print_exc(); log('ERR r10a', x['a'], x['b'], repr(e)[:200])
    log('R10A_ENS_DONE')


def part_stats():
    import pandas as pd
    R10 = [json.load(open(f)) for f in sorted(glob.glob(OUTR + '/*.json'))]; R10 = [r for r in R10 if 'skip' not in r]
    B6 = [json.load(open(f)) for f in glob.glob(ROOT + '/derived/r38/b6/*.json')]; hb = [r for r in B6 if 'skip' not in r and r['level'] == 'H']
    hp = pd.read_csv(RES + '/pilot_mdcath_homolog_pairs.tsv', sep='\t'); idn = {}
    for t in hp.itertuples(): idn[(t.a, t.b)] = idn[(t.b, t.a)] = float(t.ident)
    tmc = lambda r: r['TM'] * max(r['n'], r['m']) / max(r['n_pibar'], 1)
    def fit(rows, xs):
        A = np.column_stack([np.ones(len(rows))] + [np.array([f(r) for r in rows]) for f in xs]); y = np.array([r['pair_share_per_edge'] for r in rows]); be = np.linalg.lstsq(A, y, rcond=None)[0]
        res = y - A @ be; return be, float(np.median(np.abs(res - np.median(res))))
    def resid(rows, xs, be):
        A = np.column_stack([np.ones(len(rows))] + [np.array([f(r) for r in rows]) for f in xs]); return np.array([r['pair_share_per_edge'] for r in rows]) - A @ be
    def tost(rows, xs, ref, seed, B=2000):
        be, mad = fit(ref, xs); margin = 0.5 * mad; e = resid(rows, xs, be); Hs = np.array([r['Ha'] for r in rows]); sf = sorted(set(Hs))
        smed = lambda w: np.median([np.median(e[Hs == h]) for h in w]); est = float(smed(sf)); rng = np.random.default_rng(seed)
        bs = [smed([sf[i] for i in rng.integers(0, len(sf), len(sf))]) for _ in range(B)]; lo90, hi90 = np.percentile(bs, [5, 95]); lo95, hi95 = np.percentile(bs, [2.5, 97.5])
        se = float(np.std(bs)); v = ('동등' if (-margin < lo90 and hi90 < margin) else ('다르다' if ((lo95 > 0 or hi95 < 0) and abs(est) > margin) else '가를 수 없다'))
        return dict(n=len(rows), n_sf=len(sf), margin=margin, mad=mad, coef=be.tolist(), est=est, ci90=[float(lo90), float(hi90)], ci95=[float(lo95), float(hi95)], se=se, mde80=2.8 * se, verdict=v)
    main = [r for r in R10 if r.get('main')]; out = dict(written=time.strftime('%F %T'), n_main=len(main), n_sens=len(R10), n_ref=len(hb))
    band = R0_BAND; p1 = [bool(r.get('tmalign') and band[0] <= r['tmalign']['d2n'] <= band[1] and r['D']['d2n'] > band[1]) for r in main]
    out['pred1'] = dict(frac=float(np.mean(p1)) if p1 else None, n=len(p1), band=list(band), note='서술: 관행 열(TM-align 짝 잔기) Δ²/√n̄ 이 R0 띠 안이고 Δ²_D/√n̄ 이 띠 위')
    rows2 = []
    for r in main:
        ex = r['gap_cross'] - (r['gap_self_A'] + r['gap_self_B']) / 2; tgt = 0.75 * r['ln']; gr = r.get('gap_regions') or {}
        rows2.append(dict(ab=[r['a'], r['b']], ln=r['ln'], excess=ex, target=tgt, rel_err=abs(ex - tgt) / tgt, regions={k: gr.get(k) for k in ('overhang', 'jump_window', 'core', 'forced_min', 'n_overhang_edges', 'jumps', 'n_pibar')}))
    out['pred2'] = dict(rows=rows2, median_rel_err=float(np.median([x['rel_err'] for x in rows2])) if rows2 else None, frac_within_10pct=float(np.mean([x['rel_err'] < 0.10 for x in rows2])) if rows2 else None,
                        passed=bool(rows2 and np.median([x['rel_err'] for x in rows2]) < 0.10))
    # (고침 1, 사후·서술; 판정은 위 등록 판) 빈칸 초과를 자리로 나눈 판: 끝 겹침을 뺀 안쪽 초과와 건너뜀 창만의 초과를 0.75 × ln 과 견준다
    ok2 = [x for x in rows2 if x['regions'].get('overhang') is not None]
    if ok2:
        inn = [abs(x['regions']['jump_window'] + x['regions']['core'] - x['target']) / x['target'] for x in ok2]; jw = [abs(x['regions']['jump_window'] - x['target']) / x['target'] for x in ok2]
        out['pred2_regions_posthoc'] = dict(n=len(ok2), frac_excess_in_overhang=float(np.median([x['regions']['overhang'] / x['excess'] for x in ok2 if x['excess'] > 0])) if any(x['excess'] > 0 for x in ok2) else None,
                                            median_rel_err_internal=float(np.median(inn)), median_rel_err_jump_window=float(np.median(jw)),
                                            note='사후 서술(2026-10-20 고침 1): 등록 판의 빈칸 초과에는 끝 겹침 변(길이 차의 대부분)이 들어가 0.75 × ln 과 정의상 맞을 수 없다. 판정은 등록 판 그대로.')
    xtm = [lambda r: r['TM']]; xid = [lambda r: r['TM'], lambda r: r.get('ident', idn.get((r['a'], r['b']), np.nan))]
    hbi = [r for r in hb if (r['a'], r['b']) in idn]
    out['pred3'] = dict(main=tost(main, xtm, hb, 20261057), sens40=tost(R10, xtm, hb, 20261058), ident_cov=tost(main, xid, hbi, 20261059), tm_corrected=tost(main, [tmc], hb, 20261060))
    # (고침 1, 민감도) 기준 회귀선에서 R10a′ 쌍과 겹치는 B6 쌍을 뺀 판(같은 쌍이 기준과 시험에 함께 들어가지 않게)
    tset = {tuple(sorted((r['a'], r['b']))) for r in R10}; hb_lo = [r for r in hb if tuple(sorted((r['a'], r['b']))) not in tset]
    out['pred3']['main_ref_leaveout'] = dict(n_ref_removed=len(hb) - len(hb_lo), **tost(main, xtm, hb_lo, 20261064))
    # (고침 1b, 민감도) 짝 몫을 n̄ 대신 π̄ 짝 수로 나눈 판 — 끝 겹침이 큰 쌍에서 n̄ 이 짝 몫을 묽히는 몫을 뺌(시험·기준 모두)
    rem = lambda rows: [dict(r, pair_share_per_edge=r['pair_share_per_edge'] * r['nbar'] / max(r['n_pibar'], 1)) for r in rows]
    out['pred3']['per_matched_edge'] = tost(rem(main), xtm, rem(hb), 20261065)
    out['dlen'] = dict(test_median=float(np.median([abs(r['nA'] - r['nB']) for r in main])), ref_median=float(np.median([abs(r['nA'] - r['nB']) for r in hb])), note='변 수 차의 중앙값(시험 주 판 대 기준)')
    out['pred1']['note2'] = '상동체 쌍에서는 짝 잔기만의 Δ² 에도 구조 차가 들어가 R0 띠(같은 서열 반복 잡음) 안이기 어렵다. 이 예측의 틀은 같은 단백질 사건(R6)에 맞는다. 값은 서술로 함께 낸다.'
    out['pred1'].update(tmalign_d2n_median=float(np.median([r['tmalign']['d2n'] for r in main if r.get('tmalign')])) if any(r.get('tmalign') for r in main) else None, D_d2n_median=float(np.median([r['D']['d2n'] for r in main])) if main else None)
    out['rule'] = 'prereg 2 절: 판정 기준 = 예측 3 (주 판). TOST 여유 = 0.5 × MAD(상동체 회귀 잔차), 상과 부트스트랩 2,000 의 90% 구간'
    jdump(out, RES + '/r38_r10a.json', bk=True); log('R10A', out['pred1'], out['pred2']['median_rel_err'], out['pred3']['main']['verdict']); log('R10A_STATS_DONE')


if __name__ == '__main__':
    for part in PARTS_R:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
