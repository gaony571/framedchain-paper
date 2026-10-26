"""(prereg round36 1 절, 트랙 B6) 길이가 다른 줄기: 연속 축 281 쌍(results/r37_trunk_pairs.tsv)의 앙상블 점, 그리고 1b 의 281 쌍 r_dyn 판.
앙상블: r38_s0 과 같다 — 320 K 전체 창, 접힌 프레임(Q ≥ 0.75), 반복마다 N = 100 (도메인·온도·창 이름 seed), 실행 쌍마다 프레임 쌍 2,000, 줄인 계산.
쌍마다 낸다:
  Δ²_D, C, B_P, B_Q, w̄, 비, Δ²_D/√n̄ (n̄ = 두 도메인 변 수의 평균); 짝 몫·빈칸 몫(변당) = 교차 프레임 쌍 D² 의 짝·빈칸 분할 평균 ÷ n̄;
  Δ²_{D_matched} (성분) = d_matched = √(짝 몫) 으로 낸 에너지 거리;
  고정 대응 판 둘: π̄(두 평균 점 사이 128 출발 정렬, MeanStart 의 짝)와 관행 열(TM-align 정렬 잔기: 메도이드 Cα, d₀ = d₀(긴 쪽 길이);
    잔기 i↔j 와 i+1↔j+1 이 모두 짝인 변 i↔j). 교차는 그 대응에서 회전만 맞춘 거리, 바닥은 같은 도메인의 같은 변 집합에서 회전만 맞춘 거리 → 에너지 거리.
  J D J: 쪽마다 50 프레임(쪽의 모든 실행을 합친 데서 무작위)을 합친 100 × 100 거리 행렬(줄인 계산) → η₊ = 양의 고유값 합 ÷ 고유값 절댓값 합.
  정확한 W₂ (TM 긴 쪽 ≥ 0.5): 쪽마다 100 프레임 무작위, 지상 비용 D², 정확한 최적 수송(fs.optimal_plan).
  정적 메도이드 D: 메도이드(앙상블 안 d_fix² 합 최소) 사이 BnB(rel_tol 0.005)가 찾은 정렬의 D(유사도 하한 S_lb 로 냄 = D 의 상계; 인증이면 참값과 상대 0.5% 안) → 환산식 상대 잔차 |Δ²_D − (2 D_정적 − B_P − B_Q)| / Δ²_D.
  r_dyn 281 쌍 판: r38_p3 의 도메인 앙상블(320 K 뒤 절반)·평균 점 π̄·공유 변 GPA 윤곽·사영 평균 점 통제·동일성(r38_p3 함수 그대로).
줄인 계산의 매개변수: TM 긴 쪽 < 0.3 쌍은 r38_s0_bias2 의 결정(출발 수·처음 반폭; prereg 0.2 의 다음 단계), 나머지는 r38_s0 그대로(출발 1, 반폭 5). 쌍마다 red_params 로 남긴다.
판정(stats): 예측 1 (정상 점검) Spearman(Δ²_D/√n̄, TM 긴 쪽) < −0.5 이고 상과 짝 군집 부트스트랩 상한 < 0; 2 포화(서술); 3 환산식 잔차 중앙 TM < 0.5 에서 < 15%,
  TM ≥ 0.7 에서 그보다 큼(사후 해석 후보: 퍼짐 합과 잔차의 순위 상관); 4 짝·빈칸 몫과 TM (서술); 5 같은 상과 쌍의 ≥ 80% 에서 η₊ ≤ 0.10; 6 Spearman(W₂, Δ²_D) ≥ 0.8 (TM ≥ 0.5).
  의존 구조: 도메인 재사용 수, 도메인당 한 번 판과 두 도메인 교차 군집 부트스트랩(Exp(1) 도메인 무게의 곱; 예측 1·6 의 민감도). r_dyn 281 판: r38_p3.analyse(통제 + 같은 상과 여부 공변량; 무게 1/상과 짝 쌍 수, 도메인·상과 짝 교차 군집 부트스트랩).
부분: ens (GPU) → stats (CPU).   사용: python r38_b6.py <ROOT> <부분들 쉼표 | all>
출력: derived/r38/b6/<a>__<b>.json, results/r38_b6.json (사본 results/_backup/<날짜>/)"""
import sys, os, json, glob, time, itertools
ROOT = sys.argv[1]; PARTS_B6 = sys.argv[2].split(',') if len(sys.argv) > 2 else ['all']
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)                                      # r38_s0 의 함수 (DEV, side, MeanStart, runpair_means, estats, medoid, ...)
P3NS = {'__name__': 'p3ns'}; _p3 = open(ROOT + '/code/r38_p3.py').read(); _p3 = _p3[:_p3.index("ALL = ['pairs'")]
sys.argv = ['r38_p3', ROOT, 'none', '1']; exec(_p3, P3NS)                              # r38_p3 의 함수는 따로 (이름 충돌 없게)
from stage0.rulers import tm_align, d0_tm
from framespace.transport import optimal_plan
OUTB = ROOT + '/derived/r38/b6'; os.makedirs(OUTB, exist_ok=True)
NJ = 50; NW2 = 100
# prereg 0.2 의 다음 단계(r38_s0_bias2.py)의 결정: TM 긴 쪽 < 0.3 쌍은 그 출발 수·처음 반폭으로 줄인 계산을 한다. 결정 파일이 없으면 (1, K0) = r38_s0 그대로.
_bf = RES + '/r38_s0_bias2.json'; BIAS2 = json.load(open(_bf))['decision'] if os.path.exists(_bf) else dict(n_local=1, k0=K0)
_red_s0 = red_pairs; RED = dict(n_local=1, k0=K0)


def red_pairs(ms, ia, ib, k0=None, iters=10, want_gap=False, membudget=3e9):
    """r38_s0.red_pairs 와 같은 규칙에 출발 수 RED['n_local'](첫째 = 평균 점 출발, 나머지 무작위 seed 0)와 처음 반폭 RED['k0'] 을 건 판. (1, K0) 이면 r38_s0 과 같은 값."""
    nl = RED['n_local']; k0 = RED['k0'] if k0 is None else k0
    EA, EB = ms.EA, ms.EB; n, m = EA.shape[1], EB.shape[1]
    ia = torch.as_tensor(np.asarray(ia), device=EA.device); ib = torch.as_tensor(np.asarray(ib), device=EA.device)
    nA, nB = (EA ** 2).sum((-1, -2)), (EB ** 2).sum((-1, -2))
    S = EA.new_zeros(len(ia)); G = EA.new_zeros(len(ia)) if want_gap else None; kused = np.zeros(len(ia), int); widened = 0
    todo = torch.arange(len(ia), device=EA.device); kk = k0
    while len(todo):
        lo, w = band_windows(ms.center, m, kk); step = max(1, int(membudget / (n * w * 9 * 8 * 4 * nl)))
        edge_all = torch.zeros(len(todo), dtype=torch.bool, device=EA.device)
        for c in range(0, len(todo), step):
            sel = todo[c:c + step]; a, b = EA[ia[sel]], EB[ib[sel]]
            sc, _, M, edge = align_band(a, b, lo, w, ms.starts(ia[sel], ib[sel], n_local=nl), iters=iters)
            S[sel] = sc; edge_all[c:c + step] = edge
            if want_gap: G[sel] = (nA[ia[sel]] * (~M.any(2))).sum(1) + (nB[ib[sel]] * (~M.any(1))).sum(1)
        kused[todo.cpu().numpy()] = kk
        if w >= m: break
        todo = todo[edge_all]
        if kk == k0: widened = int(len(todo))
        kk *= 2
    d2 = (sq_norm(EA)[ia] + sq_norm(EB)[ib] - 2 * S).clamp_min(0)
    return d2.sqrt().cpu().numpy(), (G.cpu().numpy() if want_gap else None), S.cpu().numpy(), kused, widened


def fixed_d(XA, XB, ia, ib, ea, eb):
    """대응 고정(ea ↔ eb) 회전만 맞춘 거리. XA (I, n, 3, 3), XB (J, m, 3, 3) GPU 텐서, 프레임 쌍 (ia, ib)."""
    out = np.empty(len(ia)); ea_t = torch.as_tensor(ea, device=XA.device); eb_t = torch.as_tensor(eb, device=XA.device)
    for c in range(0, len(ia), 20000):
        a = XA[torch.as_tensor(ia[c:c + 20000], device=XA.device)][:, ea_t]; b = XB[torch.as_tensor(ib[c:c + 20000], device=XA.device)][:, eb_t]
        H = torch.einsum('pkcb,pkab->pca', b, a); sv = torch.linalg.svdvals(H); S = sv[:, 0] + sv[:, 1] + torch.sign(torch.det(H)) * sv[:, 2]
        out[c:c + 20000] = (((a ** 2).sum((-1, -2, -3)) + (b ** 2).sum((-1, -2, -3)) - 2 * S).clamp_min(0)).sqrt().cpu().numpy()
    return out


def runpair_fixed(XA, XB, rr, ea, eb, rng, k=NPAIR):
    IA, IB, GR = [], [], []
    for g, (a, b) in enumerate(rr):
        ia, ib = upairs(a, b, rng, k); IA.append(ia); IB.append(ib); GR.append(np.full(len(ia), g))
    IA, IB, GR = map(np.concatenate, (IA, IB, GR)); d = fixed_d(XA, XB, IA, IB, ea, eb)
    return [float(d[GR == g].mean()) for g in range(len(rr))]


def tm_edges(medA, medB):
    """TM-align(재구현, d₀ = 긴 쪽 길이) 의 잔기 짝 → 변 짝(잔기 i↔j 와 i+1↔j+1 이 모두 짝)."""
    PA, PB = np.asarray(fs.ca_trace(medA)), np.asarray(fs.ca_trace(medB)); r = tm_align(PA, PB, d0_tm(max(len(PA), len(PB)))); M = np.asarray(r['match'])
    e = [(i, j) for i, j in zip(*np.nonzero(M)) if i + 1 < len(PA) and j + 1 < len(PB) and M[i + 1, j + 1]]
    return np.array([x[0] for x in e], int), np.array([x[1] for x in e], int), float(r['S'])


def eta_plus(D):
    n = len(D); J = np.eye(n) - np.ones((n, n)) / n; ev = np.linalg.eigvalsh(J @ D @ J); return float(ev[ev > 0].sum() / np.abs(ev).sum())


def p3_pair(a, b):
    """r38_p3 의 함수로 281 쌍 판 r_dyn·통제·동일성."""
    out = {}
    for d in (a, b):
        f = f'{P3NS["OUT"]}/dom/{d}.npz'
        if os.path.exists(f):
            s_, info = P3NS['dom_small'](d); z = P3NS['dom_get'](d) if s_ is not None else None
        else:                                                                   # H-풀 밖 도메인: 같은 함수로 한 번 짓고 캐시에 둔다
            z, info = P3NS['dom_build'](d)
            if z is not None: np.savez(f'{P3NS["CACHE"]}/{d}.npz', **z); s_ = dict(m=z['m'], Pm=z['Pm'], sbar=z['sbar'], degen=z['degen'], **info)
        if z is None: return dict(skip=f'{d} P3 앙상블 없음')
        out[d] = (s_, z['A'])
    (Da, Fa), (Db, Fb) = out[a], out[b]
    al = fs.align(torch.as_tensor(Da['m']).double(), torch.as_tensor(Db['m']).double());          # 기본 출발 = 항등 + 무작위 127 (P3 와 같은 128 출발)
    bi, bj = np.nonzero(al.match.numpy())
    rec = dict(S=float(al.score), n_near=1, n_distinct=1, best=[[int(x), int(y)] for x, y in zip(bi, bj)], inter=[[int(x), int(y)] for x, y in zip(bi, bj)])
    o = P3NS['one_pair'](a, b, rec, Da, Db, Fa, Fb); return o['main']


def one(pr):
    a, b = pr['a'], pr['b']; t0 = time.time(); rng = seeded('r38_b6', a, b); _runs.clear()
    RED['n_local'], RED['k0'] = (int(BIAS2['n_local']), int(BIAS2['k0'])) if pr['TM'] < 0.3 else (1, K0)
    EA, na, La = side(a); EB, nb, Lb = side(b)
    if EA is None or EB is None or na < 2 or nb < 2: return dict(pr, skip='실행 < 2')
    nA, nB = La - 1, Lb - 1; nbar = (nA + nB) / 2
    msAB, msAA, msBB = MeanStart(EA, EB), MeanStart(EA, EA), MeanStart(EB, EB)
    cr = runpair_means(msAB, list(itertools.product(range(na), range(nb))), rng, want_gap=True)
    fa = runpair_means(msAA, list(itertools.combinations(range(na), 2)), rng, want_gap=True)
    fb = runpair_means(msBB, list(itertools.combinations(range(nb), 2)), rng, want_gap=True)
    res = dict(pr, nA=nA, nB=nB, nbar=nbar, runsA=na, runsB=nb, red_params=[RED['n_local'], RED['k0']])
    res['D'] = estats(np.mean(cr['d']), np.mean(fa['d']), np.mean(fb['d']), nbar)
    res['D_matched'] = estats(np.mean(cr['dpair']), np.mean(fa['dpair']), np.mean(fb['dpair']), nbar)
    d2 = cr['dv'] ** 2; res['gap_share_per_edge'] = float(cr['gv'].mean() / nbar); res['pair_share_per_edge'] = float((d2 - cr['gv']).mean() / nbar)
    res['widened'] = dict(cross=cr['widened'], A=fa['widened'], B=fb['widened'])
    # 고정 대응: π̄ (평균 점 정렬) 와 TM-align
    pi = msAB.match.cpu().numpy(); pi_a, pi_b = np.nonzero(pi); res['n_pibar'] = int(len(pi_a))
    medA, medB = medoid(a), medoid(b); ta, tb, Stm = tm_edges(medA.cpu().numpy(), medB.cpu().numpy()); res['n_tm_edges'] = int(len(ta)); res['S_tm'] = Stm
    for tag, ea, eb in (('pibar', pi_a, pi_b), ('tmalign', ta, tb)):
        if len(ea) < 10: res[tag] = None; continue
        C_ = runpair_fixed(EA, EB, list(itertools.product(range(na), range(nb))), ea, eb, rng)
        BP = runpair_fixed(EA, EA, list(itertools.combinations(range(na), 2)), ea, ea, rng); BQ = runpair_fixed(EB, EB, list(itertools.combinations(range(nb), 2)), eb, eb, rng)
        res[tag] = estats(np.mean(C_), np.mean(BP), np.mean(BQ), nbar)
    # J D J
    ja = np.sort(rng.choice(len(EA), min(NJ, len(EA)), replace=False)); jb = np.sort(rng.choice(len(EB), min(NJ, len(EB)), replace=False))
    iu = np.triu_indices(len(ja), 1); daa = red_pairs(msAA, ja[iu[0]], ja[iu[1]])[0]; iv = np.triu_indices(len(jb), 1); dbb = red_pairs(msBB, jb[iv[0]], jb[iv[1]])[0]
    ga, gb = np.meshgrid(ja, jb, indexing='ij'); dab = red_pairs(msAB, ga.ravel(), gb.ravel())[0].reshape(len(ja), len(jb))
    Daa = np.zeros((len(ja), len(ja))); Daa[iu] = daa; Daa = Daa + Daa.T; Dbb = np.zeros((len(jb), len(jb))); Dbb[iv] = dbb; Dbb = Dbb + Dbb.T
    res['eta_plus'] = eta_plus(np.block([[Daa, dab], [dab.T, Dbb]]))
    # 정확한 W₂ (TM ≥ 0.5)
    if pr['TM'] >= 0.5:
        wa_ = np.sort(rng.choice(len(EA), min(NW2, len(EA)), replace=False)); wb_ = np.sort(rng.choice(len(EB), min(NW2, len(EB)), replace=False))
        ga, gb = np.meshgrid(wa_, wb_, indexing='ij'); dw = red_pairs(msAB, ga.ravel(), gb.ravel())[0].reshape(len(wa_), len(wb_))
        Cw = dw ** 2; Pl = optimal_plan(Cw); res['W2'] = float(np.sqrt(max((Pl * Cw).sum(), 0)))
    # 정적 메도이드 D (BnB)
    cert = fs.bnb(medA, medB, rel_tol=REL, device=DEV); nx, ny = float(sq_norm(medA)), float(sq_norm(medB))
    res['static'] = dict(D=float(np.sqrt(max(nx + ny - 2 * cert['S_lb'], 0))), certified=bool(cert['certified']), sec=cert['sec'])
    res['conv_rel_resid'] = float(abs(res['D']['D2'] - (2 * res['static']['D'] - res['D']['B_P'] - res['D']['B_Q'])) / abs(res['D']['D2'])) if res['D']['D2'] != 0 else None
    # r_dyn 281 쌍 판
    try: res['p3'] = p3_pair(a, b)
    except Exception as e: res['p3'] = dict(err=repr(e)[:150])
    res['sec'] = time.time() - t0
    return res


def part_ens():
    import pandas as pd
    T = pd.read_csv(RES + '/r37_trunk_pairs.tsv', sep='\t'); prs = [dict(a=r.a, b=r.b, n=int(r.n), m=int(r.m), TM=float(r.TM), level=str(r.level), Ha=str(r.Ha), Hb=str(r.Hb), bin=str(r.bin)) for r in T.itertuples()]
    log('B6 쌍', len(prs))
    for pr in prs:
        fo = f"{OUTB}/{pr['a']}__{pr['b']}.json"
        if os.path.exists(fo): continue
        try:
            r = one(pr); jdump(r, fo); log('b6', pr['bin'], pr['a'], pr['b'], r.get('skip', ''), round(r.get('sec', 0)), 's')
        except Exception as e:
            traceback.print_exc(); log('ERR b6', pr['a'], pr['b'], repr(e)[:200])
    log('B6_ENS_DONE')


def part_stats():
    from scipy.stats import spearmanr
    J = [json.load(open(f)) for f in sorted(glob.glob(OUTB + '/*.json'))]; J = [j for j in J if 'skip' not in j and j.get('D')]
    tm = np.array([j['TM'] for j in J]); d2n = np.array([j['D']['d2n'] for j in J]); D2 = np.array([j['D']['D2'] for j in J])
    cl = np.array(['|'.join(sorted([j['Ha'], j['Hb']])) for j in J]); same = np.array([j['Ha'] == j['Hb'] for j in J]); rng = np.random.default_rng(20261032)
    def clboot(x, y, B=2000):
        u = sorted(set(cl)); idx = {c: np.where(cl == c)[0] for c in u}; bs = []
        for _ in range(B):
            pick = np.concatenate([idx[u[i]] for i in rng.integers(0, len(u), len(u))]); bs.append(spearmanr(x[pick], y[pick]).correlation)
        return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
    doms = sorted({j['a'] for j in J} | {j['b'] for j in J}); di = {d: i for i, d in enumerate(doms)}; iA = np.array([di[j['a']] for j in J]); iB = np.array([di[j['b']] for j in J])
    def xboot(x, y, sel=None, B=2000, seed=20261035):
        # 두 도메인 교차 군집 부트스트랩: 도메인마다 Exp(1) 무게, 쌍 무게 = 두 도메인 무게의 곱, 무게 Spearman(순위의 무게 Pearson)
        sel = np.arange(len(J)) if sel is None else np.asarray(sel); x, y = np.asarray(x, float), np.asarray(y, float); r_ = np.random.default_rng(seed); rx, ry = P3NS['ranks'](x), P3NS['ranks'](y); bs = []
        for _ in range(B):
            w_ = r_.exponential(1, len(doms)); bs.append(P3NS['wcorr'](rx, ry, w_[iA[sel]] * w_[iB[sel]]))
        return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
    reuse = {}
    for j in J:
        for d in (j['a'], j['b']): reuse[d] = reuse.get(d, 0) + 1
    order = np.random.default_rng(20261033).permutation(len(J)); used = set(); once = []
    for k in order:
        if J[k]['a'] not in used and J[k]['b'] not in used: once.append(k); used |= {J[k]['a'], J[k]['b']}
    once = np.array(sorted(once))
    p1 = dict(spearman=float(spearmanr(d2n, tm).correlation), ci=clboot(d2n, tm), once_domain=dict(n=len(once), spearman=float(spearmanr(d2n[once], tm[once]).correlation)), ci_domain_xcluster=xboot(d2n, tm))
    p1['passed'] = bool(p1['spearman'] < -0.5 and p1['ci'][1] < 0)
    bins = sorted(set(j['bin'] for j in J)); p2 = {b: dict(n=int(sum(j['bin'] == b for j in J)), d2n_median=float(np.median([j['D']['d2n'] for j in J if j['bin'] == b]))) for b in bins}
    cr = np.array([j['conv_rel_resid'] if j['conv_rel_resid'] is not None else np.nan for j in J]); lo_, hi_ = tm < 0.5, tm >= 0.7
    spread = np.array([j['D']['B_P'] + j['D']['B_Q'] for j in J]); okc = np.isfinite(cr)
    p3 = dict(median_tm_lt05=float(np.nanmedian(cr[lo_])), median_tm_ge07=float(np.nanmedian(cr[hi_])), n_lt05=int(lo_.sum()), n_ge07=int(hi_.sum()),
              posthoc_spread_vs_resid=float(spearmanr(spread[okc], cr[okc]).correlation), uncertified=int(sum(not j['static']['certified'] for j in J)))
    p3['passed'] = bool(p3['median_tm_lt05'] < 0.15 and p3['median_tm_ge07'] > p3['median_tm_lt05'])
    gs = np.array([j['gap_share_per_edge'] for j in J]); ps = np.array([j['pair_share_per_edge'] for j in J])
    p4 = dict(gap_vs_tm=float(spearmanr(gs, tm).correlation), pair_vs_tm=float(spearmanr(ps, tm).correlation))
    eta = np.array([j['eta_plus'] for j in J]); p5 = dict(frac_same_sf_le_010=float(np.mean(eta[same] <= 0.10)), n_same=int(same.sum()), frac_other_le_010=float(np.mean(eta[~same] <= 0.10)) if (~same).any() else None,
                                                          flagged_same=[f"{J[k]['a']}__{J[k]['b']}" for k in np.where(same & (eta > 0.10))[0]])
    p5['passed'] = bool(p5['frac_same_sf_le_010'] >= 0.80)
    w = [(j['W2'], j['D']['D2']) for j in J if 'W2' in j]; p6 = dict(n=len(w), spearman=float(spearmanr([x[0] for x in w], [x[1] for x in w]).correlation) if len(w) > 2 else None)
    kw = np.array([k for k, j in enumerate(J) if 'W2' in j], int); ow = [k for k in once if 'W2' in J[k]]
    if len(kw) > 2: p6.update(ci_domain_xcluster=xboot([J[k]['W2'] for k in kw], [J[k]['D']['D2'] for k in kw], sel=kw), once_domain=dict(n=len(ow), spearman=float(spearmanr([J[k]['W2'] for k in ow], [J[k]['D']['D2'] for k in ow]).correlation) if len(ow) > 2 else None))
    p6['passed'] = bool(p6['spearman'] is not None and p6['spearman'] >= 0.8)
    cols = {k: dict(spearman_vs_D=float(spearmanr([j[k]['D2'] for j in J if j.get(k)], [j['D']['D2'] for j in J if j.get(k)]).correlation),
                    median_ratio_vs_D=float(np.median([j[k]['D2'] / j['D']['D2'] for j in J if j.get(k) and j['D']['D2'] != 0]))) for k in ('D_matched', 'pibar', 'tmalign')}
    rows = [dict(ab=[j['a'], j['b']], H='|'.join(sorted([j['Ha'], j['Hb']])), main=dict(j['p3'], same_sf=float(j['Ha'] == j['Hb']))) for j in J if isinstance(j.get('p3'), dict) and j['p3'].get('r_dyn') is not None]
    rdyn281 = P3NS['analyse'](rows, z=('control', 'same_sf'), seed=20261034) if rows else None
    out = dict(written=time.strftime('%F %T'), n=len(J), n_skipped=len(glob.glob(OUTB + '/*.json')) - len(J), domain_reuse=dict(max=max(reuse.values()), n_reused=sum(v > 1 for v in reuse.values())),
               pred1_sanity=p1, pred2_saturation=p2, pred3_conversion=p3, pred4_shares=p4, pred5_eta=p5, pred6_w2=p6, columns=cols, r_dyn_281=rdyn281,
               verdict=dict(p3='맞음' if p3['passed'] else '틀림', p5='맞음' if p5['passed'] else '틀림', p6='맞음' if p6['passed'] else '틀림', p1_sanity='섬' if p1['passed'] else '서지 않음 — 계산부터 의심'),
               rule='prereg 1 절: 판정 기준 = 예측 3·5·6 (2차), 1 은 정상 점검, 2·4 는 서술')
    jdump(out, RES + '/r38_b6.json', bk=True); log('B6', out['verdict']); log('B6_STATS_DONE')


ALL_B6 = ['ens', 'stats']
if __name__ == '__main__':
    for part in (ALL_B6 if PARTS_B6 == ['all'] else PARTS_B6):
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('B6_DONE', ','.join(PARTS_B6))
