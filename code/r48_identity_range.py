"""(등록 = Drive results/prereg_2026-10-30_addendum11.json "r48_identity_range") 확인 실험 C 의 연관이 동일성 전 범위에 걸치는가, 동일성이 높은 소수의 쌍이 끄는가.
자료: 등록 판 C 의 쌍별 값(r38_c.merged 의 main: r_dyn, π̄ 동일성, control, r_anm), r46 의 쌍별 TM-align 값(derived/r46/align_pairs.jsonl: ident_tm_res, prec0, n_res),
      그리고 이 스크립트가 도메인 파일(derived/r38/c/dom/<d>.npz 의 서열·퇴화 표시)과 π̄ 짝(pibar.jsonl)에서 새로 계산하는 조성 기대 동일성 Σ_a p_a q_a (공유 결합의 첫 잔기 조성).
순위 = 평균 순위(scipy rankdata average), 부분집합마다 그 안에서 r_dyn·동일성·통제의 순위를 다시 매김. 무게 = 1/(그 집합 안의 상과 쌍 수). 편상관 식은 r38_p3 그대로(r38_c 이름 공간).
판정 집합 S15 = π̄ 동일성 < 0.15 (문턱은 원고 표 3 의 기존 동일성 묶음 < 15% 에서 가져옴). 판정 = C-1·C-2 하한 > 0 → holds, C-1 만 → partly, 아니면 fails.
서술: 전체 판; S15 의 상과만으로 낸 전체 판(SFbase); S10 = π̄ 동일성 < 0.10; T15 = TM 정렬 잔기 동일성 < 0.15 로 고르고 x 는 π̄ 동일성;
      S15 안의 A (x = TM 정렬 잔기 동일성) 와 B0 (통제 + 정밀도 ±0) — 쌍은 S15 ∩ TM 정렬 잔기 ≥ 30 (r46 규칙);
      전체 − S15 짝지은 차(서술만: 범위를 좁히면 연관이 고르게 퍼져 있어도 양수로 나옴);
      범위 제한 기대치 = 전체 자료에 순위 척도에서 같은 연관을 둔 모형(상과 무게 WLS: r_dyn 순위 ~ 동일성 순위 + 통제 순위)의 적합값에,
        상과 × 동일성 층(≥ 0.15 / < 0.15) 안에서 복원 추출한 잔차를 더한 가상 r_dyn 으로 S15 추정을 1,000 번 반복한 평균; 두 층의 잔차 분산;
      비교 규칙: S15 상한 < 기대 평균 → below, S15 하한 > 기대 평균 → above, 그 밖 → consistent;
      동일성 십분위별 r_dyn 잔차(C-2 통제 뒤)의 가중 평균과 구간; 관측·조성 기대 동일성의 중앙값.
구간 = 상과 Exp(1) 군집 부트스트랩 10,000: seed 20261131 한 흐름의 상과 무게(242 상과)를 모든 판이 함께 씀(짝지은 차가 짝지어짐). 모사 seed 20261132 (SeedSequence 로 1,000 갈래).
점검: (1) 전체 판의 평균 순위 점추정 = r46_align_fix1 고침 기록 값(C-1 0.22460683061158487, C-2 0.18094574136027425), |차| < 1e-9.
      (2) 새로 센 관측 동일성 = 저장 ident, 모든 쌍에서 |차| < 1e-12. 어긋나면 판정을 내지 않음.
부분: exp → stats   사용: python r48_identity_range.py <ROOT> exp,stats   출력: derived/r48/exp_ident.jsonl, results/r48_identity_range.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
import sys, json, time, traceback
import numpy as np
ROOT = sys.argv[1]; PARTS = sys.argv[2].split(',')
RES, DER = ROOT + '/results', ROOT + '/derived'; CD = DER + '/r38/c'; OUT48 = DER + '/r48'; os.makedirs(OUT48, exist_ok=True)
FE = OUT48 + '/exp_ident.jsonl'; FO = RES + '/r48_identity_range.json'
SEED, SEEDSIM, NBOOT, NSIM, NW, MIN_RES, CUT, CUT2 = 20261131, 20261132, 10000, 1000, 12, 30, 0.15, 0.10
REC_AVG = dict(C1=0.22460683061158487, C2=0.18094574136027425)
Z1, Z2 = ('control',), ('control', 'r_anm')


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)
def rk(v):
    from scipy.stats import rankdata
    return rankdata(v, method='average').astype(float)
def jl(f): return [json.loads(l) for l in open(f) if l.strip()]


def load_ns(path, cut, argv):
    src = open(path).read(); src = src[:src.index(cut)]
    g = {'__name__': 'r48_' + os.path.basename(path)}; old = sys.argv; sys.argv = argv
    try: exec(src, g)
    finally: sys.argv = old
    return g


G = load_ns(ROOT + '/code/r38_c.py', "ALLC = ['pairs'", ['r38_c', ROOT, 'none', '1']); G['ranks'] = rk
pcorr, wls_resid = G['pcorr'], G['wls_resid']


def part_exp():
    from concurrent.futures import ThreadPoolExecutor
    P = json.load(open(CD + '/pairs.json'))['pairs']; PB = {tuple(r['ab']): r for r in jl(CD + '/pibar.jsonl')}; PO = {tuple(r['ab']): r for r in jl(CD + '/pairs_out.jsonl')}
    doms = sorted({p[0] for p in P} | {p[1] for p in P})
    def ld(d):
        f = f'{CD}/dom/{d}.npz'
        if not os.path.exists(f): return d, None
        z = np.load(f); info = json.loads(str(z['info']))
        return d, (info.get('seq'), z['degen'] if 'degen' in z.files else None)
    with ThreadPoolExecutor(16) as ex: D = dict(ex.map(ld, doms))
    n = 0
    with open(FE, 'x') as fh:
        for p in P:
            k = (p[0], p[1]); rec = PB.get(k); po = (PO.get(k) or {}).get('main') or {}
            if rec is None or po.get('ident') is None or D.get(p[0]) is None or D.get(p[1]) is None: continue
            (sa, ga), (sb, gb) = D[p[0]], D[p[1]]
            if sa is None or ga is None or gb is None: continue
            prs = np.array(rec['best'], int).reshape(-1, 2); keep = ~(ga[prs[:, 0]] | gb[prs[:, 1]]) if len(prs) else np.zeros(0, bool); sh = prs[keep]
            if len(sh) == 0: continue
            ra = [sa[i] for i in sh[:, 0]]; rb = [sb[j] for j in sh[:, 1]]
            aa = sorted(set(ra) | set(rb)); pa = {a: ra.count(a) / len(ra) for a in aa}; pb = {a: rb.count(a) / len(rb) for a in aa}
            fh.write(json.dumps(dict(ab=list(k), H=p[2], ident_obs=float(np.mean([x == y for x, y in zip(ra, rb)])), ident_stored=po['ident'],
                                     ident_exp=float(sum(pa[a] * pb[a] for a in aa)), n_shared=int(len(sh)))) + '\n'); n += 1
    log('exp 끝', n)


SP = {}   # 부트스트랩 일꾼이 쓰는 판들(fork 로 물려받음)


def _boot_chunk(ix):
    V = SP['V']; out = {nm: np.empty(len(ix)) for nm in SP['specs']}
    for t, b in enumerate(ix):
        for nm, (Y, X, Z, w0, si) in SP['specs'].items(): out[nm][t] = pcorr(Y, X, Z, w0 * V[b][si])
    return ix, out


def _sim_one(args):
    b, key = args; S = SP['sim'][key]; rng = np.random.default_rng(SP['ss'][b])
    e = S['e']; es = np.empty_like(e)
    for grp in S['groups']: es[grp] = e[rng.choice(grp, len(grp), replace=True)]
    Ys = S['F'] + es; m = S['mask']
    return key, b, pcorr(rk(Ys[m]), S['Xs'], S['Zs'], S['ws'])


def part_stats():
    from multiprocessing import Pool
    base = G['merged'](CD, CD + '/anm_pairs.jsonl')
    c1 = G['analyse_sf'](base, 'main', z=Z1, seed=20261071); c2 = G['analyse_sf'](base, 'main', z=Z2, seed=20261074)
    dev = max(abs(c1['est'] - REC_AVG['C1']), abs(c2['est'] - REC_AVG['C2']))
    E = {tuple(r['ab']): r for r in jl(FE)}; AL = {tuple(r['ab']): r for r in jl(DER + '/r46/align_pairs.jsonl')}
    dev_id = max(abs(r['ident_obs'] - r['ident_stored']) for r in E.values())
    out = dict(written=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), check=dict(avg_dev=float(dev), avg_ok=bool(dev < 1e-9), ident_dev=float(dev_id), ident_ok=bool(dev_id < 1e-12), n_exp=len(E)))
    if not (dev < 1e-9 and dev_id < 1e-12):
        out['stopped'] = '점검 실패 — 판정을 내지 않음'; json.dump(out, open(FO, 'x'), ensure_ascii=False, indent=1); log('R48 STOP', out['check']); return
    rows = [r for r in base if all(r['main'].get(c) is not None for c in ('r_dyn', 'ident', 'control', 'r_anm'))]
    H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); sig = np.array([sfl.index(h) for h in H])
    idn = np.array([r['main']['ident'] for r in rows]); tmr = [AL.get(tuple(r['ab'])) for r in rows]
    tm_ok = np.array([a is not None and a['n_res'] >= MIN_RES and a.get('prec0') is not None for a in tmr])
    tm_id = np.array([a['ident_tm_res'] if a is not None and a['n_res'] >= MIN_RES else np.nan for a in tmr])
    s15 = idn < CUT; sfs = set(H[s15])
    sets = dict(full=np.ones(len(rows), bool), S15=s15, SFbase=np.array([h in sfs for h in H]), S10=idn < CUT2, T15=np.nan_to_num(tm_id, nan=9) < CUT, S15A=s15 & tm_ok)
    def arr(mask, xkey, z):
        ix = np.where(mask)[0]; Hs = sig[ix]; cnt = np.bincount(Hs, minlength=len(sfl)); w0 = 1.0 / cnt[Hs]
        get = lambda c: np.array([tmr[i]['prec0'] if c == 'prec0' else (tmr[i]['ident_tm_res'] if c == 'ident_tm' else rows[i]['main'][c]) for i in ix], float)
        return (rk(get('r_dyn')), rk(get(xkey)), [rk(get(c)) for c in z], w0, Hs)
    specs = {}
    for nm in ('full', 'S15', 'SFbase', 'S10', 'T15'):
        specs[f'{nm}_C1'] = arr(sets[nm], 'ident', Z1); specs[f'{nm}_C2'] = arr(sets[nm], 'ident', Z2)
    specs['A15_C1'] = arr(sets['S15A'], 'ident_tm', Z1); specs['A15_C2'] = arr(sets['S15A'], 'ident_tm', Z2)
    specs['B015_C1'] = arr(sets['S15A'], 'ident', Z1 + ('prec0',)); specs['B015_C2'] = arr(sets['S15A'], 'ident', Z2 + ('prec0',))
    est = {nm: float(pcorr(Y, X, Z, w0)) for nm, (Y, X, Z, w0, si) in specs.items()}
    V = np.random.default_rng(SEED).exponential(1, (NBOOT, len(sfl)))   # seed 한 흐름: 반복마다 상과 242 개의 무게
    SP.update(V=V, specs=specs); bs = {nm: np.empty(NBOOT) for nm in specs}
    with Pool(NW) as pool:
        for ix, o in pool.imap_unordered(_boot_chunk, np.array_split(np.arange(NBOOT), NW * 4)):
            for nm in o: bs[nm][ix] = o[nm]
    ci = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    res = {nm: dict(n=int(len(specs[nm][0])), n_sf=int(len(set(specs[nm][4]))), est=est[nm], lo=ci(bs[nm])[0], hi=ci(bs[nm])[1]) for nm in specs}
    diff = {c: dict(est=est[f'full_{c}'] - est[f'S15_{c}'], lo=ci(bs[f'full_{c}'] - bs[f'S15_{c}'])[0], hi=ci(bs[f'full_{c}'] - bs[f'S15_{c}'])[1]) for c in ('C1', 'C2')}
    log('부트스트랩 끝', {k: round(v, 4) for k, v in est.items()})
    # 범위 제한 기대치 (순위 척도에서 같은 연관)
    sim = {}; strata = (idn >= CUT).astype(int); resvar = {}
    for c, z in (('C1', Z1), ('C2', Z2)):
        Y, X, Zf, w0, _ = specs[f'full_{c}']; Xd = np.column_stack([np.ones(len(Y)), X] + Zf); W = np.sqrt(w0)
        beta = np.linalg.lstsq(Xd * W[:, None], Y * W, rcond=None)[0]; F = Xd @ beta; e = Y - F
        groups = [np.where((sig == s) & (strata == t))[0] for s in range(len(sfl)) for t in (0, 1)]; groups = [g_ for g_ in groups if len(g_)]
        resvar[c] = {('ge15' if t else 'lt15'): float(np.average(e[strata == t] ** 2, weights=w0[strata == t])) for t in (0, 1)}
        Ys_, Xs_, Zs_, ws_, _ = specs[f'S15_{c}']
        sim[c] = dict(F=F, e=e, groups=groups, mask=s15, Xs=Xs_, Zs=Zs_, ws=ws_, beta_x=float(beta[1]))
    SP.update(sim=sim, ss=np.random.SeedSequence(SEEDSIM).spawn(NSIM)); vals = {c: np.empty(NSIM) for c in sim}
    with Pool(NW) as pool:
        for c, b, v in pool.imap_unordered(_sim_one, [(b, c) for c in sim for b in range(NSIM)], chunksize=20): vals[c][b] = v
    rr = {}
    for c in ('C1', 'C2'):
        m = float(vals[c].mean()); s = res[f'S15_{c}']
        rr[c] = dict(expected_mean=m, sim_sd=float(vals[c].std()), beta_x_rank=sim[c]['beta_x'], reading='below' if s['hi'] < m else 'above' if s['lo'] > m else 'consistent')
    # 동일성 십분위별 r_dyn 잔차 (C-2 통제 뒤)
    Y, X, Zf, w0, si = specs['full_C2']; rY = wls_resid(Y, Zf, w0); qs = np.quantile(idn, np.linspace(0, 1, 11)); bins = np.clip(np.searchsorted(qs, idn, side='right') - 1, 0, 9)
    dec = []
    for k in range(10):
        m_ = bins == k; vb = [float((w0[m_] * V[b][si[m_]] * rY[m_]).sum() / (w0[m_] * V[b][si[m_]]).sum()) for b in range(2000)]
        dec.append(dict(ident_lo=float(qs[k]), ident_hi=float(qs[k + 1]), n=int(m_.sum()), mean_rY=float((w0[m_] * rY[m_]).sum() / w0[m_].sum()), ci95=ci(vb)))
    ex = lambda mask: dict(n=int(mask.sum()), obs_median=float(np.median([E[tuple(rows[i]['ab'])]['ident_obs'] for i in np.where(mask)[0] if tuple(rows[i]['ab']) in E])),
                           exp_median=float(np.median([E[tuple(rows[i]['ab'])]['ident_exp'] for i in np.where(mask)[0] if tuple(rows[i]['ab']) in E])),
                           diff_median=float(np.median([E[tuple(rows[i]['ab'])]['ident_obs'] - E[tuple(rows[i]['ab'])]['ident_exp'] for i in np.where(mask)[0] if tuple(rows[i]['ab']) in E])))
    def verdict(a, b):
        return 'holds' if a['lo'] > 0 and b['lo'] > 0 else 'partly' if a['lo'] > 0 else 'fails'
    vS = verdict(res['S15_C1'], res['S15_C2'])
    A_ok = bool(res['A15_C1']['lo'] > 0 and res['A15_C2']['lo'] > 0); B0_ok = bool(res['B015_C1']['lo'] > 0 and res['B015_C2']['lo'] > 0)
    case = (('holds_AB' if A_ok and B0_ok else 'holds_Afail' if not A_ok else 'holds_Bfail') if vS == 'holds' else vS)
    out.update(sets={k: int(v.sum()) for k, v in sets.items()}, n_sf=dict(full=len(sfl), S15=len(sfs)), est=res, diff_full_minus_S15=diff, range_restriction=rr, residual_var=resvar,
               deciles_C2=dec, composition=dict(full=ex(sets['full']), S15=ex(s15)),
               verdict=dict(S15=vS, A15_ok=A_ok, B015_ok=B0_ok, case=case),
               rule='등록 r48_identity_range: S15 의 C-1·C-2 하한 > 0 → holds, C-1 만 → partly, 아니면 fails. A15·B015 는 서술(각각 C-1·C-2 하한 둘 다 > 0 이면 ok)이고 문장 갈래만 가름. '
                    '범위 제한 기대치 비교: S15 상한 < 기대 → below, 하한 > 기대 → above, 그 밖 consistent. 전체 − S15 차는 서술만.')
    json.dump(out, open(FO, 'x'), ensure_ascii=False, indent=1); log('R48', out['verdict'], rr); log('R48_DONE')


if __name__ == '__main__':
    for part in PARTS:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR', part, repr(e)[:300])
