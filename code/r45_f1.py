"""(round45, 등록 = Drive results/prereg_2026-10-29_addendum7.json "r45_f1", 2026-10-29 UTC) 확인 실험 C 의 C-2 가 Cα 모델이 못 보는 펩타이드 면 요동에서 오는가.
r_dyn 의 MD 프로파일을 세 판으로 낸다. 쌍·공유 결합·앙상블·GPA 겹침·통제(control)·r_ANM·동일성은 등록 판(r38_c)과 같다.
  full  : 등록 판 그대로. GPA 뒤 결합마다 ‖x_i − 평균‖² 의 프레임 평균 (r38_p3.gpa_profile 과 같은 식).
  tilt  : 같은 GPA 뒤 단위 벡터 f₁ = x_i 첫째 열 / 그 노름 의 평균 제곱 편차 (결합 방향의 기울기; ANM 예측과 같은 각도 분산).
  twist : 단위 f₂·f₃ 의 평균 제곱 편차 합에서 tilt 를 빼고 2 로 나눔 = 결합 축을 도는 비틀림 (작은 회전 δf_k = ω × f_k 에서 정확, 큰 회전에서는 근사).
부분: corr K k → stats, sanity
  corr  : C 쌍(r44_corr_v2.c_pairs 와 같은 쌍)마다 π̄ 공유 결합(퇴화 뺌)에서 세 판의 r_dyn. 점검: full 판 r_dyn 이 저장된 pairs_out.jsonl 의 main 과 같은지(모든 쌍).
  stats : r38_c.merged 행(저장된 control·ident·r_ANM)에 세 판 r_dyn 을 붙여, 같은 쌍 집합·같은 상과 Exp(1) 무게(analyse_sf 와 같은 난수열)로
          C-1 (seed 20261071) · C-2 (seed 20261074) 를 세 판과 짝 차이(full − tilt, full − twist)로 한 번에 낸다(2,000). 경계면 둘째 seed 10,000 (서술).
  sanity: C 도메인마다 ANM 예측과 MD 프로파일(전체 결합)의 Spearman, 세 판(서술). full 판이 저장값 rho 와 같은지 점검.
사용: python r45_f1.py <ROOT> corr <K> <k> | stats | sanity   출력: derived/r45/f1/pairs_s<k>.jsonl, derived/r45/f1/sanity.jsonl, results/r45_f1.json (사본 results/_backup/<날짜>/)"""
import os, sys
ROOT = sys.argv[1]; PARTS45 = sys.argv[2].split(','); K45 = int(sys.argv[3]) if len(sys.argv) > 3 else 1; S45 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
_src = open(ROOT + '/code/r38_c.py').read(); _src = _src[:_src.index("ALLC = ['pairs'")]
sys.argv = ['r38_c', ROOT, 'none', '1']; exec(_src)   # r38_c (그 안에서 r38_p3) 의 함수: dom_get, dom_small, gpa_profile, spearman, ranks, pcorr, analyse_sf, merged, MINSH
OUT45 = ROOT + '/derived/r45/f1'; os.makedirs(OUT45, exist_ok=True)
VERS = ('full', 'tilt', 'twist'); SEED_C1, SEED_C2, SEED_2ND = 20261071, 20261074, 20261113; EDGE = 0.01


def profiles(F, idx, iters=3):
    """gpa_profile 와 같은 GPA(자기 평균 둘레) 뒤 세 판의 결합별 퍼짐 (len(idx),) 셋."""
    X = F[:, idx].astype(np.float64); mu = X.mean(0)
    for _ in range(iters):
        A = np.einsum('kiab,icb->kac', X, mu); W, S, Zt = np.linalg.svd(A); dt = np.sign(np.linalg.det(np.einsum('kba,kcb->kac', Zt, W)))
        Dg = np.ones((len(X), 3)); Dg[:, 2] = dt; Uk = np.einsum('kba,kb,kcb->kac', Zt, Dg, W)
        X = np.einsum('kac,kicb->kiab', Uk, X); mu = X.mean(0)
    full = ((X - mu) ** 2).sum((-1, -2)).mean(0)
    Fu = X / np.linalg.norm(X, axis=-2, keepdims=True)                 # 열마다 단위 벡터 f_k (k 는 마지막 축)
    msd = ((Fu - Fu.mean(0)) ** 2).sum(-2).mean(0)                     # (len(idx), 3): 열 f₁ f₂ f₃ 의 평균 제곱 편차
    return dict(full=full, tilt=msd[:, 0], twist=(msd[:, 1] + msd[:, 2] - msd[:, 0]) / 2)


def c_pairs():
    P = json.load(open(OUT + '/pairs.json'))['pairs']; PB = {}
    for l in open(OUT + '/pibar.jsonl'):
        r = json.loads(l); PB[tuple(r['ab'])] = r
    P = [p for p in P if (p[0], p[1]) in PB and dom_small(p[0])[0] is not None and dom_small(p[1])[0] is not None]
    return P, PB


# ---------- corr ----------
def part_corr():
    P, PB = c_pairs(); fo = OUT45 + f'/pairs_s{S45}.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    stored = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(OUT + '/pairs_out.jsonl')}; small, big = {}, {}
    def get(d):
        if d not in small: small[d] = dom_small(d)[0]
        if d not in big:
            if len(big) > 40: big.pop(next(iter(big)))
            z = dom_get(d); big[d] = None if z is None else z['A']
        return small[d], big[d]
    allp = sorted(P, key=lambda p: (p[0], p[1])); blk = np.array_split(np.arange(len(allp)), K45)[S45]
    todo = [allp[i] for i in blk if (allp[i][0], allp[i][1]) not in done]; log('r45 corr 조각', S45, '/', K45, len(todo), '쌍 남음'); t0 = time.time()
    with open(fo, 'a') as fh:
        for k, p in enumerate(todo):
            try:
                (Da, Fa), (Db, Fb) = get(p[0]), get(p[1]); o = dict(ab=[p[0], p[1]], H=p[2])
                prs = np.array(PB[(p[0], p[1])]['best'], int).reshape(-1, 2); keep = ~(Da['degen'][prs[:, 0]] | Db['degen'][prs[:, 1]]) if len(prs) else np.zeros(0, bool); sh = prs[keep]
                o['n_shared'] = int(len(sh))
                if len(sh) >= MINSH and Fa is not None and Fb is not None:
                    pa, pb = profiles(Fa, sh[:, 0]), profiles(Fb, sh[:, 1])
                    for v in VERS: o['r_' + v] = spearman(pa[v], pb[v])
                    o['twist_neg_frac'] = float(np.mean(np.concatenate([pa['twist'], pb['twist']]) < 0))
                s = stored.get((p[0], p[1]), {}).get('main', {})
                if s.get('r_dyn') is not None and o.get('r_full') is not None: o['check'] = o['r_full'] - s['r_dyn']
                fh.write(json.dumps(o) + '\n')
            except Exception as e:
                log('ERR r45 corr', p[:2], repr(e)[:150])
            if k % 300 == 0: log('r45 corr', k, '/', len(todo), round(time.time() - t0), 's')
    log('R45_CORR_DONE')


# ---------- stats ----------
def paired(rows, z, seed, B):
    """analyse_sf 와 같은 점추정·같은 난수열의 상과 Exp(1) 무게로 세 판과 짝 차이를 한 번에."""
    Y = {v: ranks(np.array([r['main']['r_' + v] for r in rows])) for v in VERS}
    X_ = ranks(np.array([r['main']['ident'] for r in rows])); Zs = [ranks(np.array([r['main'][c] for r in rows])) for c in z]
    H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); cnt = np.bincount(si); base = 1.0 / cnt[si]
    est = {v: float(pcorr(Y[v], X_, Zs, base)) for v in VERS}; rng = np.random.default_rng(seed); bs = {v: [] for v in VERS}
    for _ in range(B):
        w = base * rng.exponential(1, len(sfl))[si]
        for v in VERS: bs[v].append(pcorr(Y[v], X_, Zs, w))
    bs = {v: np.array(bs[v]) for v in VERS}; q = lambda a: [float(x) for x in np.percentile(a, [2.5, 97.5])]
    out = dict(n=len(rows), n_sf=len(sfl), seed=seed, B=B)
    for v in VERS: out[v] = dict(est=est[v], lo=q(bs[v])[0], hi=q(bs[v])[1])
    for v in ('tilt', 'twist'): d = bs['full'] - bs[v]; out['full_minus_' + v] = dict(est=est['full'] - est[v], lo=q(d)[0], hi=q(d)[1])
    return out


def part_stats():
    new = {}
    for f in sorted(glob.glob(OUT45 + '/pairs_s*.jsonl')):
        for l in open(f):
            if l.strip(): r = json.loads(l); new.setdefault(tuple(r['ab']), r)
    base_rows = merged(OUT, OUT + '/anm_pairs.jsonl'); out = dict(written=time.strftime('%F %T'), n_new=len(new))
    ch = [r['check'] for r in new.values() if 'check' in r]; out['check_r_dyn'] = dict(n=len(ch), max_abs=float(max(abs(c) for c in ch)) if ch else None)
    def rows_for(z):   # analyse_sf 와 같은 거르기(r_dyn·ident·z) + 새 두 판; full 판은 저장 r_dyn 그대로(재계산 값과의 차는 check_r_dyn)
        rs = []
        for r in base_rows:
            m, n_ = r['main'], new.get(tuple(r['ab']), {})
            if all(m.get(c) is not None for c in ('r_dyn', 'ident') + z) and all(n_.get('r_' + v) is not None for v in VERS):
                rs.append(dict(ab=r['ab'], H=r['H'], main=dict(m, r_full=m['r_dyn'], r_tilt=n_['r_tilt'], r_twist=n_['r_twist'])))
        return rs
    rows1, rows2 = rows_for(('control',)), rows_for(('control', 'r_anm'))
    # 구현 점검: 등록 함수(analyse_sf, 저장 r_dyn)와 같은 행 수·같은 값, 그 값이 원고 값(셋째 자리)과 같음, 쌍별 재계산 r_dyn 차 < 1e-9
    ref1 = analyse_sf(base_rows, 'main', z=('control',), seed=SEED_C1); ref2 = analyse_sf(base_rows, 'main', z=('control', 'r_anm'), seed=SEED_C2)
    c1 = paired(rows1, ('control',), SEED_C1, NB); c2 = paired(rows2, ('control', 'r_anm'), SEED_C2, NB)
    pub = dict(C1=(0.225, 0.119, 0.333), C2=(0.181, 0.079, 0.283)); r3 = lambda d: (round(d['est'], 3), round(d['lo'], 3), round(d['hi'], 3))
    rep = dict(ref_C1=ref1, ref_C2=ref2, n_rows_C1=len(rows1), n_rows_C2=len(rows2), published=pub,
               same_rows=bool(ref1.get('n') == len(rows1) and ref2.get('n') == len(rows2)),
               same_values=bool(all(abs(c['full'][k] - r_[k]) < 1e-9 for c, r_ in ((c1, ref1), (c2, ref2)) for k in ('est', 'lo', 'hi'))),
               matches_published=bool(r3(ref1) == pub['C1'] and r3(ref2) == pub['C2']),
               r_dyn_recomputed=bool(out['check_r_dyn']['max_abs'] is not None and out['check_r_dyn']['max_abs'] < 1e-9))
    rep['reproduces'] = bool(rep['same_rows'] and rep['same_values'] and rep['matches_published'] and rep['r_dyn_recomputed'])
    out['implementation_check'] = rep; out['C1'] = c1; out['C2'] = c2
    t = c2['tilt']; dlt = c2['full_minus_tilt']; edge = bool(abs(t['lo']) <= EDGE)
    if edge: out['C2_second_seed'] = paired(rows2, ('control', 'r_anm'), SEED_2ND, 10000)
    if not rep['reproduces']: br, sent = 'not_interpreted', '구현 점검(전체 틀 판이 0.225 · 0.181 을 재현)이 서지 않아 나머지를 해석하지 않는다'
    elif t['lo'] > 0: br, sent = 'a', '연관은 결합 방향 요동에도 있고, 그 요동을 예측하는 탄성망을 통제한 뒤에도 남는다'
    elif dlt['lo'] > 0: br, sent = 'b', '연관의 상당 몫이 결합 방향 밖의 성분(펩타이드 면 비틀림)에 실려 있다'
    else: br, sent = 'c', '결합 방향 성분만으로는 판별할 해상도가 없다'
    out['verdict'] = dict(branch=br, boundary=edge, sentence=sent + (' (경계: f₁ 판 C-2 하한이 0 에서 ±0.01 안; 둘째 seed 는 서술)' if edge else ''))
    out['twist_neg_frac_median'] = float(np.median([r.get('twist_neg_frac', np.nan) for r in new.values() if 'twist_neg_frac' in r]))
    sp = OUT45 + '/sanity.jsonl'
    if os.path.exists(sp):
        S_ = [json.loads(l) for l in open(sp) if l.strip()]; S_ = [s for s in S_ if s.get('rho_full') is not None]
        out['sanity'] = {v: dict(n=len(S_), median=float(np.median([s['rho_' + v] for s in S_])), iqr=[float(np.percentile([s['rho_' + v] for s in S_], p)) for p in (25, 75)]) for v in VERS}
        out['sanity']['check_full_vs_stored_max_abs'] = float(max(abs(s['rho_full'] - s['rho_stored']) for s in S_ if s.get('rho_stored') is not None))
        out['sanity']['tilt_above_0_60'] = bool(out['sanity']['tilt']['median'] > 0.60)
    out['rule'] = ('등록 r45_f1: 구현 점검(전체 틀 판 = 등록 C-1 0.225 · C-2 0.181, 쌍별 r_dyn 차 < 1e-9)이 서지 않으면 해석하지 않음. 판정은 C-2: '
                   '(가) f₁ 판 하한 > 0; (나) f₁ 판 하한 ≤ 0 이고 짝 차이(full − tilt) 하한 > 0; (다) 그 밖. f₁ 판 하한이 0 에서 ±0.01 안이면 "경계" 로 보고, '
                   '둘째 seed 20261113 · 10,000 은 서술. 비틀림 판·탄성망 일치도는 서술.')
    jdump(out, RES + '/r45_f1.json', bk=True); log('R45_F1', out['verdict']); log('R45_STATS_DONE')


# ---------- sanity ----------
def part_sanity():
    cdom = json.load(open(OUT + '/pairs.json'))['domains']; fo = OUT45 + '/sanity.jsonl'; done = set()
    if os.path.exists(fo): done = {json.loads(l)['d'] for l in open(fo) if l.strip()}
    with open(fo, 'a') as fh:
        for k, d in enumerate(sorted(cdom)):
            if d in done: continue
            a = f'{OUT}/anm/{d}.npz'
            if not os.path.exists(a): continue
            try:
                za = np.load(a); pred = za['pred']; z = dom_get(d)
                if z is None: continue
                A = z['A']; pr = profiles(A, np.arange(A.shape[1])); st = float(za['rho']) if np.isfinite(float(za['rho'])) else None
                fh.write(json.dumps(dict(d=d, rho_stored=st, **{'rho_' + v: spearman(pred, pr[v]) for v in VERS})) + '\n')
            except Exception as e:
                log('ERR r45 sanity', d, repr(e)[:150])
            if k % 200 == 0: log('r45 sanity', k, '/', len(cdom))
    log('R45_SANITY_DONE')


if __name__ == '__main__':
    for part in PARTS45:
        log('부분 시작', part)
        try: globals()['part_' + part]()
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
        log('부분 끝', part)
    log('R45_F1_ALL_DONE', ','.join(PARTS45))
