"""(등록 = Drive results/prereg_2026-10-30_addendum9.json "r46_align") 확인 실험 C: 서열 동일성과 r_dyn 이 같은 대응(π̄)을 공유해 정렬 정확도가 공통 원인일 수 있는가.
r_dyn·control·r_ANM 은 등록 판 그대로(π̄ 공유 결합 = π̄ 짝에서 퇴화 결합을 뺀 것; r38_c.merged 의 main). 독립 기준 = 메도이드 Cα 의 TM-align(실제 서열)을
다시 돌림(r44_corr_v2.tm_map 그대로; 정렬 잔기 = TM-align 출력에서 두 쪽 모두 틈이 아닌 모든 자리, 거리 문턱 없음).
  A  (판정; 같은 대응 오류가 동일성과 r_dyn 에 함께 들어가는 기계적 결합을 끊음): 동일성만 TM-align 정렬 잔기 전체의 동일성으로. TM 정렬 잔기 ≥ 30.
  B  (판정; "분기 → π̄ 오류 → r_dyn" 경로를 통제): 정밀도 = r_dyn 이 쓴 π̄ 공유 결합 가운데 TM-align 이 같은 결합 i 를 ±k 결합 안의 j′ 에 짝짓는 비율 (k = 0, 2)
      을 C-1·C-2 통제에 더함. B 판정 = B0·B2 가운데 약한 쪽. 0 쪽 치우침: π̄ 는 요동만큼 줄어든 평균 점의 정렬이라 정밀도에 동역학 닮음도 들어감.
  서술: A′ (자리는 π̄ 공유 결합, 짝만 TM-align 잔기 대응; TM 짝 있는 결합 ≥ 30), 재현율 판 B (r44_corr2 cover), 결합 경로 편상관, 둘째 seed.
구간 = 상과 Exp(1) 군집 부트스트랩 10,000, seed 20261121 한 흐름에서 모든 판을 같은 무게로 함께(짝지은 차 여섯 = {A, B0, B2} × {C-1, C-2}). 판정은 이 흐름 한 번.
점검 (1, tm 부분 끝): 다시 낸 TM-align 이 r44_corr2 와 정수 셈으로 완전히 같음 — TM 결합 짝 수, 동일 잔기 수(= ident_tm_res × 정렬 잔기 수), π̄ 가 담는 TM 짝 수(±0, ±2).
      모든 공통 쌍이 같으면 derived/r46/tm_check.json 에 pass = true. stats 는 pass = true 일 때만 돈다(어긋나면 통계 전에 고침 기록을 따로 씀).
점검 (2, stats 처음): 전체 C 쌍의 기준 판을 등록 seed(20261071, 20261074)로 다시 내어 results/r38_c.json C1·C2 main (est, lo, hi)와 |차| < 1e-9.
부분: tm → stats   사용: python r46_align.py <ROOT> tm|stats|tm,stats   출력: derived/r46/align_pairs.jsonl, derived/r46/tm_check.json, results/r46_align.json"""
import sys, os, json, glob, time, traceback
import numpy as np
ROOT = sys.argv[1]; PARTS46 = sys.argv[2].split(',')
_s = open(ROOT + '/code/r44_corr_v2.py').read(); _s = _s[:_s.index("# ---------- stats ----------")]
sys.argv = ['r44_corr_v2', ROOT, 'none', '1', '0']; exec(_s)   # r44_corr_v2 의 tm_map·edge_pairs·c_pairs 와 (그 안의) r38_c·r38_p3 함수: merged, analyse_sf, ranks, pcorr, spearman, dom_small, log
OUT46 = ROOT + '/derived/r46'; os.makedirs(OUT46, exist_ok=True); FP = OUT46 + '/align_pairs.jsonl'; FC = OUT46 + '/tm_check.json'
SEED = 20261121; SEED2 = 20261122; NBOOT = 10000; MIN_RES = 30; KS = (0, 2)
ORDER = {'fails': 0, 'partly': 1, 'holds': 2}


def corr2_rows():
    out = {}
    for f in sorted(glob.glob(ROOT + '/derived/r44/corr2/pairs_c2_s*.jsonl')):
        for l in open(f):
            if l.strip():
                r = json.loads(l); out.setdefault(tuple(r['ab']), r)
    return out


def load_fp():
    out = {}
    for l in open(FP):
        if l.strip():
            r = json.loads(l); out.setdefault(tuple(r['ab']), r)
    return out


def part_tm():
    P, PB = c_pairs(); done = set()
    if os.path.exists(FP): done = set(load_fp())
    small, med = {}, {}
    def get(d):
        if d not in small:
            small[d] = dom_small(d)[0]
            med[d] = np.asarray(fs.ca_trace(torch.as_tensor(np.load(f'{OUT}/medoid/{d}.npy').astype(np.float64))))
        return small[d], med[d]
    todo = [p for p in sorted(P, key=lambda p: (p[0], p[1])) if (p[0], p[1]) not in done]; log('r46 tm', len(todo), '쌍 남음'); t0 = time.time()
    with open(FP, 'a') as fh:
        for k, p in enumerate(todo):
            try:
                (Da, ca_a), (Db, ca_b) = get(p[0]), get(p[1]); sa, sb = Da['seq'], Db['seq']
                M, tm = tm_map(ca_a, ca_b, sa, sb); TP = edge_pairs(M, len(sa)); T = dict(TP)
                allp = [tuple(x) for x in np.array(PB[(p[0], p[1])]['best'], int).reshape(-1, 2).tolist()]; PM = dict(allp)
                prs = np.array(allp, int).reshape(-1, 2)
                keep = ~(Da['degen'][prs[:, 0]] | Db['degen'][prs[:, 1]]) if len(prs) else np.zeros(0, bool); sh = prs[keep]
                o = dict(ab=[p[0], p[1]], H=p[2], tm=tm, n_res=len(M), k_ident=int(sum(sa[i] == sb[j] for i, j in M.items())), n_tm_bonds=len(TP), n_shared=int(len(sh)),
                         cnt_rec={str(t): int(sum(i in PM and abs(PM[i] - j) <= t for i, j in TP)) for t in KS})
                o['ident_tm_res'] = o['k_ident'] / o['n_res'] if o['n_res'] else None
                for kk in KS: o[f'prec{kk}'] = float(np.mean([i in T and abs(T[i] - j) <= kk for i, j in sh])) if len(sh) else None
                ap = [sa[i] == sb[M[i]] for i, j in sh if i in M]; o['n_aprime'] = len(ap); o['ident_aprime'] = float(np.mean(ap)) if ap else None
                fh.write(json.dumps(o) + '\n')
            except Exception as e:
                log('ERR r46 tm', p[:2], repr(e)[:150])
            if k % 200 == 0: log('r46 tm', k, '/', len(todo), round(time.time() - t0), 's')
    # 점검 (1): 정수 셈의 완전 일치
    AL, c2 = load_fp(), corr2_rows(); bad, n_chk, miss = [], 0, 0
    for key, a in AL.items():
        q = c2.get(key)
        if q is None: miss += 1; continue
        n_chk += 1; cv = q.get('cover') or {}
        ok = a['n_tm_bonds'] == q['n_ep']['tm_med']
        ok &= (q.get('ident_tm_res') is None and a['n_res'] == 0) or (a['n_res'] > 0 and q.get('ident_tm_res') is not None and round(q['ident_tm_res'] * a['n_res']) == a['k_ident'])
        for t in KS:
            f = cv.get(f'tm_med_by_pibar|t{t}')
            ok &= (f is None and a['n_tm_bonds'] == 0) or (f is not None and round(f * a['n_tm_bonds']) == a['cnt_rec'][str(t)])
        if not ok: bad.append(list(key))
    chk = dict(written=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), n_target=len(P), n_pairs=len(AL), n_checked=n_chk, n_without_corr2=miss, n_mismatch=len(bad), mismatch=bad[:50],
               max_abs_tm_diff=float(max([abs(AL[k]['tm'] - c2[k]['tm_score_med']) for k in AL if k in c2] or [0.0])), pass_=bool(n_chk > 0 and not bad))
    chk['pass'] = chk.pop('pass_'); json.dump(chk, open(FC, 'w'), ensure_ascii=False, indent=1); log('R46_TM_CHECK', {k: chk[k] for k in ('n_pairs', 'n_checked', 'n_mismatch', 'pass')}); log('R46_TM_DONE')


def boot(rows, specs, seed, nboot):
    """specs: [(이름, key, y, x, z)] — 같은 쌍에서 seed 한 흐름의 상과 Exp(1) 무게로 모든 판을 함께 뽑는다."""
    H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); cnt = np.bincount(si); w0 = 1.0 / cnt[si]
    A = {nm: (ranks(np.array([r[key][y] for r in rows])), ranks(np.array([r[key][x] for r in rows])), [ranks(np.array([r[key][c] for r in rows])) for c in z])
         for nm, key, y, x, z in specs}
    est = {nm: float(pcorr(Y, X, Z, w0)) for nm, (Y, X, Z) in A.items()}
    rng = np.random.default_rng(seed); bs = {nm: np.empty(nboot) for nm in A}
    for b in range(nboot):
        w = w0 * rng.exponential(1, len(sfl))[si]
        for nm, (Y, X, Z) in A.items(): bs[nm][b] = pcorr(Y, X, Z, w)
    return est, bs, len(sfl)


def ci(v): return [float(x) for x in np.percentile(v, [2.5, 97.5])]


def verdict(c1, c2):
    if c1['lo'] > 0 and c2['lo'] > 0: return 'holds'
    if c1['lo'] > 0: return 'partly'
    return 'fails'


def part_stats():
    if not os.path.exists(FC) or not json.load(open(FC)).get('pass'):
        log('R46 STOP: tm_check 통과 파일 없음 — 통계를 내지 않음'); return
    base = merged(OUT, OUT + '/anm_pairs.jsonl'); reg = json.load(open(RES + '/r38_c.json'))
    rc1 = analyse_sf(base, 'main', z=('control',), seed=20261071); rc2 = analyse_sf(base, 'main', z=('control', 'r_anm'), seed=20261074)
    dev = max([abs(rc1[k] - reg['C1']['main'][k]) for k in ('est', 'lo', 'hi')] + [abs(rc2[k] - reg['C2']['main'][k]) for k in ('est', 'lo', 'hi')])
    out = dict(written=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), check=dict(registered_dev=float(dev), registered_ok=bool(dev < 1e-9), n_base=len(base), tm_check=json.load(open(FC))))
    if dev >= 1e-9:
        out['stopped'] = '등록값 재현 실패 — 판정을 내지 않음'; json.dump(out, open(RES + '/r46_align.json', 'x'), ensure_ascii=False, indent=1); log('R46 STOP', out['check']['registered_dev']); return
    AL, c2r = load_fp(), corr2_rows(); rows = []
    for r in base:
        m = r['main']; key = tuple(r['ab']); a = AL.get(key)
        if a is None or key not in c2r or not all(m.get(c) is not None for c in ('r_dyn', 'ident', 'control', 'r_anm')): continue
        if a['n_res'] < MIN_RES or a['prec0'] is None: continue
        cv = c2r[key].get('cover') or {}
        rows.append(dict(ab=r['ab'], H=r['H'], main=dict(m, prec0=a['prec0'], prec2=a['prec2'], rec0=cv.get('tm_med_by_pibar|t0'), rec2=cv.get('tm_med_by_pibar|t2')),
                         alt=dict(r_dyn=m['r_dyn'], control=m['control'], r_anm=m['r_anm'], ident=a['ident_tm_res']),
                         ap=dict(r_dyn=m['r_dyn'], control=m['control'], r_anm=m['r_anm'], ident=a['ident_aprime'], n=a['n_aprime'])))
    Z1, Z2 = ('control',), ('control', 'r_anm')
    specs = [('base_C1', 'main', 'r_dyn', 'ident', Z1), ('base_C2', 'main', 'r_dyn', 'ident', Z2),
             ('A_C1', 'alt', 'r_dyn', 'ident', Z1), ('A_C2', 'alt', 'r_dyn', 'ident', Z2),
             ('B0_C1', 'main', 'r_dyn', 'ident', Z1 + ('prec0',)), ('B0_C2', 'main', 'r_dyn', 'ident', Z2 + ('prec0',)),
             ('B2_C1', 'main', 'r_dyn', 'ident', Z1 + ('prec2',)), ('B2_C2', 'main', 'r_dyn', 'ident', Z2 + ('prec2',)),
             ('cpl_rdyn_prec0', 'main', 'r_dyn', 'prec0', Z2), ('cpl_ident_prec0', 'main', 'ident', 'prec0', Z2),
             ('cpl_rdyn_prec2', 'main', 'r_dyn', 'prec2', Z2), ('cpl_ident_prec2', 'main', 'ident', 'prec2', Z2)]
    est, bs, nsf = boot(rows, specs, SEED, NBOOT)
    res = {nm: dict(est=est[nm], lo=ci(bs[nm])[0], hi=ci(bs[nm])[1]) for nm in est}
    diff = {f'base-{v}_{c}': dict(est=est[f'base_{c}'] - est[f'{v}_{c}'], lo=ci(bs[f'base_{c}'] - bs[f'{v}_{c}'])[0], hi=ci(bs[f'base_{c}'] - bs[f'{v}_{c}'])[1])
            for v in ('A', 'B0', 'B2') for c in ('C1', 'C2')}
    est2, bs2, _ = boot(rows, specs[:8], SEED2, NBOOT)
    seed2 = {nm: dict(lo=ci(bs2[nm])[0], hi=ci(bs2[nm])[1]) for nm in est2}
    rows_ap = [r for r in rows if r['ap']['n'] >= MIN_RES and r['ap']['ident'] is not None]
    rows_rc = [r for r in rows if r['main']['rec0'] is not None and r['main']['rec2'] is not None]
    desc_an = dict(Aprime_C1=analyse_sf(rows_ap, 'ap', z=Z1, seed=SEED + 10, B=NBOOT), Aprime_C2=analyse_sf(rows_ap, 'ap', z=Z2, seed=SEED + 11, B=NBOOT),
                   recall0_C1=analyse_sf(rows_rc, 'main', z=Z1 + ('rec0',), seed=SEED + 12, B=NBOOT), recall0_C2=analyse_sf(rows_rc, 'main', z=Z2 + ('rec0',), seed=SEED + 13, B=NBOOT),
                   recall2_C1=analyse_sf(rows_rc, 'main', z=Z1 + ('rec2',), seed=SEED + 14, B=NBOOT), recall2_C2=analyse_sf(rows_rc, 'main', z=Z2 + ('rec2',), seed=SEED + 15, B=NBOOT))
    med = lambda k, key='main': float(np.median([r[key][k] for r in rows]))
    desc = dict(n_pairs=len(rows), n_sf=nsf, n_aprime=len(rows_ap), n_recall=len(rows_rc),
                ident_pibar_median=med('ident'), ident_tm_res_median=med('ident', 'alt'), prec0_median=med('prec0'), prec2_median=med('prec2'),
                spearman_ident_pibar_tm=float(spearman(np.array([r['main']['ident'] for r in rows]), np.array([r['alt']['ident'] for r in rows]))))
    vA = verdict(res['A_C1'], res['A_C2']); vB0 = verdict(res['B0_C1'], res['B0_C2']); vB2 = verdict(res['B2_C1'], res['B2_C2'])
    vB = min((vB0, vB2), key=ORDER.get)
    reading = 'both_hold' if vA == vB == 'holds' else 'both_fail' if vA == vB == 'fails' else 'mixed'
    not_all = any(d['lo'] > 0 for d in diff.values())   # 짝지은 차 여섯 모두가 대상
    if reading == 'both_hold': case = 'hold_not_all' if not_all else 'hold'
    elif reading == 'both_fail': case = 'fail'
    elif vB == 'fails': case = 'A_ok_B_fail'          # A 는 적어도 partly
    elif vA == 'fails': case = 'A_fail_B_ok'          # B 는 적어도 partly
    else: case = 'mixed_partly'                       # 둘 다 적어도 partly, 하나 이상 partly
    out.update(desc=desc, est=res, paired_diff=diff, seed2=seed2, descriptive=desc_an,
               verdict=dict(A=vA, B0=vB0, B2=vB2, B=vB, reading=reading, not_all=bool(not_all), sentence_case=case, abstract_change=bool(vA == 'fails')),
               rule='등록 r46_align: 판마다 C-1·C-2 하한 > 0 → holds, C-1 만 → partly, 아니면 fails; B = B0·B2 중 약한 쪽. 문장 갈래: hold (짝지은 차 여섯 가운데 하나라도 하한 > 0 이면 hold_not_all), '
                    'fail, A_ok_B_fail (B 의 0 쪽 치우침 문장), A_fail_B_ok, mixed_partly. 초록은 A 가 fails 일 때만 단서(B 와 무관). 판정은 seed 20261121 의 10,000 회 한 흐름. 둘째 seed·A′·재현율 판·결합 편상관은 서술.')
    json.dump(out, open(RES + '/r46_align.json', 'x'), ensure_ascii=False, indent=1)
    log('R46', out['verdict']); log('R46_STATS_DONE')


if __name__ == '__main__':
    for part in PARTS46:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
