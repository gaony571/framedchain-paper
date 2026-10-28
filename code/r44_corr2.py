"""(round44, 등록 = Drive results/prereg_2026-10-28_addendum4.json "r44_agree"·"r44_tmmean", 2026-10-28 06:46:46 UTC) 확인 실험 C 의 대응 일치(독립 기준·허용 폭)와 평균 점 TM-align 대응.
r44_corr_v2.py (f364e242) 의 함수(seq_map, tm_map, edge_pairs, dom_small, dom_get, gpa_profile, rot_share, spearman, analyse_sf, MINSH)를 그대로 쓴다.
대응: pibar (저장), tm_med (메도이드 Cα 의 TM-align — r44_corr_v2 의 tm 과 같은 계산), tm_mean (사영 평균 점 Pm 의 Cα, fs.ca_trace), seq (전역 BLOSUM62).
포함률 = 기준 짝 (i, j) 가운데 시험 대응이 결합 i 를 |j′ − j| ≤ t 인 j′ 에 짝짓는 비율 (t = 0, 1, 2; 퇴화 결합을 빼기 전의 짝). 동일성 = tm_med 잔기 정렬의 동일성.
tm_mean 대응 위에서 r_dyn·control·ident·r_ANM 을 다시 내고 C-1·C-2 (r38_c.analyse_sf, seed 20261099). 같은 쌍 판의 π̄·tm_med 값은 derived/r44/corr/pairs_v2_s*.jsonl.
부분: corr K k → stats   사용: python r44_corr2.py <ROOT> corr <K> <k> | stats   출력: derived/r44/corr2/pairs_c2_s<k>.jsonl, results/r44_corr2.json"""
import os, sys
ROOT = sys.argv[1]; PARTS2 = sys.argv[2].split(','); K2 = int(sys.argv[3]) if len(sys.argv) > 3 else 1; S2 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
_s = open(ROOT + '/code/r44_corr_v2.py').read(); _s = _s[:_s.index("# ---------- stats ----------")]
sys.argv = ['r44_corr_v2', ROOT, 'none', '1', '0']; exec(_s)   # r44_corr_v2 의 대응 함수와 (그 안의) r38_c·r38_p3 함수
OUT2 = ROOT + '/derived/r44/corr2'; os.makedirs(OUT2, exist_ok=True); SEED2 = 20261099; TOL = (0, 1, 2)
COMP = (('tm_med', 'seq'), ('tm_med', 'pibar'), ('tm_med', 'tm_mean'), ('tm_mean', 'seq'), ('pibar', 'seq'), ('pibar', 'tm_med'), ('pibar', 'tm_mean'))   # (기준, 시험)
BINS = (('lt15', 0.0, 0.15), ('15_20', 0.15, 0.20), ('20_30', 0.20, 0.30), ('ge30', 0.30, 1.01), ('all', 0.0, 1.01))


def cover(ref, test, t):
    """기준 짝 가운데 시험 대응이 같은 결합 i 를 ±t 안의 결합에 짝짓는 비율."""
    if not ref: return None
    T = dict(test)
    return float(np.mean([i in T and abs(T[i] - j) <= t for i, j in ref]))


def part_corr():
    P, PB = c_pairs(); fo = OUT2 + f'/pairs_c2_s{S2}.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    small, big, med, anm = {}, {}, {}, {}
    def get(d):
        if d not in small: small[d] = dom_small(d)[0]
        if d not in big:
            if len(big) > 40: big.pop(next(iter(big)))
            z = dom_get(d); big[d] = None if z is None else z['A']
        if d not in med:
            med[d] = np.asarray(fs.ca_trace(torch.as_tensor(np.load(f'{OUT}/medoid/{d}.npy').astype(np.float64))))
            a = f'{OUT}/anm/{d}.npz'; anm[d] = np.load(a)['pred'] if os.path.exists(a) else None
        return small[d], big[d], med[d], anm[d]
    allp = sorted(P, key=lambda p: (p[0], p[1])); blk = np.array_split(np.arange(len(allp)), K2)[S2]
    todo = [allp[i] for i in blk if (allp[i][0], allp[i][1]) not in done]; log('corr2 조각', S2, '/', K2, len(todo), '쌍 남음'); t0 = time.time()
    with open(fo, 'a') as fh:
        for k, p in enumerate(todo):
            try:
                (Da, Fa, ca_a, pa), (Db, Fb, ca_b, pb) = get(p[0]), get(p[1]); sa, sb = Da['seq'], Db['seq']
                Mmed, tm_med = tm_map(ca_a, ca_b, sa, sb)
                mca_a = np.asarray(fs.ca_trace(torch.as_tensor(np.asarray(Da['Pm'], dtype=np.float64)))); mca_b = np.asarray(fs.ca_trace(torch.as_tensor(np.asarray(Db['Pm'], dtype=np.float64))))
                Mmean, tm_mean = tm_map(mca_a, mca_b, sa, sb)
                EP = dict(pibar=[tuple(x) for x in np.array(PB[(p[0], p[1])]['best'], int).reshape(-1, 2).tolist()],
                          seq=edge_pairs(seq_map(sa, sb), len(sa)), tm_med=edge_pairs(Mmed, len(sa)), tm_mean=edge_pairs(Mmean, len(sa)))
                o = dict(ab=[p[0], p[1]], H=p[2], tm_score_med=tm_med, tm_score_mean=tm_mean, n_ep={c: len(v) for c, v in EP.items()},
                         ident_tm_res=float(np.mean([sa[i] == sb[j] for i, j in Mmed.items()])) if Mmed else None,
                         cover={f'{rf}_by_{ts}|t{t}': cover(EP[rf], EP[ts], t) for rf, ts in COMP for t in TOL})
                prs = np.array(EP['tm_mean'], int).reshape(-1, 2); keep = ~(Da['degen'][prs[:, 0]] | Db['degen'][prs[:, 1]]) if len(prs) else np.zeros(0, bool); sh = prs[keep]
                q = dict(n_pairs=int(len(prs)), n_shared=int(len(sh)))
                if len(sh) >= MINSH:
                    q.update(r_dyn=spearman(gpa_profile(Fa, sh[:, 0]), gpa_profile(Fb, sh[:, 1])), control=rot_share(Da['Pm'][sh[:, 0]], Db['Pm'][sh[:, 1]]),
                             ident=float(np.mean([sa[i] == sb[j] for i, j in sh])),
                             r_anm=spearman(pa[sh[:, 0]], pb[sh[:, 1]]) if pa is not None and pb is not None else None)
                o['tm_mean'] = q
                fh.write(json.dumps(o) + '\n')
            except Exception as e:
                log('ERR corr2', p[:2], repr(e)[:150])
            if k % 100 == 0: log('corr2', k, '/', len(todo), round(time.time() - t0), 's')
    log('R44_CORR2_DONE')


def part_stats():
    import glob
    rows = [json.loads(l) for f in sorted(glob.glob(OUT2 + '/pairs_c2_s*.jsonl')) for l in open(f) if l.strip()]
    seen = set(); rows = [r for r in rows if tuple(r['ab']) not in seen and not seen.add(tuple(r['ab']))]
    v2 = {tuple(r['ab']): r for f in sorted(glob.glob(OUT44 + '/pairs_v2_s*.jsonl')) for l in open(f) if l.strip() for r in [json.loads(l)]}
    out = dict(written=time.strftime('%F %T'), n_rows=len(rows))
    out['check_tm_med'] = dict(n=sum(tuple(r['ab']) in v2 for r in rows),
                               same_n_pairs=float(np.mean([r['n_ep']['tm_med'] == v2[tuple(r['ab'])]['tm']['n_pairs'] for r in rows if tuple(r['ab']) in v2])))
    ag = {}
    for rf, ts in COMP:
        for t in TOL:
            key = f'{rf}_by_{ts}|t{t}'; ag[key] = {}
            for nm, lo, hi in BINS:
                v = np.array([r['cover'][key] for r in rows if r['cover'].get(key) is not None and r.get('ident_tm_res') is not None and lo <= r['ident_tm_res'] < hi])
                ag[key][nm] = dict(n=int(len(v)), median=float(np.median(v)) if len(v) else None,
                                   q25=float(np.percentile(v, 25)) if len(v) else None, q75=float(np.percentile(v, 75)) if len(v) else None)
    out['agree'] = ag
    out['bin_counts'] = {nm: ag['tm_med_by_seq|t0'][nm]['n'] for nm, _, _ in BINS}
    sq, pb = ag['tm_med_by_seq|t2']['lt15']['median'], ag['tm_med_by_pibar|t2']['lt15']['median']
    out['predictions_agree'] = dict(pa_seq_covers_tm_lt02=bool(sq is not None and sq < 0.2), pb_pibar_covers_tm_more=bool(sq is not None and pb is not None and pb > sq))
    full = dict(C1=analyse_sf(rows, key='tm_mean', z=('control',), seed=SEED2), C2=analyse_sf(rows, key='tm_mean', z=('control', 'r_anm'), seed=SEED2 + 1),
                n_shared_median=float(np.median([r['tm_mean']['n_shared'] for r in rows])))
    ok = lambda r, c: r.get(c) and r[c].get('r_dyn') is not None and r[c].get('r_anm') is not None and r[c].get('control') is not None
    com = []
    for r in rows:
        w = v2.get(tuple(r['ab']))
        if w is not None and ok(w, 'pibar') and ok(w, 'tm') and ok(r, 'tm_mean'):
            com.append(dict(ab=r['ab'], H=r['H'], pibar=w['pibar'], tm_med=w['tm'], tm_mean=r['tm_mean']))
    cm = {c: dict(C1=analyse_sf(com, key=c, z=('control',), seed=SEED2 + 10 + 2 * i), C2=analyse_sf(com, key=c, z=('control', 'r_anm'), seed=SEED2 + 11 + 2 * i))
          for i, c in enumerate(('pibar', 'tm_med', 'tm_mean'))}
    a, b, m = cm['pibar']['C1'].get('est'), cm['tm_med']['C1'].get('est'), cm['tm_mean']['C1'].get('est')
    share = float((m - b) / (a - b)) if None not in (a, b, m) and abs(a - b) > 1e-9 else None
    reading = None if share is None else ('주로 정렬 방법' if share < 1 / 3 else '주로 정렬한 구조' if share > 2 / 3 else '둘 다')
    out['tmmean'] = dict(full=full, common=dict(n=len(com), **cm), share=share, reading=reading,
                         robust=bool(full['C1'].get('lo', -1) > 0 and full['C2'].get('lo', -1) > 0))
    out['rule'] = ('등록 r44_agree: PA = TM_med 동일성 < 15% 에서 seq 가 TM_med 를 담는 비율(t = 2) 중앙 < 0.2, PB = 같은 쌍에서 π̄ 의 것 > seq 의 것. '
                   '등록 r44_tmmean: 몫 = (C1_TMmean − C1_TMmed)/(C1_π̄ − C1_TMmed) (같은 쌍) < 1/3 정렬 방법, > 2/3 정렬한 구조, 사이 둘 다; 섬 = C-1·C-2 하한 > 0.')
    jdump(out, RES + '/r44_corr2.json', bk=True); log('R44_CORR2', out['predictions_agree'], out['tmmean']['reading'], out['tmmean']['robust']); log('R44_CORR2_STATS_DONE')


if __name__ == '__main__':
    for part in PARTS2:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('R44_CORR2_ALL_DONE', ','.join(PARTS2))
