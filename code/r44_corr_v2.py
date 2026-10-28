"""(round44, 등록 = Drive results/prereg_2026-10-28_addendum.json "r44_corr" + addendum2 "r44_corr_v2"·"r44_gp", 2026-10-28 UTC) 확인 실험 C 의 대응 바꾸기 — 나눠 돌리는 판.
r44_corr.py (4c23a80c) 와 계산 식은 같다. 바뀐 것: (1) corr 를 K 조각으로 나눠 돌림(성능만), (2) 넷째 대응 pibar_noGP (π̄ 공유 결합 가운데 두 쪽 네 잔기에 Gly·Pro 가 없는 것; 등록 r44_gp),
(3) 대응마다 Gly/Pro 공변량 gp_match (같은 자리 둘 다 G 또는 둘 다 P 인 비율)·gp_any (한쪽이라도 G/P 인 비율), (4) stats 가 조각 파일을 합치고 r44_gp 판정을 더함.
공유 결합을 정하는 대응을 셋으로 둔다:
  pibar : 등록 판. 평균 점 D 정렬(r38_c 의 pibar.jsonl 'best').
  seq   : 서열 정렬. Biopython PairwiseAligner 전역, BLOSUM62, 틈 열기 −10 · 늘이기 −0.5, 끝 틈 0 (EMBOSS needle 기본과 같은 꼴).
  tm    : TM-align(tmtools, 메도이드 Cα, 실제 서열). 두 쪽 모두 틈이 아닌 자리를 정렬로 본다.
잔기 정렬 → 결합 짝: 결합 i (잔기 i → i+1) 와 결합 j 는 잔기 i↔j, i+1↔j+1 이 모두 정렬될 때 짝이다. 퇴화 결합은 r38_p3.one_pair 처럼 뺀다.
쌍마다 r38_p3.one_pair 와 같은 r_dyn·control(d_fix 짝 몫)·ident, r38_c.anm_pairs 와 같은 r_ANM·BLOSUM62, 그리고 대응 사이 일치(π̄ 짝 가운데 다른 대응에도 있는 비율, Jaccard).
부분: prep → corr → stats
  prep : C 도메인 앙상블을 r38_c.dom_build(같은 seed)로 다시 만들어 로컬 캐시(/content/ccache)에 둔다(앞 런타임의 캐시는 없음). 작업자 여럿.
  corr : C 쌍(derived/r38/c/pairs.json 가운데 pibar 가 있는 것)마다 세 대응. 점검: π̄ 판 r_dyn·control 이 저장된 pairs_out.jsonl 의 main 과 같은지(모든 쌍).
  stats: 대응마다 C-1(z = control), C-2(z = control, r_ANM) — r38_c.analyse_sf 그대로(상과 Exp(1) 군집 부트스트랩 2,000). 세 대응 모두 공유 결합 ≥ 30 인 같은 쌍 집합에서도.
사용: python r44_corr_v2.py <ROOT> corr <K> <k> | stats   출력: derived/r44/corr/pairs_v2_s<k>.jsonl, results/r44_corr.json (사본 results/_backup/<날짜>/)"""
import os, sys
ROOT = sys.argv[1]; PARTS44 = sys.argv[2].split(','); K44 = int(sys.argv[3]) if len(sys.argv) > 3 else 1; S44 = int(sys.argv[4]) if len(sys.argv) > 4 else 0; NW44 = 1
_src = open(ROOT + '/code/r38_c.py').read(); _src = _src[:_src.index("ALLC = ['pairs'")]
sys.argv = ['r38_c', ROOT, 'none', str(NW44)]; exec(_src)   # r38_c (그 안에서 r38_p3) 의 함수: dom_get, dom_small, gpa_profile, rot_share, spearman, analyse_sf, blosum, MINSH
OUT44 = ROOT + '/derived/r44/corr'; os.makedirs(OUT44, exist_ok=True)
CORRS = ('pibar', 'seq', 'tm', 'pibar_noGP'); SEED44 = 20261094; GP = set('GP')


# ---------- 대응 ----------
_aligner = None
def seq_map(sa, sb):
    global _aligner
    from Bio import Align
    from Bio.Align import substitution_matrices
    if _aligner is None:
        a = Align.PairwiseAligner(); a.mode = 'global'; a.substitution_matrix = substitution_matrices.load('BLOSUM62')
        a.open_gap_score = -10; a.extend_gap_score = -0.5
        for k in ('target_end_gap_score', 'query_end_gap_score'): setattr(a, k, 0.0)
        _aligner = a
    ok = set(_aligner.substitution_matrix.alphabet); fix = lambda s: ''.join(c if c in ok else 'X' for c in s)   # 표준 밖 글자는 X (자리 수는 그대로)
    aln = _aligner.align(fix(sa), fix(sb))[0]; tb, qb = aln.aligned; M = {}
    for (t0, t1), (q0, q1) in zip(tb, qb):
        for k in range(int(t1) - int(t0)): M[int(t0) + k] = int(q0) + k
    return M


def tm_map(ca_a, ca_b, sa, sb):
    from tmtools import tm_align
    r = tm_align(np.ascontiguousarray(ca_a, dtype=np.float64), np.ascontiguousarray(ca_b, dtype=np.float64), sa, sb); M = {}; i = j = 0
    if not hasattr(r, 'seqxA'): raise RuntimeError('tmtools 결과에 정렬 문자열(seqxA)이 없음')
    for x, y in zip(r.seqxA, r.seqyA):
        if x != '-' and y != '-': M[i] = j
        i += x != '-'; j += y != '-'
    return M, float(max(r.tm_norm_chain1, r.tm_norm_chain2))


def edge_pairs(M, La):
    return [(i, M[i]) for i in range(La - 1) if i in M and M.get(i + 1) == M[i] + 1]


# ---------- prep ----------
def _prep_worker(d):
    try:
        z = dom_get(d); return d, 'ok' if z is not None else 'none'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:120]


def c_pairs():
    P = json.load(open(OUT + '/pairs.json'))['pairs']; PB = {}
    for l in open(OUT + '/pibar.jsonl'):
        r = json.loads(l); PB[tuple(r['ab'])] = r
    P = [p for p in P if (p[0], p[1]) in PB and dom_small(p[0])[0] is not None and dom_small(p[1])[0] is not None]
    return P, PB


def part_prep():
    from multiprocessing import Pool
    P, _ = c_pairs(); doms = sorted({x for p in P for x in p[:2]}); cnt = {}
    with Pool(NW44) as pool:
        for k, (d, st) in enumerate(pool.imap_unordered(_prep_worker, doms)):
            cnt[st[:12]] = cnt.get(st[:12], 0) + 1
            if k % 100 == 0 or not st.startswith('ok'): log('prep', k, '/', len(doms), d, st)
    log('R44_PREP', cnt); log('R44_PREP_DONE')


# ---------- corr ----------
def part_corr():
    P, PB = c_pairs(); fo = OUT44 + f'/pairs_v2_s{S44}.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    stored = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(OUT + '/pairs_out.jsonl')}
    B = blosum(); small, big, med, anm = {}, {}, {}, {}
    def get(d):
        if d not in small: small[d] = dom_small(d)[0]
        if d not in big:
            if len(big) > 40: big.pop(next(iter(big)))
            z = dom_get(d); big[d] = None if z is None else z['A']
        if d not in med:
            med[d] = np.asarray(fs.ca_trace(torch.as_tensor(np.load(f'{OUT}/medoid/{d}.npy').astype(np.float64))))
            a = f'{OUT}/anm/{d}.npz'; anm[d] = np.load(a)['pred'] if os.path.exists(a) else None
        return small[d], big[d], med[d], anm[d]
    allp = sorted(P, key=lambda p: (p[0], p[1])); blk = np.array_split(np.arange(len(allp)), K44)[S44]
    todo = [allp[i] for i in blk if (allp[i][0], allp[i][1]) not in done]; log('corr 조각', S44, '/', K44, len(todo), '쌍 남음'); t0 = time.time()
    with open(fo, 'a') as fh:
        for k, p in enumerate(todo):
            try:
                (Da, Fa, ca_a, pa), (Db, Fb, ca_b, pb) = get(p[0]), get(p[1]); sa, sb = Da['seq'], Db['seq']
                tmM, tm = tm_map(ca_a, ca_b, sa, sb)
                EP = dict(pibar=[tuple(x) for x in np.array(PB[(p[0], p[1])]['best'], int).reshape(-1, 2).tolist()],
                          seq=edge_pairs(seq_map(sa, sb), len(sa)), tm=edge_pairs(tmM, len(sa)))
                EP['pibar_noGP'] = [(i, j) for i, j in EP['pibar'] if not ({sa[i], sa[i + 1], sb[j], sb[j + 1]} & GP)]
                o = dict(ab=[p[0], p[1]], H=p[2], tm_med=tm)
                for c in CORRS:
                    prs = np.array(EP[c], int).reshape(-1, 2); keep = ~(Da['degen'][prs[:, 0]] | Db['degen'][prs[:, 1]]) if len(prs) else np.zeros(0, bool); sh = prs[keep]
                    q = dict(n_pairs=int(len(prs)), n_shared=int(len(sh)))
                    if len(sh) >= MINSH:
                        q.update(r_dyn=spearman(gpa_profile(Fa, sh[:, 0]), gpa_profile(Fb, sh[:, 1])), control=rot_share(Da['Pm'][sh[:, 0]], Db['Pm'][sh[:, 1]]),
                                 ident=float(np.mean([sa[i] == sb[j] for i, j in sh])),
                                 gp_match=float(np.mean([sa[i] == sb[j] and sa[i] in GP for i, j in sh])), gp_any=float(np.mean([sa[i] in GP or sb[j] in GP for i, j in sh])),
                                 r_anm=spearman(pa[sh[:, 0]], pb[sh[:, 1]]) if pa is not None and pb is not None else None,
                                 blosum=float(np.mean([B[sa[i], sb[j]] if sa[i] in B.alphabet and sb[j] in B.alphabet else 0.0 for i, j in sh])) if B is not False else None)
                    if c not in ('pibar', 'pibar_noGP'):
                        Sp, Sx = set(EP['pibar']), set(EP[c]); q.update(frac_pibar_in=len(Sp & Sx) / len(Sp) if Sp else None, jaccard=len(Sp & Sx) / len(Sp | Sx) if Sp | Sx else None)
                    o[c] = q
                s = stored.get((p[0], p[1]), {}).get('main', {})
                if s.get('r_dyn') is not None and o['pibar'].get('r_dyn') is not None:
                    o['check'] = dict(r_dyn=o['pibar']['r_dyn'] - s['r_dyn'], control=o['pibar']['control'] - s['control'])
                fh.write(json.dumps(o) + '\n')
            except Exception as e:
                log('ERR corr', p[:2], repr(e)[:150])
            if k % 300 == 0: log('corr', k, '/', len(todo), round(time.time() - t0), 's')
    log('R44_CORR_DONE')


# ---------- stats ----------
def part_stats():
    rows = [json.loads(l) for f in sorted(glob.glob(OUT44 + '/pairs_v2_s*.jsonl')) for l in open(f) if l.strip()]
    seen = set(); rows = [r for r in rows if tuple(r['ab']) not in seen and not seen.add(tuple(r['ab']))]; out = dict(written=time.strftime('%F %T'), n_rows=len(rows), source='pairs_v2_s*.jsonl')
    ch = [r['check'] for r in rows if 'check' in r]
    out['check'] = dict(n=len(ch), max_abs_r_dyn=float(max(abs(c['r_dyn']) for c in ch)) if ch else None, max_abs_control=float(max(abs(c['control']) for c in ch)) if ch else None)
    full = {}; common = [r for r in rows if all(r.get(c, {}).get('r_dyn') is not None and r[c].get('r_anm') is not None for c in CORRS)]
    for i, c in enumerate(CORRS):
        full[c] = dict(C1=analyse_sf(rows, key=c, z=('control',), seed=SEED44 + 10 * i), C2=analyse_sf(rows, key=c, z=('control', 'r_anm'), seed=SEED44 + 10 * i + 1),
                       C1_common=analyse_sf(common, key=c, z=('control',), seed=SEED44 + 10 * i + 2), C2_common=analyse_sf(common, key=c, z=('control', 'r_anm'), seed=SEED44 + 10 * i + 3),
                       n_shared_median=float(np.median([r[c]['n_shared'] for r in rows if c in r])), ident_median=float(np.median([r[c]['ident'] for r in rows if r.get(c, {}).get('ident') is not None])))
        if c in ('seq', 'tm'):
            fr = [r[c]['frac_pibar_in'] for r in rows if r.get(c, {}).get('frac_pibar_in') is not None]; jc = [r[c]['jaccard'] for r in rows if r.get(c, {}).get('jaccard') is not None]
            full[c]['agree'] = dict(frac_pibar_in_median=float(np.median(fr)), frac_pibar_in_iqr=[float(np.percentile(fr, 25)), float(np.percentile(fr, 75))], jaccard_median=float(np.median(jc)))
    out['by_corr'] = full; out['n_common'] = len(common)
    ok = {c: bool(full[c]['C1'].get('lo', -1) > 0 and full[c]['C2'].get('lo', -1) > 0) for c in ('seq', 'tm')}
    out['robust'] = dict(seq=ok['seq'], tm=ok['tm'], both=bool(ok['seq'] and ok['tm']))
    out['predictions'] = dict(p1_tm_agrees_more=bool(full['tm']['agree']['frac_pibar_in_median'] > full['seq']['agree']['frac_pibar_in_median']), p2_tm_robust=ok['tm'], p3_seq_robust=ok['seq'],
                              p3b_seq_smaller=bool(full['seq']['C1'].get('est', 9) < full['pibar']['C1'].get('est', -9)))
    out['sentence'] = ('C 의 효과는 공유 결합을 서열 정렬·TM-align 대응으로 정해도 남는다' if ok['seq'] and ok['tm'] else
                       f"대응을 바꾸면: 서열 정렬 {'섬' if ok['seq'] else '안 섬'}, TM-align {'섬' if ok['tm'] else '안 섬'}")
    gp = dict(cov_C1=analyse_sf(rows, key='pibar', z=('control', 'gp_match', 'gp_any'), seed=SEED44 + 100), cov_C2=analyse_sf(rows, key='pibar', z=('control', 'r_anm', 'gp_match', 'gp_any'), seed=SEED44 + 101),
              excl_C1=full['pibar_noGP']['C1'], excl_C2=full['pibar_noGP']['C2'], ref_C1=full['pibar']['C1'], ref_C2=full['pibar']['C2'],
              gp_match_median=float(np.median([r['pibar']['gp_match'] for r in rows if r.get('pibar', {}).get('gp_match') is not None])))
    gp['remains'] = bool(all(gp[k].get('lo', -1) > 0 for k in ('cov_C1', 'cov_C2', 'excl_C1', 'excl_C2')))
    gp['halved'] = bool(any(gp[k].get('est', 9) <= 0.5 * gp[r_].get('est', 0) for k, r_ in (('cov_C1', 'ref_C1'), ('excl_C1', 'ref_C1'))))
    gp['sentence'] = ('서열 동일성의 효과는 Gly·Pro 일치를 통제하거나 그 결합을 빼도 남는다' if gp['remains'] and not gp['halved'] else
                      'Gly·Pro 일치가 효과의 상당 부분을 나른다' if gp['halved'] else 'Gly·Pro 통제 뒤 효과의 하한이 0 에 닿는다')
    out['r44_gp'] = gp
    out['rule'] = ('등록 r44_corr: 대응마다 C-1·C-2 는 r38_c.analyse_sf 그대로. 선다 = 하한 > 0 (C-1·C-2 둘 다). 대응 사이 일치는 서술. '
                   '등록 r44_gp: 공변량 판(gp_match·gp_any 더함)·제외 판(pibar_noGP) 모두 하한 > 0 이면 "남는다", C-1 점추정이 등록 판의 절반 이하로 줄면 "상당 부분을 나른다".')
    jdump(out, RES + '/r44_corr.json', bk=True); log('R44_CORR', out['robust'], out['sentence']); log('R44_CORR_STATS_DONE')


if __name__ == '__main__':
    for part in PARTS44:
        log('부분 시작', part)
        try: (part_corr() if part == 'corr' else globals()['part_' + part]()); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('R44_CORR_ALL_DONE', ','.join(PARTS44))
