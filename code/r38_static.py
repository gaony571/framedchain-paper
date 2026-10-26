"""(prereg round36 1 절 정적 축, 범위 = 7g; 2026-10-20 UTC) 그림 4 배경: 메도이드 D.
H-풀 쌍: r38_p3_sens 의 메도이드 정렬(derived/r38/p3/medalign.jsonl, 항등 + 무작위 127 출발)에서 D² = ‖x‖² + ‖y‖² − 2S.
다른 상과 후보 6,000 쌍(results/r37_trunk_candidates.tsv 의 level T·A·C·X 각 1,500): 같은 메도이드 규칙(r38_p3 앙상블 → 앙상블 안 d_fix² 합 최소 프레임;
H-풀 도메인은 dom2 의 메도이드 그대로) · 같은 정렬(fs.align 128 출발, iters 10, 시작 회전 seed = 묶음 첫 순번). 도메인이 길어 묶음 크기만 메모리에 맞춤(r38_c.part_pibar 와 같은 규칙).
BnB(rel_tol 0.005) 인증: 두 묶음을 합친 쌍 목록(정렬 순)에서 무작위 500 (default_rng(20261061)) → 휴리스틱 D 와 S_ub 로 낸 D 하한의 차 = 휴리스틱 오차 상한.
등록은 "모두 BnB" 였으나 7g 에서 이렇게 줄였다(투명성 표).
부분: dom (CPU) → align (GPU) → bnb (GPU) → stats 사용: python r38_static.py <ROOT> <부분들 쉼표> [작업자=3]
출력: derived/r38/static/{medoid/<d>.npy, align.jsonl, bnb.jsonl, table.jsonl}, results/r38_static.json"""
import os, sys
ROOT = sys.argv[1]; PARTSS = sys.argv[2].split(','); NWS = int(sys.argv[3]) if len(sys.argv) > 3 else 3
_src = open(ROOT + '/code/r38_p3.py').read(); _src = _src[:_src.index("ALL = ['pairs'")]
sys.argv = ['r38_p3', ROOT, 'none', str(NWS)]; exec(_src)
import pandas as pd
SOUT = ROOT + '/derived/r38/static'; os.makedirs(SOUT + '/medoid', exist_ok=True)
CAND = pd.read_csv(RES + '/r37_trunk_candidates.tsv', sep='\t'); OTH = CAND[CAND.level.isin(['T', 'A', 'C', 'X'])].sort_values(['level', 'a', 'b']).reset_index(drop=True)
REL = 0.005


def medoid_of(d):
    f2 = f'{OUT}/dom2/{d}.npz'
    if os.path.exists(f2):
        z = np.load(f2)
        if 'medoid' in z.files: return z['medoid'].astype(np.float64)
    f = f'{SOUT}/medoid/{d}.npy'
    return np.load(f).astype(np.float64) if os.path.exists(f) else None


def _dom_worker_s(d):
    fo = f'{SOUT}/medoid/{d}.npy'
    if os.path.exists(fo) or os.path.exists(fo + '.skip'): return d, 'skip'
    f2 = f'{OUT}/dom2/{d}.npz'
    if os.path.exists(f2) and 'medoid' in np.load(f2).files: return d, 'dom2'
    try:
        r, info = dom_build(d)
        if r is None: open(fo + '.skip', 'w').write(info['skip']); return d, 'no ens'
        A = r['A']; k = int(np.asarray(fs.fixed_sq_distance(A, A)).sum(1).argmin()); np.save(fo, A[k]); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


def part_dom():
    from multiprocessing import Pool
    doms = sorted(set(OTH.a) | set(OTH.b)); cnt = {}
    with Pool(NWS) as p:
        for k, (d, st) in enumerate(p.imap_unordered(_dom_worker_s, doms)):
            cnt[st[:20]] = cnt.get(st[:20], 0) + 1
            if st not in ('skip', 'ok', 'dom2', 'no ens') or k % 200 == 0: log('static dom', k, '/', len(doms), d, st)
    log('static dom 상태', cnt); log('STATIC_DOM_DONE')


def part_align(cap=None):
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'; fo = SOUT + '/align.jsonl'; done = set()
    if cap is None: cap = (torch.cuda.get_device_properties(0).total_memory if dev == 'cuda' else 40e9) / (128 * 8 * 16)
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    M = {d: medoid_of(d) for d in sorted(set(OTH.a) | set(OTH.b))}; M = {d: v for d, v in M.items() if v is not None}
    P = [[r.a, r.b, r.level, float(r.TM)] for r in OTH.itertuples()]
    todo = sorted([p for p in P if (p[0], p[1]) not in done and p[0] in M and p[1] in M], key=lambda p: (len(M[p[0]]) // 16, len(M[p[1]]) // 16)); log('static align', len(todo), '쌍', dev, '묶음 상한', int(cap))
    b0 = 0; nb = 0; t0 = time.time()
    while b0 < len(todo):
        k_ = 1
        while k_ < 48 and b0 + k_ < len(todo) and (k_ + 1) * max(len(M[p[0]]) for p in todo[b0:b0 + k_ + 1]) * max(len(M[p[1]]) for p in todo[b0:b0 + k_ + 1]) <= cap: k_ += 1
        B = todo[b0:b0 + k_]; n = max(len(M[p[0]]) for p in B); m = max(len(M[p[1]]) for p in B)
        x = fs.pad_stack([M[p[0]] for p in B], n).to(dev, torch.float64); y = fs.pad_stack([M[p[1]] for p in B], m).to(dev, torch.float64)
        al = fs.align(x, y, n_random=127, iters=10, seed=b0); Ms = al.match.cpu().numpy(); Us = al.rotation.cpu().numpy(); Sv = al.score.cpu().numpy()
        with open(fo, 'a') as fh:
            for q, p in enumerate(B):
                na, nb_ = len(M[p[0]]), len(M[p[1]]); ii, jj = np.nonzero(Ms[q][:na, :nb_]); a_, b_ = M[p[0]][ii], M[p[1]][jj]
                share = float(((a_ - np.einsum('xy,iyz->ixz', Us[q], b_)) ** 2).sum() / max(len(ii), 1))
                fh.write(json.dumps(dict(ab=[p[0], p[1]], level=p[2], TM=p[3], S=float(Sv[q]), n_match=int(len(ii)), share=share)) + '\n')
        b0 += len(B); nb += 1
        if nb % 20 == 1: log('static align', b0, '/', len(todo), round(time.time() - t0), 's')
        del al, x, y
        if dev == 'cuda': torch.cuda.empty_cache()
    log('STATIC_ALIGN_DONE')


def union():
    H = [json.loads(l) for l in open(OUT + '/medalign.jsonl')]; O = [json.loads(l) for l in open(SOUT + '/align.jsonl')]
    S = {('H', r['ab'][0], r['ab'][1]): r['S'] for r in H}; S.update({('O', r['ab'][0], r['ab'][1]): r['S'] for r in O}); return S, H, O


def part_bnb():
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'; fo = SOUT + '/bnb.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: r = json.loads(l); done.add((r['src'], r['ab'][0], r['ab'][1]))
            except Exception: pass
    S, _, _ = union(); keys = sorted(S); pick = [keys[i] for i in sorted(np.random.default_rng(20261061).choice(len(keys), 500, replace=False))]
    log('static bnb', len(pick), '쌍 (합친 목록', len(keys), ')', len([k for k in pick if k not in done]), '남음', dev)
    for k, key in enumerate(pick):
        if key in done: continue
        src, a, b = key; x, y = medoid_of(a), medoid_of(b)
        try:
            c = fs.bnb(x, y, rel_tol=REL, device=dev); nx, ny = float((x ** 2).sum()), float((y ** 2).sum())
            Dh = float(np.sqrt(max(nx + ny - 2 * S[key], 0))); Dlo = float(np.sqrt(max(nx + ny - 2 * c['S_ub'], 0)))
            o = dict(src=src, ab=[a, b], n=len(x), m=len(y), S_heur=S[key], S_bnb_heur=float(c['S_heur']), S_lb=float(c['S_lb']), S_ub=float(c['S_ub']), certified=bool(c['certified']),
                     D_heur=Dh, D_lo=Dlo, err_ub=Dh - Dlo, rel_err_ub=(Dh - Dlo) / Dh if Dh > 0 else None, sec=float(c['sec']))
        except Exception as e:
            o = dict(src=src, ab=[a, b], err=repr(e)[:200])
        with open(fo, 'a') as fh: fh.write(json.dumps(o) + '\n')
        if k % 25 == 0: log('static bnb', k, '/', len(pick), a, b, o.get('certified'), o.get('rel_err_ub'))
        if dev == 'cuda': torch.cuda.empty_cache()
    log('STATIC_BNB_DONE')


def part_stats():
    S, H, O = union(); Pj = json.load(open(OUT + '/pairs.json')); tmH = {(p[0], p[1]): (p[6], p[3]) for p in Pj['pairs']}; nrm = {}
    def nr(d):
        if d not in nrm: x = medoid_of(d); nrm[d] = (float((x ** 2).sum()), len(x))
        return nrm[d]
    tab = []
    for src, rows in (('H', H), ('O', O)):
        for r in rows:
            a, b = r['ab']; (xa, na), (xb, nb_) = nr(a), nr(b); D2 = max(xa + xb - 2 * r['S'], 0.0); nbar = (na + nb_) / 2
            tm, idn = (tmH.get((a, b), (None, None)) if src == 'H' else (r['TM'], None))
            tab.append(dict(src=src, level='H' if src == 'H' else r['level'], ab=[a, b], TM=tm, ident=idn, ne_a=na, ne_b=nb_, D2=D2, D=float(np.sqrt(D2)), D2n=D2 / np.sqrt(nbar), share=r['share'], n_match=r['n_match']))
    with open(SOUT + '/table.jsonl', 'w') as fh:
        for t in tab: fh.write(json.dumps(t) + '\n')
    bins = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.01]; byb = {}
    for lv in ('H', 'T', 'A', 'C', 'X'):
        for lo, hi in zip(bins[:-1], bins[1:]):
            v = [t['D2n'] for t in tab if t['level'] == lv and t['TM'] is not None and lo <= t['TM'] < hi]
            if v: byb[f'{lv}|{lo}-{hi}'] = dict(n=len(v), median_D2n=float(np.median(v)), q25=float(np.percentile(v, 25)), q75=float(np.percentile(v, 75)))
    bn = [json.loads(l) for l in open(SOUT + '/bnb.jsonl')] if os.path.exists(SOUT + '/bnb.jsonl') else []; ok = [r for r in bn if 'err' not in r]
    re_ = np.array([r['rel_err_ub'] for r in ok if r['rel_err_ub'] is not None])
    bsum = dict(n=len(bn), n_err=len(bn) - len(ok), n_certified=sum(r['certified'] for r in ok), rel_err_ub_median=float(np.median(re_)) if len(re_) else None,
                rel_err_ub_p95=float(np.percentile(re_, 95)) if len(re_) else None, rel_err_ub_max=float(re_.max()) if len(re_) else None,
                err_ub_max=float(max(r['err_ub'] for r in ok)) if ok else None, n_heur_below_lb=sum(r['S_heur'] < r['S_lb'] - 1e-9 for r in ok),
                by_src={s: dict(n=sum(r['src'] == s for r in ok), rel_err_ub_max=float(max([r['rel_err_ub'] or 0 for r in ok if r['src'] == s], default=0))) for s in ('H', 'O')})
    out = dict(written=time.strftime('%F %T'), n_H=len(H), n_O=len(O), n_O_registered=len(OTH), by_level_TM=byb, bnb=bsum,
               rule='prereg 1 절 정적 축(범위 7g): 메도이드 D² = ‖x‖² + ‖y‖² − 2S (128 출발), D2n = D²/√n̄. BnB 무작위 500 (seed 20261061, rel_tol 0.005) 으로 휴리스틱 오차 상한(D_heur − D_lo, D_lo = S_ub 로 낸 D).')
    jdump(out, RES + '/r38_static.json', bk=True); log('STATIC', len(H), len(O), 'BnB', bsum); log('STATIC_STATS_DONE')


if __name__ == '__main__':
    for part in PARTSS:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
