"""(prereg round36 1c, 확인 실험 C) 독립 상과에서 P3 다시(C-1), 그리고 구조가 예측하는 동역학을 넘어서(C-2). r38_p3.py 의 함수를 그대로 쓴다(그 파일은 바꾸지 않음).
부분: pairs → dom → anm → pibar(GPU) → pair → anmpair → stats 사용: python r38_c.py <ROOT> <부분들 쉼표 | all> [작업자 수=3]
  pairs : mdCATH 전체 목록(raw/mdcath_all_domains_listing.json)과 CATH H 로, H-풀 상과 밖·H-풀 도메인 안 쓰는 같은 상과 쌍(1c 의 4,395 쌍·268 상과).
  dom : r38_p3.dom_build 와 같은 앙상블(320 K 뒤 절반, Q ≥ 0.75, 반복마다 ≥ 40 이면 최대 100, 반복 ≥ 2; 같은 seed). 새로 받은 도메인은 raw/mdcath_c/<d>.npz
           (같은 프레임: 반복 0 첫 프레임 + 반복마다 뒤 절반)에서, 이미 있던 도메인은 h5 에서. 메도이드(앙상블 안 d_fix² 합 최소)도 남김.
  anm : H-풀 도메인(메도이드 = derived/r38/p3/dom2)과 C 도메인(메도이드 = 위)마다 ANM 예측 변 방향 퍼짐 윤곽(1c 정의: Cα 15 Å 같은 용수철,
           Hessian 의사역행렬(0 고유값 6 개 뺌), 변 i = tr(P_i Cov(δr_{i+1} − δr_i) P_i) / |b_i|²)과 바름 점검(MD 전체 변 GPA 윤곽과의 Spearman).
  pibar : r38_p3.part_pibar 와 같은 계산(출력만 derived/r38/c; 긴 도메인이라 묶음 크기만 메모리에 맞춤). pair : r38_p3.part_pair 그대로.
  anmpair: 쌍마다 r38_p3.one_pair 와 같은 공유 변(주 = π̄ − 퇴화, 교집합 판)에서 r_ANM, BLOSUM62 평균 점수, 긴 쪽 길이 — H-풀(derived/r38/p3)과 C 둘 다.
  stats : 1c 판정. C-1 = 편 Spearman(r_dyn, 동일성 | 통제), C-2 = 통제에 r_ANM 더함. 구간 = 상과 Exp(1) 군집 부트스트랩 2,000 (쌍 무게 1/상과 쌍 수).
           선다(C-1·C-2 같은 기준) = 주 판 95% 하한 > 0, 그리고 같은 쌍 집합(교집합 ≥ 30 변)의 교집합 판 하한 > 0. 문지기: C-1 이 서면 C-2.
           C-2 는 새 자료(주)와 H-풀 둘에서. 보조: BLOSUM62 판, log 긴 쪽 길이 공변량 판, 상과 크기별, 등록 1b 구간 방법의 C-1.
출력: derived/r38/c/{pairs.json, dom/, medoid/, anm/, pibar.jsonl, pairs_out.jsonl, anm_pairs.jsonl, anm_pairs_hpool.jsonl}, results/r38_c.json"""
import os, sys
ROOT = sys.argv[1]; PARTSC = sys.argv[2].split(',') if len(sys.argv) > 2 else ['all']; NWC = int(sys.argv[3]) if len(sys.argv) > 3 else 3
_src = open(ROOT + '/code/r38_p3.py').read(); _src = _src[:_src.index("ALL = ['pairs'")]
sys.argv = ['r38_p3', ROOT, 'none', str(NWC)]; exec(_src) # r38_p3 의 함수들 (OUT·CACHE 를 아래에서 C 로 바꿈)
P3OUT = OUT; P3CACHE = CACHE
OUT = ROOT + '/derived/r38/c'; CACHE = '/content/ccache'; RAWC = ROOT + '/raw/mdcath_c'
for p_ in (OUT, OUT + '/dom', OUT + '/medoid', OUT + '/anm', CACHE): os.makedirs(p_, exist_ok=True)
ANM_CUT = 15.0; dom_build_h5 = dom_build


# ---------- pairs ----------
def part_pairs():
    import pandas as pd, collections
    alld = json.load(open(ROOT + '/raw/mdcath_all_domains_listing.json'))['domains']; cath = {}
    for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
        if l.startswith('#') or not l.strip(): continue
        f = l.split(); cath[f[0]] = ('.'.join(f[1:5]), int(f[10]))
    hp = pd.read_csv(RES + '/pilot_mdcath_homolog_pairs.tsv', sep='\t'); hsf = set(hp.H); hdom = set(hp.a) | set(hp.b); by = collections.defaultdict(list)
    for d in alld:
        if d in cath and d not in hdom: by[cath[d][0]].append(d)
    prs = sorted([[a, b, H, None, cath[a][1], cath[b][1], None] for H, ds in by.items() if H not in hsf for a, b in itertools.combinations(sorted(ds), 2)], key=lambda p: (p[2], p[0], p[1]))
    doms = sorted({x for p in prs for x in p[:2]})
    jdump(dict(written=time.strftime('%F %T'), n_pairs=len(prs), n_sf=len({p[2] for p in prs}), n_domains=len(doms), domains=doms, pairs=prs, cap50=[],
               rule='1c: mdCATH 전체에서 H-풀 상과 밖, H-풀 도메인 안 쓰는 같은 CATH H 쌍. 열 = a, b, H, (동일성 표 없음), La, Lb, (TM 없음)'), OUT + '/pairs.json', bk=True)
    log('C pairs', len(prs), '쌍', len({p[2] for p in prs}), '상과', len(doms), '도메인')


# ---------- dom ----------
def dom_build_npz(d):
    """r38_p3.dom_build 와 같은 계산; 입력만 raw/mdcath_c/<d>.npz (반복 0 첫 프레임 + 반복마다 뒤 절반)."""
    z = np.load(f'{RAWC}/{d}.npz'); smd = ''.join(z['smd'].tolist()); L = len(smd); first = z['first']
    ca0 = first[:, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
    rng = np.random.default_rng(zlib.crc32(f'r38_p3|{d}'.encode())); runs, lab = [], []
    for r, rk in enumerate(sorted(z['reps'].tolist())): # r = 반복 순번 (h5 판의 enumerate 와 같음)
        Xh = z[f'run{rk}']; idx = np.where(Qf(Xh) >= QT)[0]
        if len(idx) >= NMIN:
            sel = np.sort(rng.choice(idx, min(NMAX, len(idx)), replace=False)); runs.append(Xh[sel]); lab += [r] * len(sel)
    if len(runs) < 2: return None, dict(dom=d, L=L, skip=f'접힌 프레임 ≥ {NMIN} 인 반복 {len(runs)} 개')
    E = frames(np.concatenate(runs)); m = np.asarray(fs.mean_point(E)); A = np.asarray(fs.align_to_reference(E, torch.as_tensor(m)), dtype=np.float32)
    s = 2 * np.linalg.norm(A[:, :, :, 0], axis=-1); sbar = s.mean(0)
    U, S, Vt = np.linalg.svd(m); dt = np.sign(np.linalg.det(U @ Vt)); Dg = np.ones((len(m), 3)); Dg[:, 2] = dt
    Rbar = np.einsum('iab,ib,ibc->iac', U, Dg, Vt); Pm = Rbar * (sbar / 2)[:, None, None]
    degen = S[:, 2] < 0.5 * sbar / 2
    info = dict(dom=d, L=L, ne=L - 1, n_frames=int(len(A)), n_runs=len(runs), seq=smd)
    return dict(A=A, rep=np.array(lab), m=m, Pm=Pm, sbar=sbar, degen=degen, sv=S), info


def dom_build(d):
    return dom_build_npz(d) if os.path.exists(f'{RAWC}/{d}.npz') else dom_build_h5(d)


def _dom_worker_c(d):
    fo = f'{OUT}/dom/{d}.npz'
    if os.path.exists(fo): return d, 'skip'
    if not os.path.exists(f'{RAWC}/{d}.npz') and not os.path.exists(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5'): return d, 'wait'
    try:
        r, info = dom_build(d)
        if r is None: np.savez(fo, info=json.dumps(info, ensure_ascii=False)); return d, info['skip']
        np.savez(f'{CACHE}/{d}.npz', **r); A = r['A']; k = int(np.asarray(fs.fixed_sq_distance(A, A)).sum(1).argmin()); np.save(f'{OUT}/medoid/{d}.npy', A[k])
        np.savez(fo, m=r['m'], Pm=r['Pm'], sbar=r['sbar'], degen=r['degen'], sv=r['sv'], info=json.dumps(info, ensure_ascii=False)); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


def part_dom():
    from multiprocessing import Pool
    doms = json.load(open(OUT + '/pairs.json'))['domains']; cnt = {}
    with Pool(NWC) as p:
        for k, (d, st) in enumerate(p.imap_unordered(_dom_worker_c, doms)):
            cnt[st[:20]] = cnt.get(st[:20], 0) + 1
            if st not in ('skip', 'ok', 'wait') or k % 100 == 0: log('C dom', k, d, st)
    log('C dom 상태', cnt); log('C_DOM_DONE' if cnt.get('wait', 0) == 0 else 'C_DOM_PARTIAL')


# ---------- pibar (GPU) ----------
part_pibar_p3 = part_pibar
def part_pibar(cap=None):
    """r38_p3.part_pibar 와 같은 계산(128 출발, 11 회, EPS, 출력 꼴, 시작 회전 seed = 묶음 첫 순번 b0). 다른 점은 묶음 크기뿐:
    긴 도메인이라 고정 48 대신 묶음마다 Pn·n·m ≤ cap (기본 GPU 메모리 / (128·8·16)), Pn ≤ 48. cap 이 크면 r38_p3 와 같은 묶음이다."""
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'; P = json.load(open(OUT + '/pairs.json'))['pairs']; fo = OUT + '/pibar.jsonl'
    if cap is None: cap = (torch.cuda.get_device_properties(0).total_memory if dev == 'cuda' else 40e9) / (128 * 8 * 16)
    done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    D = {}
    for d in sorted({p[0] for p in P} | {p[1] for p in P}):
        if not os.path.exists(f'{OUT}/dom/{d}.npz'): continue
        s_, _ = dom_small(d)
        if s_ is not None: D[d] = s_['m']
    todo = [p for p in P if (p[0], p[1]) not in done and p[0] in D and p[1] in D]
    todo.sort(key=lambda p: (len(D[p[0]]) // 16, len(D[p[1]]) // 16)); log('pibar', len(todo), '쌍 남음', dev, '묶음 상한', int(cap))
    t0 = time.time(); b0 = 0; nb = 0
    while b0 < len(todo):
        k_ = 1
        while k_ < 48 and b0 + k_ < len(todo) and (k_ + 1) * max(len(D[p[0]]) for p in todo[b0:b0 + k_ + 1]) * max(len(D[p[1]]) for p in todo[b0:b0 + k_ + 1]) <= cap: k_ += 1
        B = todo[b0:b0 + k_]; n = max(len(D[p[0]]) for p in B); m = max(len(D[p[1]]) for p in B)
        x = fs.pad_stack([D[p[0]] for p in B], n).to(dev, torch.float64); y = fs.pad_stack([D[p[1]] for p in B], m).to(dev, torch.float64); Pn = len(B)
        I = torch.eye(3, dtype=x.dtype, device=dev).expand(Pn, 1, 3, 3)
        U = torch.cat([I, random_rotations(Pn * 127, seed=b0, dtype=x.dtype, device=dev).view(Pn, 127, 3, 3)], 1)
        best = torch.full((Pn, 128), -1e30, dtype=x.dtype, device=dev); bG = None
        for it in range(11):
            with torch.no_grad(): K = torch.einsum('pnxz,psmxz->psnm', x, torch.einsum('psxy,pmyz->psmxz', U, y))
            K = K.detach().requires_grad_(True)
            with torch.enable_grad(): Sc = dp_score(K); G, = torch.autograd.grad(Sc.sum(), K)
            with torch.no_grad():
                imp = Sc > best + 1e-7 * (1 + best.abs()); best = torch.where(imp, Sc, best)
                bG = G if bG is None else torch.where(imp[..., None, None], G, bG)
                if it == 10 or not imp.any(): break
                U = torch.where(imp[..., None, None], procrustes(torch.einsum('psnm,pmxz,pnyz->psxy', G, y, x)), U)
        Mt = (bG > 0.5).cpu().numpy(); bs = best.cpu().numpy()
        with open(fo, 'a') as fh:
            for q, p in enumerate(B):
                na, nb_ = len(D[p[0]]), len(D[p[1]]); Smax = float(bs[q].max()); k = int(bs[q].argmax())
                near = np.where(bs[q] >= (1 - EPS) * Smax)[0]; Mq = Mt[q][:, :na, :nb_]
                inter = np.logical_and.reduce(Mq[near], axis=0); bi, bj = np.nonzero(Mq[k]); ii_, jj_ = np.nonzero(inter)
                nd = len({Mq[s_].tobytes() for s_ in near})
                fh.write(json.dumps(dict(ab=[p[0], p[1]], S=Smax, n_near=int(len(near)), n_distinct=nd, best=[[int(a), int(c)] for a, c in zip(bi, bj)], inter=[[int(a), int(c)] for a, c in zip(ii_, jj_)])) + '\n')
        b0 += len(B); nb += 1
        if nb % 20 == 1: log('pibar', b0, '/', len(todo), round(time.time() - t0), 's')
        del K, G, bG, Sc
        if dev == 'cuda': torch.cuda.empty_cache()
    log('PIBAR_DONE')


# ---------- ANM ----------
def anm_edge_spread(ca, cut=ANM_CUT):
    """1c 의 ANM 예측 변 방향 퍼짐 (len(ca) − 1,)."""
    ca = np.asarray(ca, np.float64); N = len(ca); dvec = ca[None] - ca[:, None]; r2 = (dvec ** 2).sum(-1); K = (r2 < cut ** 2) & ~np.eye(N, dtype=bool)
    Hb = np.where(K[..., None, None], -dvec[..., :, None] * dvec[..., None, :] / np.where(K, r2, 1.0)[..., None, None], 0.0)
    Hb[np.arange(N), np.arange(N)] = -Hb.sum(1)
    H = Hb.transpose(0, 2, 1, 3).reshape(3 * N, 3 * N); w, V = np.linalg.eigh(H)
    Cv = (V[:, 6:] / w[6:]) @ V[:, 6:].T; Cb = Cv.reshape(N, 3, N, 3).transpose(0, 2, 1, 3)
    i = np.arange(N - 1); b = ca[i + 1] - ca[i]; bn = b / np.linalg.norm(b, axis=1, keepdims=True)
    cov = Cb[i + 1, i + 1] + Cb[i, i] - Cb[i, i + 1] - Cb[i + 1, i]; Pp = np.eye(3)[None] - bn[:, :, None] * bn[:, None, :]
    return np.einsum('kab,kbc,kca->k', Pp, cov, Pp) / (b ** 2).sum(1), w[:8]


def _anm_worker(args):
    d, src = args; fo = f'{OUT}/anm/{d}.npz'
    if os.path.exists(fo): return d, 'skip'
    try:
        if src == 'hpool':
            f = f'{P3OUT}/dom2/{d}.npz'
            if not os.path.exists(f) or 'medoid' not in np.load(f).files: return d, 'no medoid'
            med = np.load(f)['medoid'].astype(np.float64); Afile = f'{P3CACHE}/{d}.npz'
        else:
            f = f'{OUT}/medoid/{d}.npy'
            if not os.path.exists(f): return d, 'no medoid'
            med = np.load(f).astype(np.float64); Afile = f'{CACHE}/{d}.npz'
        ca = np.asarray(fs.ca_trace(torch.as_tensor(med))); pred, w8 = anm_edge_spread(ca)
        A = None # MD 윤곽(바름 점검): C 는 모두, H-풀은 캐시가 있거나 고정 무작위 200 도메인만(Drive h5 를 다시 다 열지 않음)
        if os.path.exists(Afile): A = np.load(Afile)['A'].astype(np.float64)
        elif src == 'c' or d in HSAN:
            r, _ = (dom_build_h5(d) if src == 'hpool' else dom_build(d)); A = r['A'].astype(np.float64) if r is not None else None
        md = gpa_profile(A, np.arange(A.shape[1])) if A is not None else None; rho = spearman(pred, md) if md is not None else None
        np.savez(fo, pred=pred, md=md if md is not None else np.zeros(0), rho=np.array(np.nan if rho is None else rho), w8=w8); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


HSAN = set()
def part_anm():
    global HSAN
    from multiprocessing import Pool
    hd = json.load(open(P3OUT + '/pairs.json'))['domains']; HSAN = set(np.random.default_rng(20261072).choice(sorted(hd), 200, replace=False).tolist()); cd = json.load(open(OUT + '/pairs.json'))['domains'] if os.path.exists(OUT + '/pairs.json') else []
    jobs = [(d, 'hpool') for d in hd] + [(d, 'c') for d in cd]; cnt = {}
    with Pool(NWC) as p:
        for k, (d, st) in enumerate(p.imap_unordered(_anm_worker, jobs)):
            cnt[st[:20]] = cnt.get(st[:20], 0) + 1
            if st not in ('skip', 'ok', 'no medoid') or k % 200 == 0: log('C anm', k, '/', len(jobs), d, st)
    log('C anm 상태', cnt); log('C_ANM_DONE')


# ---------- anmpair ----------
B62 = None
def blosum():
    global B62
    if B62 is None:
        try:
            from Bio.Align import substitution_matrices
            B62 = substitution_matrices.load('BLOSUM62')
        except Exception:
            B62 = False
    return B62


def anm_pairs(outdir, fo):
    P = json.load(open(outdir + '/pairs.json'))['pairs']; PB = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(outdir + '/pibar.jsonl')}; B = blosum(); sm = {}
    def small(d):
        if d not in sm:
            z = np.load(f'{outdir}/dom/{d}.npz'); info = json.loads(str(z['info'])); a = f'{OUT}/anm/{d}.npz'
            sm[d] = (z['degen'] if 'degen' in z.files else None, info.get('seq'), np.load(a)['pred'] if os.path.exists(a) else None)
        return sm[d]
    with open(fo, 'w') as fh:
        for p in P:
            rec = PB.get((p[0], p[1]))
            if rec is None: continue
            (ga, sa, pa), (gb, sb_, pb) = small(p[0]), small(p[1])
            if ga is None or gb is None: continue
            o = dict(ab=[p[0], p[1]], H=p[2], L_long=max(len(sa), len(sb_)))
            for tag, prs in (('main', rec['best']), ('inter', rec['inter'])):
                prs = np.array(prs, int).reshape(-1, 2); keep = ~(ga[prs[:, 0]] | gb[prs[:, 1]]) if len(prs) else np.zeros(0, bool); sh = prs[keep]
                if len(sh) >= MINSH:
                    o[tag] = dict(n_shared=int(len(sh)), r_anm=spearman(pa[sh[:, 0]], pb[sh[:, 1]]) if pa is not None and pb is not None else None,
                                  blosum=float(np.mean([B[sa[i], sb_[j]] if sa[i] in B.alphabet and sb_[j] in B.alphabet else 0.0 for i, j in sh])) if B is not False else None)
            fh.write(json.dumps(o) + '\n')
    log('anm_pairs', fo)


def part_anmpair():
    anm_pairs(OUT, OUT + '/anm_pairs.jsonl'); anm_pairs(P3OUT, OUT + '/anm_pairs_hpool.jsonl'); log('C_ANMPAIR_DONE')


# ---------- stats ----------
def analyse_sf(rows, key='main', y='r_dyn', x='ident', z=('control',), seed=20261070, B=NB, sf_only=True):
    """r38_p3.analyse 와 같은 점추정, 부트스트랩은 상과 Exp(1) 무게만(1c 의 등록 구간). sf_only=False 면 1b 의 상과 × 두 도메인 곱(비교용)."""
    rows = [r for r in rows if r.get(key) and r[key].get(y) is not None and r[key].get(x) is not None and all(r[key].get(c) is not None for c in z)]
    if len(rows) < 10: return dict(n=len(rows))
    Y = ranks(np.array([r[key][y] for r in rows])); X_ = ranks(np.array([r[key][x] for r in rows])); Zs = [ranks(np.array([r[key][c] for r in rows])) for c in z]
    H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); cnt = np.bincount(si); base = 1.0 / cnt[si]
    doms = sorted({r['ab'][0] for r in rows} | {r['ab'][1] for r in rows}); di = {d: i for i, d in enumerate(doms)}
    ia = np.array([di[r['ab'][0]] for r in rows]); ib = np.array([di[r['ab'][1]] for r in rows])
    est = pcorr(Y, X_, Zs, base); rng = np.random.default_rng(seed); bs = []
    for _ in range(B):
        v = rng.exponential(1, len(sfl)); w = base * v[si]
        if not sf_only: wd = rng.exponential(1, len(doms)); w = w * wd[ia] * wd[ib]
        bs.append(pcorr(Y, X_, Zs, w))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(n=len(rows), n_sf=len(sfl), est=float(est), lo=float(lo), hi=float(hi), se=float(np.std(bs)), mde80=float(2.8 * np.std(bs)))


def merged(outdir, anmf):
    po = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(outdir + '/pairs_out.jsonl')}; an = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(anmf)}; rows = []
    for k, r in po.items():
        a = an.get(k, {}); row = dict(ab=r['ab'], H=r['H'])
        for tag in ('main', 'inter'):
            m = dict(r.get(tag, {})); m.update(a.get(tag, {})); m['logL'] = float(np.log(a.get('L_long', np.nan))) if a.get('L_long') else None; row[tag] = m
        rows.append(row)
    return rows


def judge(rows, z, seed):
    both = [r for r in rows if r['main'].get('r_dyn') is not None and r['inter'].get('r_dyn') is not None and all(r['main'].get(c) is not None for c in z) and all(r['inter'].get(c) is not None for c in z)]
    m = analyse_sf(rows, 'main', z=z, seed=seed); gm = analyse_sf(both, 'main', z=z, seed=seed + 1); gi = analyse_sf(both, 'inter', z=z, seed=seed + 2)
    ok = bool(m.get('lo', -1) > 0 and gi.get('lo', -1) > 0)
    return dict(main=m, gate_main=gm, gate_inter=gi, passed=ok, verdict=verdict(m, gi))


def part_stats():
    rc = merged(OUT, OUT + '/anm_pairs.jsonl'); rh = merged(P3OUT, OUT + '/anm_pairs_hpool.jsonl')
    cdom = set(json.load(open(OUT + '/pairs.json'))['domains']); an = {}; w7 = {}
    for f in glob.glob(OUT + '/anm/*.npz'):
        z_ = np.load(f); x = float(z_['rho']); d = os.path.basename(f)[:-4]; g_ = 'c' if d in cdom else 'hpool'
        w7.setdefault(g_, []).append(float(z_['w8'][6] / z_['w8'][7])) # 7 번째 고유값이 0 에 붙으면(망이 끊김) 의사역행렬이 터진다
        if np.isfinite(x): an.setdefault(g_, []).append(x)
    c1 = judge(rc, ('control',), 20261071); c2 = judge(rc, ('control', 'r_anm'), 20261074); c2h = judge(rh, ('control', 'r_anm'), 20261077)
    sec = dict(blosum=analyse_sf(rc, x='blosum', seed=20261080), logL=analyse_sf(rc, z=('control', 'logL'), seed=20261081),
               small_sf=analyse_sf([r for r in rc if sum(q['H'] == r['H'] for q in rc) <= 2], seed=20261082), large_sf=analyse_sf([r for r in rc if sum(q['H'] == r['H'] for q in rc) > 2], seed=20261083),
               c1_registered_1b_interval=analyse_sf(rc, seed=20261084, sf_only=False), hpool_c1_sf=analyse_sf(rh, seed=20261085), hpool_blosum=analyse_sf(rh, x='blosum', seed=20261086))
    text = ('독립 상과에서 재현되고, 구조가 예측하는 동역학 닮음을 통제해도 남는다' if c1['passed'] and c2['passed'] else
            '재현되지만 구조가 예측하는 동역학으로 설명될 수 있다' if c1['passed'] else f"독립 상과에서 재현되지 않았다(검출 가능한 최소 효과 {c1['main'].get('mde80', float('nan')):.2f}; 긴 도메인 모집단)")
    out = dict(written=time.strftime('%F %T'), C1=c1, C2=c2 if c1['passed'] else dict(c2, note='문지기: C-1 이 서지 않아 판정하지 않음(값은 보고)'), C2_hpool=c2h,
               anm_sanity={k: dict(n=len(v), median_rho=float(np.median(v)), weak=bool(np.median(v) < 0.2), n_anm=len(w7.get(k, [])),
                                   n_mode7_near0=int(sum(x < 1e-6 for x in w7.get(k, [])))) for k, v in an.items()}, secondary=sec, sentence=text,
               rule='prereg 1c: 상과 Exp(1) 군집 부트스트랩 2,000; 선다 = 주 판 하한 > 0 이고 같은 쌍 집합 교집합 판 하한 > 0; 문지기 C-1 → C-2')
    jdump(out, RES + '/r38_c.json', bk=True); log('C', dict(C1=c1['passed'], C2=c2['passed'], C2_hpool=c2h['passed']), text); log('C_STATS_DONE')


ALLC = ['pairs', 'dom', 'anm', 'pibar', 'pair', 'anmpair', 'stats']
if __name__ == '__main__':
    for part in (ALLC if PARTSC == ['all'] else PARTSC):
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('C_DONE', ','.join(PARTSC))
