"""(prereg round36 7 절, 1차 P1 — 트랙 B5) 시연 1: 열 적응 θ*. r38_p3.py 의 함수(md_atoms·frames·gpa_profile·Q 정의)를 불러 쓴다(그 파일은 바꾸지 않음).
부분: pairs → mean → align → spread → sim → pre → stats      사용: python r38_p1.py <ROOT> <부분들 쉼표 | all> [작업자 수=3]
  pairs : T-풀 = H-풀 쌍(derived/r38/p3/pairs.json) 가운데 HT = 호열(OGT ≥ 60 °C)–원핵 중온(OGT < 40 °C, prok) 쌍,
          MM = HT 가 있는 상과 안의 원핵 중온–원핵 중온 쌍(대칭 판), PRE = HT 가 없는 상과의 원핵 중온–원핵 중온 쌍(두 규약 차의 사전 측정, 기록만).
          OGT = results/pilot_mdcath_organism_temp_v2.tsv (Topt, prok). MM·PRE 의 데워지는 쪽 = (a, b) 정렬 목록 순서로 default_rng(SEED_DIR).integers(0, 2)
          (0 → a). 두 종류 쌍이 모두 있는 상과 수를 센다(20 미만이면 짝짓지 않은 판이 주).
  mean  : 도메인마다 320 K·뒤 절반·Q ≥ 0.75 칸의 평균 점(π̄ 용; spread 의 같은 칸과 같은 프레임).
  align : 쌍마다 두 320 K 평균 점 사이 128 출발(항등 + 무작위 127, 쌍 이름 seed; 번갈아 10 번) 정렬. π̄ = 최적 출발의 짝,
          alts = 거의 최적(S ≥ (1 − 0.005) S_최적) 가운데 π̄ 와 다른 대응들(S 큰 순 최대 7) — θ* 흔들림(0.4) 용.
  spread: 도메인마다 h5 를 한 번 열어 칸 (T, 창, Q) 을 만든다: T ∈ 320·348·379·413, (창, Q) ∈ (뒤 절반, 0.75) 주 판, (뒤 절반, 0.7)·(뒤 절반, 0.8)·
          (50–100 ns, 0.75) 민감도. 반복마다 그 창의 접힌 프레임(첫 320 K 프레임의 천연 접촉; r38_p3 와 같은 Q) ≥ N 이면 N 개 무작위
          (뒤 절반 N = 100, 50–100 ns N = 50; 도메인·T·창·Q 이름 seed), 그런 반복 ≥ 2 면 그 칸을 쓴다. 칸의 프레임을 전체 변 평균 점에 맞춘 뒤,
          그 도메인이 든 쌍마다 공유 변 총퍼짐(Σ 변별 퍼짐): gpa = 공유 변 GPA(자기 평균 둘레, 3 번; 주 판), full = 전체 변 중첩 그대로의 공유 변 퍼짐(민감도),
          runs = 반복마다 gpa(시뮬레이션 잡음). 주 칸에서는 alts 대응마다 gpa 도. 반복별 궤적 길이(프레임 수)도 남긴다.
  sim   : 추정량 시뮬레이션(판정 전): 실제 쌍의 상과 구조·쓸 수 있는 온도·u 와 반복 잡음(log 퍼짐의 반복 사이 SE)으로
          LR(T) = b (T − θ) + ε_h(T) − ε_r 를 만들고(b = 0.171 / 28 K, 결과 3.2), θ ∈ 320·335·350·380·420 K (모든 쌍 같음, 그리고 쌍마다 N(0, 15 K) 흩어짐)
          에서 θ* → 상과 요약 → 상과 중앙값 → 부트스트랩 2,000 을 그대로 돌려(반복 200) 치우침과 95% 구간 피복률. 피복률 < 0.90 이면 구간 방법을 고친 뒤 동결.
  pre   : PRE 쌍에서 두 규약(gpa / full)의 θ* 차 — 기록(판정·규칙 조정에 쓰지 않음).
  stats : 주 판 조건 1·2, 검열 비율, 민감도(전체 변 중첩, Q 0.7·0.8, 50–100 ns, 413 K 뺌, 성향 점수 재가중), θ* 흔들림, 궤적 길이, 비단조, 2차 용량–반응.
θ*: LR(T) = log(S_데워진 쪽(T) / S_기준 쪽(320 K)). 데워진 쪽의 쓸 수 있는 온도(320 K 포함 둘 이상; 320 K 칸이 없으면 뺌, 수 기록) 사이 선형 보간으로
   LR = 0 인 온도. LR(320) > 0 → "< 320". 처음으로 LR 이 0 이상이 되는 가장 낮은 구간에서 보간(비단조 표시). u(데워진 쪽 맨 위 쓸 수 있는 온도)까지
   못 넘으면 "> u". 기준 쪽 320 K 칸이 없으면 뺌. 공유 변 < 30 이면 뺌.
집계: 쌍 → 상과 = 아래 중앙값(호열–중온 쪽: "> u" 는 u 자리, 같은 값이면 위; 대칭 판 쪽: "> u" 는 가장 큰 값, 중앙값이 검열 쌍이면 "> u_C",
   u_C = 하한 순위 중앙값; "< 320" 은 양쪽 모두 가장 작음). 조건 1 = 호열–중온 상과 요약의 상과 중앙값 95% 하한 > 320 (상과 부트스트랩 2,000;
   상과 요약의 왼쪽 검열 > 50% 면 틀림, 오른쪽 검열 > 50% 면 맞음). 조건 2 = 짝지은 상과 차 d_s (prereg 표, 부호 미정·같은 방향 이중 검열 상과는 뺌)의
   상과 중앙값 95% 하한 > 0; 뺀 상과가 절반을 넘거나 두 종류가 다 있는 상과 < 20 이면 짝짓지 않은 판. P1 = 조건 1 과 조건 2 가 모두 섬.
출력: derived/r38/p1/{pairs.json, mean/<d>.npz, align.jsonl, spread/<d>.json}, results/r38_p1_{sim,pre}.json, results/r38_p1.json (사본 results/_backup/<날짜>/)"""
import os, sys
ROOT = sys.argv[1]; PARTS1 = sys.argv[2].split(',') if len(sys.argv) > 2 else ['all']; NW1 = int(sys.argv[3]) if len(sys.argv) > 3 else 3
_src = open(ROOT + '/code/r38_p3.py').read(); _src = _src[:_src.index("ALL = ['pairs'")]
sys.argv = ['r38_p3', ROOT, 'none', str(NW1)]; exec(_src)                 # md_atoms, frames, gpa_profile, log, jdump, backup, fs, torch, np, h5py, zlib ...
OUT1 = ROOT + '/derived/r38/p1'
for p_ in (OUT1, OUT1 + '/mean', OUT1 + '/spread'): os.makedirs(p_, exist_ok=True)
TT = ('320', '348', '379', '413'); TK = {t: float(t) for t in TT}
CELLS = [('half', 0.75), ('half', 0.7), ('half', 0.8), ('abs', 0.75)]; MAIN = ('half', 0.75)
NWIN = {'half': 100, 'abs': 50}; MINRUN = 2; MINSH1 = 30; SEED_DIR = 20261040; NALT = 7; EPS1 = 0.005
B_SLOPE = 0.171 / 28.0; SIM_THETA = (320.0, 335.0, 350.0, 380.0, 420.0); SIM_REP = 200; NB1 = 2000
ckey = lambda T, w, q: f'{T}|{w}|{q:.2f}'


# ---------- pairs ----------
def part_pairs():
    import pandas as pd
    org = pd.read_csv(RES + '/pilot_mdcath_organism_temp_v2.tsv', sep='\t'); om = {r.dom: r for r in org.itertuples()}
    def kind(d):
        r = om.get(d)
        if r is None or not np.isfinite(r.Topt): return None
        if r.Topt >= 60: return 'T'
        if r.Topt < 40 and bool(r.prok): return 'M'
        return None
    pj = json.load(open(OUT + '/pairs.json')); P = pj['pairs']; HT, MM = [], []
    for p in P:
        ka, kb = kind(p[0]), kind(p[1])
        if {ka, kb} == {'T', 'M'}:
            th, me = (p[0], p[1]) if ka == 'T' else (p[1], p[0])
            HT.append(dict(a=p[0], b=p[1], H=p[2], ident=p[3], La=p[4], Lb=p[5], tm=p[6], hot=th, ref=me, dOGT=float(om[th].Topt - om[me].Topt)))
        elif ka == 'M' and kb == 'M': MM.append(dict(a=p[0], b=p[1], H=p[2], ident=p[3], La=p[4], Lb=p[5], tm=p[6]))
    sfT = {r['H'] for r in HT}; mm_in = sorted([r for r in MM if r['H'] in sfT], key=lambda r: (r['a'], r['b'])); pre = sorted([r for r in MM if r['H'] not in sfT], key=lambda r: (r['a'], r['b']))
    rng = np.random.default_rng(SEED_DIR)
    for r in mm_in + pre:
        o = int(rng.integers(0, 2)); r['hot'], r['ref'] = (r['a'], r['b']) if o == 0 else (r['b'], r['a'])
    both = sorted(sfT & {r['H'] for r in mm_in})
    doms = sorted({r[k] for r in HT + mm_in + pre for k in ('a', 'b')})
    jdump(dict(written=time.strftime('%F %T'), HT=HT, MM=mm_in, PRE=pre, domains=doms, n=dict(HT=len(HT), HT_sf=len(sfT), MM=len(mm_in), MM_sf=len({r['H'] for r in mm_in}),
               sf_both=len(both), PRE=len(pre), PRE_sf=len({r['H'] for r in pre}), domains=len(doms)), sf_both=both, seed_dir=SEED_DIR,
               rule='HT = 호열(Topt ≥ 60)–원핵 중온(Topt < 40, prok); MM = HT 상과 안 원핵 중온 쌍; PRE = 나머지 상과의 원핵 중온 쌍. 방향 = (a, b) 정렬 순서로 default_rng(20261040).integers(0, 2)'),
          OUT1 + '/pairs.json', bk=True)
    log('P1 pairs', dict(HT=len(HT), HT_sf=len(sfT), MM=len(mm_in), sf_both=len(both), PRE=len(pre), doms=len(doms)))


def all_pairs():
    pj = json.load(open(OUT1 + '/pairs.json')); return [dict(r, kind=k) for k in ('HT', 'MM', 'PRE') for r in pj[k]]


# ---------- 칸 만들기 ----------
def build_cells(d, temps=TT, cells=CELLS):
    """{ckey: (frames A (n, ne, 3, 3) 평균에 맞춤, run 라벨)} 와 반복별 궤적 길이 {T: [..]}."""
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode()); L = len(smd)
    ld = lambda T: [g[T][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g[T].keys(), key=int)]
    ref = ld('320'); ca0 = ref[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
    out, lens = {}, {}
    for T in temps:
        if T not in g: continue
        reps = ref if T == '320' else ld(T); lens[T] = [int(len(X)) for X in reps]; qv = [Qf(X) for X in reps]
        for w, q in cells:
            N = NWIN[w]; rng = np.random.default_rng(zlib.crc32(f'r38_p1|{d}|{T}|{w}|{q:.2f}'.encode())); runs, lab = [], []
            for r, X in enumerate(reps):
                lo, hi = (len(X) // 2, len(X)) if w == 'half' else (50, 100)
                if hi > len(X): continue
                idx = lo + np.where(qv[r][lo:hi] >= q)[0]
                if len(idx) >= N: runs.append(X[np.sort(rng.choice(idx, N, replace=False))]); lab += [r] * N
            if len(runs) >= MINRUN:
                E = frames(np.concatenate(runs)); m = np.asarray(fs.mean_point(E)); A = np.asarray(fs.align_to_reference(E, torch.as_tensor(m)), dtype=np.float64)
                out[ckey(T, w, q)] = (A, np.array(lab), m)
    h.close(); return out, lens, smd


# ---------- mean ----------
def _mean_worker(d):
    fo = f'{OUT1}/mean/{d}.npz'
    if os.path.exists(fo): return d, 'skip'
    try:
        c, lens, smd = build_cells(d, temps=('320',), cells=[MAIN])
        k = ckey('320', *MAIN)
        if k not in c: np.savez(fo, none=np.array(1)); return d, 'none'
        np.savez(fo, m=c[k][2], n=np.array(len(c[k][0]))); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


def part_mean():
    from multiprocessing import Pool
    doms = json.load(open(OUT1 + '/pairs.json'))['domains']
    with Pool(NW1) as p:
        for k, (d, s_) in enumerate(p.imap_unordered(_mean_worker, doms)):
            if s_ not in ('skip', 'ok') or k % 100 == 0: log('p1 mean', k, d, s_)
    log('P1_MEAN_DONE')


def mean_of(d):
    f = f'{OUT1}/mean/{d}.npz'
    if not os.path.exists(f): return None
    z = np.load(f); return z['m'] if 'm' in z.files else None


# ---------- align ----------
def align_batch(B, Ms, dev='cpu'):
    """P3 part_pibar 와 같은 계산(쌍마다 128 출발). 반환: [(S, best 짝 [[i, j]], alts [[[i, j]...]...], n_near)]"""
    n = max(len(Ms[p[0]]) for p in B); m = max(len(Ms[p[1]]) for p in B); Pn = len(B)
    x = fs.pad_stack([Ms[p[0]] for p in B], n).to(dev, torch.float64); y = fs.pad_stack([Ms[p[1]] for p in B], m).to(dev, torch.float64)
    I = torch.eye(3, dtype=x.dtype, device=dev).expand(Pn, 1, 3, 3)
    Rr = torch.stack([random_rotations(127, seed=zlib.crc32(f'r38_p1|{p[0]}|{p[1]}'.encode()) % (2 ** 31), dtype=x.dtype, device=dev) for p in B])
    U = torch.cat([I, Rr], 1); best = torch.full((Pn, 128), -1e30, dtype=x.dtype, device=dev); bG = None
    for it in range(11):
        with torch.no_grad(): K = torch.einsum('pnxz,psmxz->psnm', x, torch.einsum('psxy,pmyz->psmxz', U, y))
        K = K.detach().requires_grad_(True)
        with torch.enable_grad(): Sc = dp_score(K); G, = torch.autograd.grad(Sc.sum(), K)
        with torch.no_grad():
            imp = Sc > best + 1e-7 * (1 + best.abs()); best = torch.where(imp, Sc, best)
            bG = G if bG is None else torch.where(imp[..., None, None], G, bG)
            if it == 10 or not imp.any(): break
            U = torch.where(imp[..., None, None], procrustes(torch.einsum('psnm,pmxz,pnyz->psxy', G, y, x)), U)
    Mt = (bG > 0.5).cpu().numpy(); bs = best.cpu().numpy(); res = []
    for q, p in enumerate(B):
        na, nb = len(Ms[p[0]]), len(Ms[p[1]]); Smax = float(bs[q].max()); k = int(bs[q].argmax()); Mq = Mt[q][:, :na, :nb]
        near = [s_ for s_ in np.argsort(-bs[q]) if bs[q][s_] >= (1 - EPS1) * Smax]; seen = {Mq[k].tobytes()}; alts = []
        for s_ in near:
            b_ = Mq[s_].tobytes()
            if b_ not in seen and len(alts) < NALT: seen.add(b_); alts.append(np.argwhere(Mq[s_]).tolist())
        res.append((Smax, np.argwhere(Mq[k]).tolist(), alts, len(near)))
    return res


def part_align():
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'; fo = OUT1 + '/align.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    Ps = {(r['a'], r['b']) for r in all_pairs()}; Ms = {}
    for d in sorted({x for p in Ps for x in p}):
        m = mean_of(d)
        if m is not None: Ms[d] = m
    todo = sorted([p for p in Ps if p not in done and p[0] in Ms and p[1] in Ms], key=lambda p: (len(Ms[p[0]]) // 16, len(Ms[p[1]]) // 16)); log('p1 align', len(todo), dev)
    BS = 48 if dev == 'cuda' else 8; t0 = time.time()
    for b0 in range(0, len(todo), BS):
        B = todo[b0:b0 + BS]; res = align_batch(B, Ms, dev)
        with open(fo, 'a') as fh:
            for p, (S, best, alts, nn) in zip(B, res): fh.write(json.dumps(dict(ab=list(p), S=S, best=best, alts=alts, n_near=nn)) + '\n')
        if (b0 // BS) % 25 == 0: log('p1 align', b0 + len(B), '/', len(todo), round(time.time() - t0), 's')
    log('P1_ALIGN_DONE')


# ---------- spread ----------
def full_spread(A, idx):
    X = A[:, idx]; return ((X - X.mean(0)) ** 2).sum((-1, -2)).mean(0)


def _spread_worker(d):
    fo = f'{OUT1}/spread/{d}.json'
    if os.path.exists(fo): return d, 'skip'
    try:
        mine = [(k, v) for k, v in AL.items() if d in k]
        if not mine: json.dump(dict(dom=d, skip='쌍 없음'), open(fo, 'w')); return d, 'none'
        c, lens, smd = build_cells(d); res = dict(dom=d, lens=lens, cells=sorted(c), pairs={})
        for (a, b), al in mine:
            side = 0 if d == a else 1; eds = lambda prs: np.array(sorted({x[side] for x in prs}), int)
            sets = [eds(al['best'])] + [eds(x) for x in al['alts']]; o = {}
            for kk, (A, lab, m) in c.items():
                e0 = sets[0]
                if len(e0) < MINSH1: continue
                g0 = gpa_profile(A, e0); rr = [float(gpa_profile(A[lab == r], e0).sum()) for r in sorted(set(lab.tolist()))]
                o[kk] = dict(gpa=float(g0.sum()), full=float(full_spread(A, e0).sum()), runs=rr)
                if kk.endswith('|' + MAIN[0] + f'|{MAIN[1]:.2f}'): o[kk]['alts'] = [float(gpa_profile(A, e).sum()) if len(e) >= MINSH1 else None for e in sets[1:]]
            res['pairs'][f'{a}__{b}'] = dict(n_shared=int(len(sets[0])), cells=o)
        json.dump(res, open(fo, 'w')); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


AL = {}
def part_spread():
    global AL
    from multiprocessing import Pool
    AL = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(OUT1 + '/align.jsonl')}
    doms = sorted({x for p in AL for x in p})
    with Pool(NW1) as p:
        for k, (d, s_) in enumerate(p.imap_unordered(_spread_worker, doms)):
            if s_ not in ('skip', 'ok') or k % 50 == 0: log('p1 spread', k, '/', len(doms), d, s_)
    log('P1_SPREAD_DONE')


# ---------- θ* 와 집계 ----------
def theta(LR, temps):
    """LR: 온도 순서대로 log 퍼짐 비, temps: 그 온도(K). 반환 (종류, 값, 비단조): ('num', θ) | ('gt', u) | ('lt', 320)."""
    LR = np.asarray(LR, float); x = np.asarray(temps, float); nonmono = bool(np.any(np.diff(LR) < 0))
    if LR[0] > 0: return ('lt', 320.0, nonmono)
    if LR[0] == 0: return ('num', float(x[0]), nonmono)
    for k in range(len(LR) - 1):
        if LR[k] < 0 <= LR[k + 1]: return ('num', float(x[k] + (0 - LR[k]) * (x[k + 1] - x[k]) / (LR[k + 1] - LR[k])), nonmono)
    return ('gt', float(x[-1]), nonmono)


def key_ht(it): return (-np.inf, 0) if it[0] == 'lt' else ((it[1], 1) if it[0] == 'gt' else (it[1], 0))
def key_mm(it): return (-np.inf, 0) if it[0] == 'lt' else ((np.inf, 0) if it[0] == 'gt' else (it[1], 0))


def wmedian(items, key, w=None):
    """아래 중앙값(무게가 있으면 누적 무게가 절반 이상이 되는 첫 항목; 같은 무게면 아래 중앙값과 같음)."""
    o = sorted(range(len(items)), key=lambda i: key(items[i])); w = np.ones(len(items)) if w is None else np.asarray(w, float)
    cw = np.cumsum(w[o]); return items[o[int(np.searchsorted(cw, cw[-1] / 2 - 1e-12))]]


def sf_summary(items, side, w=None):
    if side == 'ht': return wmedian(items, key_ht, w)
    med = wmedian(items, key_mm, w)
    if med[0] == 'gt': lo = wmedian(items, key_ht, w); return ('gt', float(lo[1]))
    return med


def dsum(T, C):
    """prereg 7 절 조건 2 표. 반환 ('num', d) | ('le0', 0.0) | ('excl', 까닭)."""
    tk, tv = T[0], T[1]; ck, cv = C[0], C[1]
    if tk == ck and tk in ('gt', 'lt'): return ('excl', '같은 방향 이중 검열')
    if tk in ('num', 'gt') and ck == 'num': return ('num', tv - cv)
    if tk in ('num', 'gt') and ck == 'lt': return ('num', tv - 320.0)
    if tk == 'lt': return ('le0', 0.0)
    if tk == 'num' and ck == 'gt': return ('le0', 0.0) if tv <= cv else ('excl', '부호 미정')
    raise ValueError((T, C))


def boot_median(items, key, rng, B=NB1, w=None):
    n = len(items); meds = []
    for _ in range(B):
        ix = rng.integers(0, n, n); meds.append(wmedian([items[i] for i in ix], key, None if w is None else [w[i] for i in ix]))
    s = sorted(meds, key=key); return s[int(np.floor(0.025 * B))], s[int(np.ceil(0.975 * B)) - 1]


def cond1(T_s, rng, B=NB1):
    items = list(T_s.values()); med = wmedian(items, key_ht); lo, hi = boot_median(items, key_ht, rng, B)
    fl = float(np.mean([i[0] == 'lt' for i in items])); fg = float(np.mean([i[0] == 'gt' for i in items]))
    ok = bool(key_ht(lo) > (320.0, 0)) if fl <= 0.5 and fg <= 0.5 else (fg > 0.5)
    return dict(n_sf=len(items), median=list(med), ci95=[list(lo), list(hi)], frac_left_cens=fl, frac_right_cens=fg, passed=ok,
                rule='하한 순위 상과 중앙값의 상과 부트스트랩 95% 하한 > 320; 왼쪽 검열 > 50% 면 틀림, 오른쪽 검열 > 50% 면 맞음')


def cond2(T_s, C_s, rng, B=NB1):
    both = sorted(set(T_s) & set(C_s)); ds = {h: dsum(T_s[h], C_s[h]) for h in both}; excl = {h: v[1] for h, v in ds.items() if v[0] == 'excl'}
    keep = [ds[h] for h in both if ds[h][0] != 'excl']; kd = lambda it: (it[1], 0)
    out = dict(n_sf_both=len(both), n_excl=len(excl), excl=excl, n_le0=sum(v[0] == 'le0' for v in keep))
    if len(both) < 20 or len(excl) > len(both) / 2 or not keep:
        Tm, Cm = wmedian(list(T_s.values()), key_ht), sf_summary(list(C_s.values()), 'mm'); d0 = dsum(Tm, Cm); bs = []
        tl, cl = list(T_s.values()), list(C_s.values())
        for _ in range(B):
            a_ = [tl[i] for i in rng.integers(0, len(tl), len(tl))]; c_ = [cl[i] for i in rng.integers(0, len(cl), len(cl))]
            v = dsum(wmedian(a_, key_ht), sf_summary(c_, 'mm'))
            if v[0] != 'excl': bs.append(v[1])
        lo, hi = (np.percentile(bs, [2.5, 97.5]) if bs else (np.nan, np.nan))
        out.update(version='짝짓지 않은 판', d=list(d0), ci95=[float(lo), float(hi)], passed=bool(lo > 0))
    else:
        med = wmedian(keep, kd); lo, hi = boot_median(keep, kd, rng, B)
        out.update(version='짝지은 판', d=list(med), ci95=[list(lo), list(hi)], passed=bool(lo[1] > 0 and lo[0] == 'num'))
    return out


def pair_theta(sp, r, cell, conv='gpa', drop413=False, alt=None):
    """sp: {dom: spread json}, r: 쌍 dict. 반환 (종류, 값, 비단조) 또는 ('skip', 까닭)."""
    key = f"{r['a']}__{r['b']}"; H_, R_ = sp.get(r['hot']), sp.get(r['ref'])
    if H_ is None or R_ is None or key not in H_['pairs'] or key not in R_['pairs']: return ('skip', '퍼짐 없음')
    if H_['pairs'][key]['n_shared'] < MINSH1: return ('skip', '공유 변 < 30')
    hc, rc = H_['pairs'][key]['cells'], R_['pairs'][key]['cells']; w, q = cell; k320 = ckey('320', w, q)
    if k320 not in rc: return ('skip', '기준 쪽 320 K 칸 없음')
    if k320 not in hc: return ('skip', '데워진 쪽 320 K 칸 없음')
    val = lambda c: (c['gpa'] if alt is None else c['alts'][alt]) if conv == 'gpa' else c['full']
    Sr = val(rc[k320]); temps, LR = [], []
    for T in TT:
        if drop413 and T == '413': continue
        kk = ckey(T, w, q)
        if kk in hc and val(hc[kk]) is not None and Sr is not None: temps.append(TK[T]); LR.append(np.log(val(hc[kk]) / Sr))
    if not temps or temps[0] != 320.0: return ('skip', '데워진 쪽 320 K 값 없음')
    if len(temps) < 2: return ('skip', '데워진 쪽 온도 점 < 2')
    return theta(LR, temps)


def summarize(pairs, thetas, side, w=None):
    by = {}
    for r, t in zip(pairs, thetas):
        if t[0] == 'skip': continue
        by.setdefault(r['H'], []).append((t[:2], 1.0 if w is None else w[(r['a'], r['b'])]))
    return {h: sf_summary([x[0] for x in v], side, [x[1] for x in v]) for h, v in by.items()}


def verdict_set(HT, MM, sp, cell=MAIN, conv='gpa', drop413=False, wmm=None, seed=20261041):
    rng = np.random.default_rng(seed)
    tH = [pair_theta(sp, r, cell, conv, drop413) for r in HT]; tM = [pair_theta(sp, r, cell, conv, drop413) for r in MM]
    T_s, C_s = summarize(HT, tH, 'ht'), summarize(MM, tM, 'mm', wmm)
    c1, c2 = cond1(T_s, rng), cond2(T_s, C_s, rng)
    cnt = lambda ts: dict(n=len(ts), used=sum(t[0] != 'skip' for t in ts), left=sum(t[0] == 'lt' for t in ts), right=sum(t[0] == 'gt' for t in ts),
                          nonmono=sum(1 for t in ts if t[0] != 'skip' and t[2]), skip={k: sum(1 for t in ts if t[0] == 'skip' and t[1] == k) for k in {t[1] for t in ts if t[0] == 'skip'}})
    return dict(cond1=c1, cond2=c2, P1=bool(c1['passed'] and c2['passed']), pairs_HT=cnt(tH), pairs_MM=cnt(tM)), tH, tM


def load_spread(doms):
    sp = {}
    for d in doms:
        f = f'{OUT1}/spread/{d}.json'
        if os.path.exists(f):
            z = json.load(open(f))
            if 'pairs' in z: sp[d] = z
    return sp


# ---------- sim ----------
def sim_noise(sp, r, cell=MAIN):
    """쌍의 쓸 수 있는 온도, u, 반복 잡음(log 퍼짐의 반복 사이 SE: 데워진 쪽 T 별, 기준 쪽 320)."""
    key = f"{r['a']}__{r['b']}"; H_, R_ = sp.get(r['hot']), sp.get(r['ref'])
    if H_ is None or R_ is None or key not in H_['pairs'] or key not in R_['pairs'] or H_['pairs'][key]['n_shared'] < MINSH1: return None
    hc, rc = H_['pairs'][key]['cells'], R_['pairs'][key]['cells']; k320 = ckey('320', *cell)
    if k320 not in rc or k320 not in hc: return None
    se = lambda runs: float(np.std(np.log(runs), ddof=1) / np.sqrt(len(runs))) if len(runs) > 1 else 0.0
    temps, ses = [], []
    for T in TT:
        kk = ckey(T, *cell)
        if kk in hc: temps.append(TK[T]); ses.append(se(hc[kk]['runs']))
    if len(temps) < 2: return None
    return dict(H=r['H'], temps=temps, se_h=ses, se_r=se(rc[k320]['runs']))


def sim_once(nz, theta0, sd_theta, rng, B=NB1):
    items = {}
    for z in nz:
        th = theta0 + (rng.normal(0, sd_theta) if sd_theta > 0 else 0.0); er = rng.normal(0, z['se_r'])
        LR = [B_SLOPE * (T - th) + rng.normal(0, s) - er for T, s in zip(z['temps'], z['se_h'])]
        items.setdefault(z['H'], []).append(theta(LR, z['temps'])[:2])
    T_s = {h: sf_summary(v, 'ht') for h, v in items.items()}; it = list(T_s.values()); med = wmedian(it, key_ht); lo, hi = boot_median(it, key_ht, rng, B)
    return med, lo, hi


def part_sim():
    pj = json.load(open(OUT1 + '/pairs.json')); sp = load_spread(pj['domains']); nz = [z for z in (sim_noise(sp, r) for r in pj['HT']) if z is not None]
    rng = np.random.default_rng(20261042); out = dict(written=time.strftime('%F %T'), n_pairs=len(nz), n_sf=len({z['H'] for z in nz}), slope=B_SLOPE, reps=SIM_REP,
                                                     se_h_median={T: float(np.median([s for z in nz for t, s in zip(z['temps'], z['se_h']) if t == TK[T]] or [np.nan])) for T in TT},
                                                     se_r_median=float(np.median([z['se_r'] for z in nz])), rows=[])
    for sdv in (0.0, 15.0):
        for th0 in SIM_THETA:
            cov, bias, lob = [], [], []
            for _ in range(SIM_REP):
                med, lo, hi = sim_once(nz, th0, sdv, rng)
                lo_ok = key_ht(lo) <= (th0, 1) if lo[0] != 'lt' else True; hi_ok = True if hi[0] == 'gt' else (hi[0] != 'lt' and hi[1] >= th0)
                cov.append(lo_ok and hi_ok); lob.append(key_ht(lo) > (320.0, 0))
                if med[0] == 'num': bias.append(med[1] - th0)
            out['rows'].append(dict(theta=th0, sd_theta=sdv, coverage=float(np.mean(cov)), frac_lo_gt_320=float(np.mean(lob)), bias_median_num=float(np.median(bias)) if bias else None, n_num=len(bias)))
            log('p1 sim', th0, sdv, out['rows'][-1])
    out['coverage_min'] = float(min(r['coverage'] for r in out['rows'] if r['theta'] < 413)); out['rule'] = '피복률(θ < 413 시나리오) < 0.90 이면 구간 방법을 고친 뒤 동결; 치우침은 서술에만'
    jdump(out, RES + '/r38_p1_sim.json', bk=True); log('P1_SIM_DONE', out['coverage_min'])


# ---------- pre ----------
def part_pre():
    pj = json.load(open(OUT1 + '/pairs.json')); sp = load_spread(pj['domains']); rows = []
    for r in pj['PRE']:
        a, b = pair_theta(sp, r, MAIN, 'gpa'), pair_theta(sp, r, MAIN, 'full')
        if a[0] != 'skip' and b[0] != 'skip': rows.append(dict(H=r['H'], gpa=a[:2], full=b[:2]))
    num = [(x['gpa'][1], x['full'][1]) for x in rows if x['gpa'][0] == 'num' and x['full'][0] == 'num']; dd = np.array([f - g for g, f in num]) if num else np.array([])
    out = dict(written=time.strftime('%F %T'), n_pairs=len(pj['PRE']), n_used=len(rows), n_both_numeric=len(num), n_sf=len({x['H'] for x in rows}),
               diff_full_minus_gpa=dict(median=float(np.median(dd)) if len(dd) else None, iqr=np.percentile(dd, [25, 75]).tolist() if len(dd) else None, mean=float(dd.mean()) if len(dd) else None),
               censor_agree=float(np.mean([x['gpa'][0] == x['full'][0] for x in rows])) if rows else None,
               rule='prereg 7 절 사전 측정: T-풀 밖 H-풀 원핵 중온–중온 쌍(방향 무작위)에서 두 규약의 θ* 차. 기록이지 판정·규칙 조정에 쓰지 않음')
    jdump(out, RES + '/r38_p1_pre.json', bk=True); log('P1_PRE_DONE', out['diff_full_minus_gpa'])


# ---------- stats ----------
def logistic_w(X, y, iters=50):
    X = np.column_stack([np.ones(len(X)), (X - X.mean(0)) / (X.std(0) + 1e-12)]); b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p); b += np.linalg.solve(X.T @ (X * W[:, None]) + 1e-8 * np.eye(len(b)), X.T @ (y - p))
    return 1 / (1 + np.exp(-X @ b))


def part_stats():
    pj = json.load(open(OUT1 + '/pairs.json')); HT, MM = pj['HT'], pj['MM']; sp = load_spread(pj['domains'])
    main, tH, tM = verdict_set(HT, MM, sp); out = dict(written=time.strftime('%F %T'), main=main)
    sens = {}
    for name, kw in (('full_superposition', dict(conv='full')), ('q070', dict(cell=('half', 0.7))), ('q080', dict(cell=('half', 0.8))), ('abs_50_100ns', dict(cell=('abs', 0.75))), ('no413', dict(drop413=True))):
        sens[name] = verdict_set(HT, MM, sp, seed=20261043, **kw)[0]
    # 성향 점수 재가중 (짝지은 상과의 두 종류 쌍; TM·동일성·길이 차·공유 변 수)
    sfb = set(pj['sf_both']); ns = lambda r: (sp.get(r['hot'], {}).get('pairs', {}).get(f"{r['a']}__{r['b']}", {}).get('n_shared', np.nan))
    rr = [(r, 1) for r in HT if r['H'] in sfb] + [(r, 0) for r in MM if r['H'] in sfb]; Xp = np.array([[r['tm'], r['ident'], abs(r['La'] - r['Lb']), ns(r)] for r, _ in rr], float)
    ok = np.isfinite(Xp).all(1); wmm = {}
    if ok.sum() > 10:
        ps = logistic_w(Xp[ok], np.array([y for _, y in rr], float)[ok]); i_ok = np.where(ok)[0]
        raw = {(rr[i][0]['a'], rr[i][0]['b']): ps[k] / (1 - ps[k]) for k, i in enumerate(i_ok) if rr[i][1] == 0}
        cap = 5 * np.median(list(raw.values())); wmm = {k: min(v, cap) for k, v in raw.items()}
    wfull = {(r['a'], r['b']): wmm.get((r['a'], r['b']), 1.0) for r in MM}
    sens['propensity_reweight'] = verdict_set(HT, MM, sp, wmm=wfull, seed=20261044)[0]; sens['propensity_reweight']['n_weighted'] = len(wmm)
    dist = {k: {v: dict(zip(('p10', 'p50', 'p90'), np.nanpercentile([f(r) for r in g], [10, 50, 90]).round(3).tolist())) for v, f in
                (('tm', lambda r: r['tm']), ('ident', lambda r: r['ident']), ('dlen', lambda r: abs(r['La'] - r['Lb'])), ('n_shared', ns))} for k, g in (('HT', HT), ('MM', MM))}
    # θ* 흔들림 (주 칸, 거의 최적 대응들)
    fl = []
    for r in HT + MM:
        key = f"{r['a']}__{r['b']}"; H_ = sp.get(r['hot'])
        if H_ is None or key not in H_['pairs']: continue
        k320 = ckey('320', *MAIN); c = H_['pairs'][key]['cells'].get(k320)
        if c is None or not c.get('alts'): continue
        vals = [pair_theta(sp, r, MAIN, alt=j) for j in range(len(c['alts']))] + [pair_theta(sp, r, MAIN)]
        nums = [v[1] for v in vals if v[0] == 'num']
        if len(nums) >= 2: fl.append(max(nums) - min(nums))
    # 궤적 길이(온도별, 데워진 쪽), 2차 용량–반응
    tl = {k: {T: float(np.median([x for r in g if r['hot'] in sp for x in sp[r['hot']]['lens'].get(T, [])] or [np.nan])) for T in TT} for k, g in (('HT', HT), ('MM', MM))}
    num = [(r, t) for r, t in zip(HT, tH) if t[0] != 'skip']; dr = None
    if num:
        val = lambda t: t[1] - 320 if t[0] == 'num' else (t[1] - 320 + 1e-3 if t[0] == 'gt' else -1e-3 - 1e3)
        H_ = np.array([r['H'] for r, _ in num]); x = np.array([r['dOGT'] for r, _ in num]); y = np.array([val(t) for _, t in num]); rx, ry = np.zeros(len(x)), np.zeros(len(y))
        for h in set(H_):
            s_ = H_ == h
            if s_.sum() >= 3: rx[s_] = ranks(x[s_]) - ranks(x[s_]).mean(); ry[s_] = ranks(y[s_]) - ranks(y[s_]).mean()
        m3 = np.isin(H_, [h for h in set(H_) if (H_ == h).sum() >= 3]); dr = dict(n=int(m3.sum()), n_sf=int(len(set(H_[m3]))), within_sf_spearman=float(np.corrcoef(rx[m3], ry[m3])[0, 1]) if m3.sum() > 3 else None)
    nH = [t[1] for t in tH if t[0] == 'num']
    out.update(sensitivity=sens, matching=dist, align_uncertainty=dict(n=len(fl), median_K=float(np.median(fl)) if fl else None, p90_K=float(np.percentile(fl, 90)) if fl else None),
               traj_len_median=tl, dose_response_2nd=dr, frac_340_350=float(np.mean([(340 <= v <= 350) for v in nH])) if nH else None,
               rule='prereg 7 절. P1 = 조건 1(상과 중앙값 95% 하한 > 320 K_mdCATH, 하한 순위) 과 조건 2(짝지은 상과 차 d_s 의 상과 중앙값 95% 하한 > 0) 모두. 치우침 보정 없음',
               tag='다시 재는 것: 앞 판 θ* − 320 = +21 K [6.7, 36.4], 대칭 판 −1.7 K (결과 3.5)')
    jdump(out, RES + '/r38_p1.json', bk=True); log('P1', dict(cond1=main['cond1']['passed'], cond2=main['cond2']['passed'], P1=main['P1'])); log('P1_STATS_DONE')


ALL1 = ['pairs', 'mean', 'align', 'spread', 'sim', 'pre', 'stats']
if __name__ == '__main__':
    for part in (ALL1 if PARTS1 == ['all'] else PARTS1):
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('P1_DONE', ','.join(PARTS1))
