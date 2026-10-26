"""(prereg 7h, 확인 실험 P1′) 새 쌍에서 열 적응을 안정성과 유연성으로 나눠 다시. r38_p1.py 의 함수를 그대로 쓴다(그 파일은 바꾸지 않음).
입력: derived/r38/p1n/pairs.json (등록 2026-10-20 13:35:54 UTC), MD = raw/mdcath_p1n/<d>.npz (r38_p1n_fetch) 가 있으면 그것, 없으면 h5 (build_cells 와 같은 프레임).
부분: mean → align (CPU 작업자 여럿; 계산은 r38_p1.align_batch 그대로) → spread → fold → stats 사용: python r38_p1n.py <ROOT> <부분들> [작업자=4]
  fold : 도메인·온도(320·348·379·413)마다 뒤 절반 프레임의 접힌 비율 f_q (Q ≥ 0.7·0.75·0.8; 320 K 첫 프레임 천연 접촉), 반복 평균.
  stats: prereg 7h — P1′-1 U2(Q 0.7) 상한 < 0, P1′-2 U2(0.7) − U2(0.8) 상한 < 0, P1′-3 접힌 비율 logit 차(379 K 주, 413 보조)의 짝지은 상과 차 하한 > 0.
         상과 부트스트랩 2,000 (seed 20261049), 한쪽 97.5% 한계. 서술: U2(0.75·0.8·50–100 ns)·U1, 1차 P1 조건 1·2, 상과 독립(indep) 부분.
출력: derived/r38/p1n/{mean/, align.jsonl, spread/, fold/}, results/r38_p1n.json"""
import os, sys, json
ROOT = sys.argv[1]; PARTSN = sys.argv[2].split(','); NWN = int(sys.argv[3]) if len(sys.argv) > 3 else 4
_s = open(ROOT + '/code/r38_p1.py').read(); _s = _s[:_s.index("ALL1 = [")]
sys.argv = ['r38_p1', ROOT, 'none', str(NWN)]; exec(_s)
OUT1 = ROOT + '/derived/r38/p1n'; RAWN = ROOT + '/raw/mdcath_p1n'
for p_ in (OUT1 + '/mean', OUT1 + '/spread', OUT1 + '/fold'): os.makedirs(p_, exist_ok=True)
build_cells_h5 = build_cells; QF3 = (0.7, 0.75, 0.8)


def _npz_q(z):
    smd = ''.join(z['smd'].tolist()); L = len(smd); ca0 = z['first'][:, 0]
    ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    return smd, (lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1))


def build_cells_npz(d, temps=TT, cells=CELLS):
    """r38_p1.build_cells 와 같은 계산; 입력만 npz (뒤 절반 h<T>_<r>, 50–100 ns a<T>_<r>). 같은 이름 seed·같은 상대 위치라 같은 프레임."""
    z = np.load(f'{RAWN}/{d}.npz'); smd, Qf = _npz_q(z); out, lens = {}, {}
    for T in temps:
        if f'reps{T}' not in z.files: continue
        reps = [int(r) for r in z[f'reps{T}']]; lens[T] = [int(n) for n in z[f'nframes{T}']]
        H_ = [z[f'h{T}_{r}'] for r in reps]; A_ = [z[f'a{T}_{r}'] if f'a{T}_{r}' in z.files else None for r in reps]
        qh = [Qf(X) for X in H_]; qa = [Qf(X) if X is not None else None for X in A_]
        for w, q in cells:
            N = NWIN[w]; rng = np.random.default_rng(zlib.crc32(f'r38_p1|{d}|{T}|{w}|{q:.2f}'.encode())); runs, lab = [], []
            for r, (Xh, Xa, vh, va) in enumerate(zip(H_, A_, qh, qa)):
                X, v = (Xh, vh) if w == 'half' else (Xa, va)
                if X is None: continue
                idx = np.where(v >= q)[0]
                if len(idx) >= N: runs.append(X[np.sort(rng.choice(idx, N, replace=False))]); lab += [r] * N
            if len(runs) >= MINRUN:
                E = frames(np.concatenate(runs)); m = np.asarray(fs.mean_point(E)); A = np.asarray(fs.align_to_reference(E, torch.as_tensor(m)), dtype=np.float64)
                out[ckey(T, w, q)] = (A, np.array(lab), m)
    return out, lens, smd


def build_cells(d, temps=TT, cells=CELLS):
    return build_cells_npz(d, temps, cells) if os.path.exists(f'{RAWN}/{d}.npz') else build_cells_h5(d, temps, cells)


# ---------- align (CPU 작업자 여럿) ----------
MSA = {}; NMCAP = 2.2e5 # 묶음 상한(쌍 수 × n × m): 작업자 하나 약 5 GB (2026-10-20 16:24 UTC 측정; 처음 판 8 쌍 고정은 긴 쌍에서 RAM 을 넘길 어림이라 바꿈)
def _align_worker(B):
    torch.set_num_threads(1); return [(list(p), r) for p, r in zip(B, align_batch(B, MSA, 'cpu'))]


def part_align():
    global MSA
    from multiprocessing import Pool
    fo = OUT1 + '/align.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    Ps = {(r['a'], r['b']) for r in all_pairs()}
    for d in sorted({x for p in Ps for x in p}):
        m = mean_of(d)
        if m is not None: MSA[d] = m
    todo = sorted([p for p in Ps if p not in done and p[0] in MSA and p[1] in MSA], key=lambda p: (len(MSA[p[0]]) // 16, len(MSA[p[1]]) // 16)); log('p1n align', len(todo), 'cpu', NWN, '작업자')
    batches, cur = [], [] # 묶음 = 최대 8 쌍이고 쌍 수 × 최대 n × 최대 m ≤ NMCAP (메모리; 쌍마다 seed 가 따로라 결과는 묶음과 무관)
    for p in todo:
        nxt = cur + [p]
        if cur and (len(nxt) > 8 or len(nxt) * max(len(MSA[q[0]]) for q in nxt) * max(len(MSA[q[1]]) for q in nxt) > NMCAP): batches.append(cur); nxt = [p]
        cur = nxt
    if cur: batches.append(cur)
    t0 = time.time(); k = 0
    with Pool(NWN) as pool, open(fo, 'a') as fh:
        for res in pool.imap_unordered(_align_worker, batches):
            for p, (S, best, alts, nn) in res: fh.write(json.dumps(dict(ab=p, S=S, best=best, alts=alts, n_near=nn)) + '\n')
            fh.flush(); k += 1
            if k % 20 == 0: log('p1n align', k, '/', len(batches), '묶음', round(time.time() - t0), 's')
    log('P1N_ALIGN_DONE')


# ---------- fold ----------
def _fold_worker(d):
    fo = f'{OUT1}/fold/{d}.json'
    if os.path.exists(fo): return d, 'skip'
    try:
        res = dict(dom=d, f={})
        if os.path.exists(f'{RAWN}/{d}.npz'):
            z = np.load(f'{RAWN}/{d}.npz'); _, Qf = _npz_q(z)
            for T in TT:
                if f'reps{T}' not in z.files: continue
                qs = [Qf(z[f'h{T}_{r}']) for r in [int(x) for x in z[f'reps{T}']]]
                res['f'][T] = {f'{q:.2f}': float(np.mean([(v >= q).mean() for v in qs])) for q in QF3}
        else:
            h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode()); L = len(smd)
            ld = lambda T: [g[T][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g[T].keys(), key=int)]
            ref = ld('320'); ca0 = ref[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
            Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
            for T in TT:
                if T not in g: continue
                reps = ref if T == '320' else ld(T); qs = [Qf(X[len(X) // 2:]) for X in reps]
                res['f'][T] = {f'{q:.2f}': float(np.mean([(v >= q).mean() for v in qs])) for q in QF3}
            h.close()
        json.dump(res, open(fo, 'w')); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


def part_fold():
    from multiprocessing import Pool
    doms = json.load(open(OUT1 + '/pairs.json'))['domains']
    with Pool(NWN) as p:
        for k, (d, s_) in enumerate(p.imap_unordered(_fold_worker, doms)):
            if s_ not in ('skip', 'ok') or k % 100 == 0: log('p1n fold', k, '/', len(doms), d, s_)
    log('P1N_FOLD_DONE')


# ---------- stats ----------
def lr320_of(SP, r, cell):
    key = f"{r['a']}__{r['b']}"; H_, R_ = SP.get(r['hot']), SP.get(r['ref'])
    if H_ is None or R_ is None or key not in H_['pairs'] or key not in R_['pairs'] or H_['pairs'][key]['n_shared'] < MINSH1: return None
    k = ckey('320', *cell); hc, rc = H_['pairs'][key]['cells'], R_['pairs'][key]['cells']
    if k not in hc or k not in rc or hc[k]['gpa'] is None or rc[k]['gpa'] is None: return None
    return float(np.log(hc[k]['gpa'] / rc[k]['gpa']))


def sf_med(prs, f):
    by = {}
    for r in prs:
        v = f(r)
        if v is not None and np.isfinite(v): by.setdefault(r['H'], []).append(v)
    return {h: float(np.median(v)) for h, v in by.items()}


def ustats(HT, MM, SP, rng, cells=(('half', 0.7), ('half', 0.75), ('half', 0.8), ('abs', 0.75)), B=2000):
    """U1·U2 (점, 2.5·97.5%) 와 P1′-2 의 짝지은 부트스트랩(같은 상과 재추출에서 U2(0.7) − U2(0.8))."""
    T = {c: sf_med(HT, lambda r: lr320_of(SP, r, c)) for c in cells}; C = {c: sf_med(MM, lambda r: lr320_of(SP, r, c)) for c in cells}; out = {}
    for c in cells:
        tv = np.array(list(T[c].values())); both = sorted(set(T[c]) & set(C[c])); dv = np.array([T[c][h] - C[c][h] for h in both])
        bt = [np.median(tv[rng.integers(0, len(tv), len(tv))]) for _ in range(B)] if len(tv) else []; bd = [np.median(dv[rng.integers(0, len(dv), len(dv))]) for _ in range(B)] if len(dv) else []
        out[f'{c[0]}|{c[1]:.2f}'] = dict(U1=dict(n_sf=len(tv), median=float(np.median(tv)) if len(tv) else None, lo=float(np.percentile(bt, 2.5)) if bt else None, hi=float(np.percentile(bt, 97.5)) if bt else None),
                                         U2=dict(n_sf=len(dv), median=float(np.median(dv)) if len(dv) else None, lo=float(np.percentile(bd, 2.5)) if bd else None, hi=float(np.percentile(bd, 97.5)) if bd else None))
    c7, c8 = ('half', 0.7), ('half', 0.8); sfs = sorted(set(T[c7]) & set(C[c7]) & set(T[c8]) & set(C[c8]))
    d7 = np.array([T[c7][h] - C[c7][h] for h in sfs]); d8 = np.array([T[c8][h] - C[c8][h] for h in sfs]); bs = []
    for _ in range(B):
        ix = rng.integers(0, len(sfs), len(sfs)); bs.append(np.median(d7[ix]) - np.median(d8[ix]))
    out['decomp'] = dict(n_sf=len(sfs), U2_07=float(np.median(d7)) if len(sfs) else None, U2_08=float(np.median(d8)) if len(sfs) else None,
                         diff=float(np.median(d7) - np.median(d8)) if len(sfs) else None, lo=float(np.percentile(bs, 2.5)) if bs else None, hi=float(np.percentile(bs, 97.5)) if bs else None)
    return out


def fold_stats(HT, MM, FD, rng, B=2000):
    lg = lambda f: float(np.log(np.clip(f, 0.005, 0.995) / (1 - np.clip(f, 0.005, 0.995))))
    def pv(r, T):
        a, b = FD.get(r['hot'], {}).get(T), FD.get(r['ref'], {}).get(T)
        return lg(a['0.75']) - lg(b['0.75']) if a and b else None
    out = {}
    for T in ('379', '413'):
        Ts, Cs = sf_med(HT, lambda r: pv(r, T)), sf_med(MM, lambda r: pv(r, T)); both = sorted(set(Ts) & set(Cs)); dv = np.array([Ts[h] - Cs[h] for h in both])
        bd = [np.median(dv[rng.integers(0, len(dv), len(dv))]) for _ in range(B)] if len(dv) else []
        out[T] = dict(n_sf=len(dv), median=float(np.median(dv)) if len(dv) else None, lo=float(np.percentile(bd, 2.5)) if bd else None, hi=float(np.percentile(bd, 97.5)) if bd else None,
                      HT_sf_median=float(np.median(list(Ts.values()))) if Ts else None, MM_sf_median=float(np.median(list(Cs.values()))) if Cs else None)
    return out


def part_stats():
    pj = json.load(open(OUT1 + '/pairs.json')); SP = load_spread(pj['domains']); FD = {}
    for d in pj['domains']:
        f = f'{OUT1}/fold/{d}.json'
        if os.path.exists(f): FD[d] = json.load(open(f))['f']
    res = dict(written=time.strftime('%F %T'))
    for nm, sel in (('all', lambda r: True), ('indep', lambda r: r['set'] == 'indep')):
        HT = [r for r in pj['HT'] if sel(r)]; MM = [r for r in pj['MM'] if sel(r)]; rng = np.random.default_rng(20261049)
        U = ustats(HT, MM, SP, rng); F = fold_stats(HT, MM, FD, rng); V = verdict_set(HT, MM, SP, seed=20261050)[0]
        res[nm] = dict(U=U, fold=F, p1_conditions=dict(cond1=V['cond1'], cond2=V['cond2'], P1=V['P1'], pairs_HT=V['pairs_HT'], pairs_MM=V['pairs_MM']))
    a = res['all']; u7 = a['U']['half|0.70']['U2']; dc = a['U']['decomp']; f379 = a['fold']['379']
    p1 = bool(u7['hi'] is not None and u7['hi'] < 0); p2 = bool(dc['hi'] is not None and dc['hi'] < 0); p3 = bool(f379['lo'] is not None and f379['lo'] > 0)
    res['verdict'] = dict(P1p_1=p1, P1p_2=p2, P1p_3=p3,
                          sentence=('호열–중온의 유연성 차는 근천연 들뜸을 넣을 때 나타나고 엄격한 천연 상태에서는 작다 — 새 쌍에서 재현된다' if p1 and p2 else
                                    '같은 방향으로 재현되지만 안정성·유연성 분해는 확인되지 않았다' if p1 else '새 쌍에서 재현되지 않았다') + (' / 같은 온도에서 호열 쪽이 접힌 상태에 더 오래 머문다' if p3 else ''),
                          rule='prereg 7h: P1′-1 U2(Q 0.7) 97.5% 상한 < 0; P1′-2 U2(0.7) − U2(0.8) 97.5% 상한 < 0; P1′-3 접힌 비율 logit 차 짝지은 상과 차(379 K) 2.5% 하한 > 0. 상과 부트스트랩 2,000 (seed 20261049)')
    jdump(res, RES + '/r38_p1n.json', bk=True); log('P1N', res['verdict']); log('P1N_STATS_DONE')


ALLN = ['mean', 'align', 'spread', 'fold', 'stats']
if __name__ == '__main__':
    for part in (ALLN if PARTSN == ['all'] else PARTSN):
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('P1N_DONE', ','.join(PARTSN))
