"""(prereg round36 0 절 첫 세션 + 0.0b E7 5 도메인) 계산의 바름. 2026-10-19.
부분(prereg 0 절 순서): pick → b0 → e7 → trunc → ens → bias → ustat → unc
  pick  : 고정 seed 로 E7 5 도메인(P2 67 − 7c 바름 점검 5), 잘라 낸 사슬 10 도메인(E7 5 + 5), 첫 세션 30 쌍
          (r37 후보 목록 − 281 쌍, TM 긴 쪽 7 구간 층 5·5·5·5·5·3·2; 두 도메인 모두 320 K 접힌 실행 ≥ 2)을 뽑는다.
  b0    : 0.0 BnB 검증. 30 쌍 가운데 긴 쪽 길이가 짧은 10 쌍의 메도이드 프레임. BnB 값이 128 출발 값 이상(d 로는 이하)이고,
          무작위 회전 격자(10만 회전의 DP 점수 + 상위 256 에서 번갈아 다듬기)와 sim_long 0.005 안에서 맞는지.
  e7    : 0.0b. 같은 표본(실행마다 N = 100, 실행 쌍마다 프레임 쌍 2,000)에서 Δ², w̄, 비를 D(줄인 계산)와 맞춘 좌표로 낸다.
          칸 = MD T 대 320 K (348·379·413·450; 전체·뒤 절반) + R0 (320 K 반복 2 대 3, 10 가지).
  trunc : 0.1 잘라 낸 사슬. 같은 궤적의 좌표에서 끝 10 변만 뺀 앙상블. 교차는 다른 반복끼리만.
          빈칸 몫(전체 대 잘라 낸) − 빈칸 몫(전체 대 전체) = 0.75 Σs²(자른 변), 짝 몫의 Δ² 는 R0 띠 안.
  ens   : 30 쌍의 Δ²_D (교차 = 모든 실행 쌍, 바닥 = 같은 쪽 실행 쌍, 실행 쌍마다 프레임 쌍 2,000). bias 의 구간 Δ²_D 로 쓴다.
  bias  : 0.2 상계 편향. TM 구간마다 교차 프레임 쌍 200·바닥 프레임 쌍 200 의 줄인 계산 d 를 BnB 값과 견준다(구간을 번갈아 돈다).
          첫 세션은 30 쌍의 7 구간만; R10a′·R6 은 본 계산에서.
  ustat : 0.3 불완전 U-통계량. 30 쌍 앞 10 쌍의 첫 교차 실행 쌍에서 2,000 쌍 평균 대 전수(100 × 100) 평균.
  unc   : 0.4 정렬 불확실성. 30 쌍의 메도이드: 128 출발 가운데 S ≥ (1 − 0.005) S_최적 인 서로 다른 대응들 사이의
          짝 몫·빈칸 몫 흔들림(최대 − 최소, 변당). θ* 흔들림은 P1 코드(트랙 B5)에서 붙인다.
줄인 계산 = framespace_wip/stage0 의 평균 점 출발(MeanStart; 평균 점끼리 128 출발) + 띠 DP(반폭 5, 가장자리에 닿은 쌍은 반폭 두 배,
  열 수까지) — stage0.w2.pair_costs 의 band 가지와 같은 규칙을 쌍 목록에 쓴 것. float64, GPU.
접힌 프레임 = 첫 320 K 프레임의 천연 접촉(Cα, |i − j| ≥ 4, < 8 Å)의 Q ≥ 0.75 (r37_ff·r38_md_ref 와 같다).
실행 = 반복. 창 안 접힌 프레임 ≥ 100 인 반복에서 100 개를 무작위로(도메인·온도·창 이름 seed, 모든 부분이 같은 실행을 쓴다).
사용: python r38_s0.py <ROOT> <부분들 쉼표 | all>
출력: derived/r38/s0/<부분>/..., results/r38_s0_<부분>.json (결과마다 results/_backup/<날짜>/ 에 사본)"""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '4'
import sys, json, glob, time, itertools, zlib, shutil, traceback
ROOT = sys.argv[1]; PARTS = sys.argv[2].split(',') if len(sys.argv) > 2 else ['all']; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])     # np, torch, fs, frames, cost, TEMPS, QT
src = open(ROOT + '/code/r37_ff.py').read(); exec(src[src.index('AA = dict('):src.index('def stats')])   # AA, TT, md_atoms
import h5py, difflib
from scipy.stats import spearmanr
from stage0.band import MeanStart, align_band, band_windows
from framespace.distance import sq_norm, dp_score
from framespace.geometry import procrustes, random_rotations

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'; DT = torch.float64
OUT = ROOT + '/derived/r38/s0'; RES = ROOT + '/results'; os.makedirs(OUT, exist_ok=True)
N = 100; NPAIR = 2000; K0 = 5; REL = 0.005; CUT = 10
SEED = dict(e7=20261020, pairs=20261021, trunc=20261022)
BINS = [('0.2-0.3', 0.2, 0.3, 5), ('0.3-0.4', 0.3, 0.4, 5), ('0.4-0.5', 0.4, 0.5, 5), ('0.5-0.6', 0.5, 0.6, 5),
        ('0.6-0.7', 0.6, 0.7, 5), ('0.7-0.8', 0.7, 0.8, 3), ('0.8-1.01', 0.8, 1.01, 2)]
R0_BAND = (-0.0033, 0.0006)     # 결과 문서 3.1: R0 replica 띠(Δ²/√변, n = 314 의 95%)


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def backup(p):
    b = ROOT + '/results/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True)
    shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + os.path.basename(p))


def jdump(o, p, bk=False):
    json.dump(o, open(p + '.tmp', 'w'), ensure_ascii=False, indent=1, default=float); os.replace(p + '.tmp', p)
    if bk: backup(p)


def seeded(*tag):
    return np.random.default_rng(zlib.crc32('|'.join(map(str, tag)).encode()))


# ---------- 자료 ----------
_md = {}
def md_load(d, temps=('320',)):
    """mdCATH h5 -> ({T: [반복 (F, L, 3, 3) Å (CA, C, O)]}, Qf, L). 320 K 은 늘 읽는다(천연 접촉)."""
    temps = tuple(sorted(set(temps) | {'320'}))
    if (d, temps) in _md: return _md[(d, temps)]
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode()); L = len(smd)
    MD = {T: [g[T][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g[T].keys(), key=int)] for T in temps if T in g}; h.close()
    ca0 = MD['320'][0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
    if len(_md) >= 3: _md.clear()
    _md[(d, temps)] = (MD, Qf, L); return _md[(d, temps)]


_runs = {}
def dom_runs(d, T='320', W='full', temps=('320',)):
    """[(반복 번호, (N, L, 3, 3) 원자)] 와 L. 창 안 접힌 프레임 ≥ N 인 반복만."""
    key = (d, T, W)
    if key not in _runs:
        MD, Qf, L = md_load(d, temps); rng = seeded('r38_s0_runs', d, T, W); out = []
        for r, X in enumerate(MD.get(T, [])):
            Xw = X if W == 'full' else X[len(X) // 2:]
            idx = np.where(Qf(Xw) >= QT)[0]
            if len(idx) >= N: out.append((r, Xw[np.sort(rng.choice(idx, N, replace=False))]))
        _runs[key] = (out, L)
    return _runs[key]


def coords(X):
    return torch.as_tensor(frames(np.ascontiguousarray(X))).to(device=DEV, dtype=DT)


def side(d, T='320', W='full', temps=('320',)):
    R, L = dom_runs(d, T, W, temps)
    return (coords(np.concatenate([X for _, X in R])) if R else None), len(R), L


# ---------- 줄인 계산 ----------
def red_pairs(ms, ia, ib, k0=K0, iters=10, want_gap=False, membudget=3e9):
    """줄인 계산을 쌍 목록 (ia[p], ib[p]) 에. ms = MeanStart(EA, EB).
    반환: d (P,), 빈칸 몫 (P,) 또는 None, S (P,), 마지막 반폭 (P,), 처음 띠에서 넓힌 쌍 수."""
    EA, EB = ms.EA, ms.EB; n, m = EA.shape[1], EB.shape[1]
    ia = torch.as_tensor(np.asarray(ia), device=EA.device); ib = torch.as_tensor(np.asarray(ib), device=EA.device)
    nA, nB = (EA ** 2).sum((-1, -2)), (EB ** 2).sum((-1, -2))
    S = EA.new_zeros(len(ia)); G = EA.new_zeros(len(ia)) if want_gap else None; kused = np.zeros(len(ia), int); widened = 0
    todo = torch.arange(len(ia), device=EA.device); kk = k0
    while len(todo):
        lo, w = band_windows(ms.center, m, kk); step = max(1, int(membudget / (n * w * 9 * 8 * 4)))
        edge_all = torch.zeros(len(todo), dtype=torch.bool, device=EA.device)
        for c in range(0, len(todo), step):
            sel = todo[c:c + step]; a, b = EA[ia[sel]], EB[ib[sel]]
            sc, _, M, edge = align_band(a, b, lo, w, ms.starts(ia[sel], ib[sel], n_local=1), iters=iters)
            S[sel] = sc; edge_all[c:c + step] = edge
            if want_gap: G[sel] = (nA[ia[sel]] * (~M.any(2))).sum(1) + (nB[ib[sel]] * (~M.any(1))).sum(1)
        kused[todo.cpu().numpy()] = kk
        if w >= m: break
        todo = todo[edge_all]
        if kk == k0: widened = int(len(todo))
        kk *= 2
    d2 = (sq_norm(EA)[ia] + sq_norm(EB)[ib] - 2 * S).clamp_min(0)
    return d2.sqrt().cpu().numpy(), (G.cpu().numpy() if want_gap else None), S.cpu().numpy(), kused, widened


def upairs(a, b, rng, k=NPAIR, nn=N):
    """실행 a (그 쪽 앙상블에서 a 번째 실행), 실행 b 의 프레임 쌍 k 개(비복원) -> 앙상블 안 색인."""
    f = rng.choice(nn * nn, k, replace=False); return a * nn + f // nn, b * nn + f % nn


def runpair_means(ms, rr, rng, want_gap=False):
    """실행 쌍 목록 rr 마다 프레임 쌍 2,000 의 평균 d (와 빈칸 몫·짝 몫 거리). 한 번에 푼다."""
    IA, IB, GR = [], [], []
    for g, (a, b) in enumerate(rr):
        ia, ib = upairs(a, b, rng); IA.append(ia); IB.append(ib); GR.append(np.full(len(ia), g))
    IA, IB, GR = map(np.concatenate, (IA, IB, GR))
    d, gap, S, ku, wid = red_pairs(ms, IA, IB, want_gap=want_gap)
    out = dict(d=[float(d[GR == g].mean()) for g in range(len(rr))], widened=wid, kmax=int(ku.max()), IA=IA, IB=IB, GR=GR, dv=d)
    if want_gap:
        pair = np.sqrt(np.maximum(d ** 2 - gap, 0))
        out.update(gap=[float(gap[GR == g].mean()) for g in range(len(rr))], dpair=[float(pair[GR == g].mean()) for g in range(len(rr))], gv=gap)
    return out


def estats(C, BP, BQ, ne=None):
    D2 = 2 * C - BP - BQ; w = (BP + BQ) / 2
    o = dict(C=float(C), B_P=float(BP), B_Q=float(BQ), D2=float(D2), wbar=float(w), ratio=float(np.sqrt(max((C / w) ** 2 - 1, 0))) if w > 0 else float('nan'))
    if ne: o['d2n'] = float(D2 / np.sqrt(ne))
    return o


def medoid(d):
    E, nr, L = side(d)
    return E[int(fs.fixed_sq_distance(E, E).sum(1).argmin())]


def align_all(x, y, n_random=127, iters=10, seed=0, improve_tol=1e-7):
    """fs.align 과 같은 계산(항등 + 무작위 127 출발, 번갈아 풀기)을 한 쌍에 하고, 출발마다의 최종 (S, U, 짝) 을 모두 돌려준다."""
    x, y = x[None], y[None]
    I = torch.eye(3, dtype=x.dtype, device=x.device).expand(1, 1, 3, 3)
    U = torch.cat([I, random_rotations(n_random, seed=seed, dtype=x.dtype, device=x.device).view(1, n_random, 3, 3)], 1)
    best = torch.full((1, U.shape[1]), -1e30, dtype=x.dtype, device=x.device); bU, bG = U.clone(), None
    for it in range(iters + 1):
        with torch.no_grad():
            K = torch.einsum('pnxz,psmxz->psnm', x, torch.einsum('psxy,pmyz->psmxz', U, y))
        K = K.detach().requires_grad_(True)
        with torch.enable_grad():
            Sc = dp_score(K); Gr, = torch.autograd.grad(Sc.sum(), K)
        with torch.no_grad():
            imp = Sc > best + improve_tol * (1 + best.abs())
            best = torch.where(imp, Sc, best); bU = torch.where(imp[..., None, None], U, bU)
            bG = Gr if bG is None else torch.where(imp[..., None, None], Gr, bG)
            if it == iters or not imp.any(): break
            U = torch.where(imp[..., None, None], procrustes(torch.einsum('psnm,pmxz,pnyz->psxy', Gr, y, x)), U)
    return best[0], bU[0], bG[0] > 0.5


# ---------- pick ----------
def avail(d, memo={}):
    if d not in memo:
        try: memo[d] = os.path.exists(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5') and len(dom_runs(d)[0]) >= 2
        except Exception as e: memo[d] = False; log('읽기 실패', d, repr(e)[:120])
    return memo[d]


def part_pick():
    fo = OUT + '/pick.json'
    if os.path.exists(fo): log('pick 있음'); return
    import pandas as pd
    p2 = open(ROOT + '/derived/r38/p2_67.txt').read().split(); s5 = open(ROOT + '/derived/r38/asamt_sanity5.txt').read().split()
    pool = sorted(d for d in p2 if d not in s5)
    e7 = sorted(np.random.default_rng(SEED['e7']).choice(pool, 5, replace=False).tolist())
    rest = [d for d in pool if d not in e7]
    tr = e7 + sorted(np.random.default_rng(SEED['trunc']).choice(rest, 5, replace=False).tolist())
    cand = pd.read_csv(RES + '/r37_trunk_candidates.tsv', sep='\t'); main = pd.read_csv(RES + '/r37_trunk_pairs.tsv', sep='\t')
    used = {frozenset((a, b)) for a, b in zip(main.a, main.b)}
    cand = cand[[frozenset((a, b)) not in used for a, b in zip(cand.a, cand.b)]]
    rng = np.random.default_rng(SEED['pairs']); pairs, skipped = [], []
    for lab, lo, hi, k in BINS:
        sub_ = cand[(cand.TM >= lo) & (cand.TM < hi)].sort_values(['a', 'b']).reset_index(drop=True); got = 0
        for o in rng.permutation(len(sub_)):
            r = sub_.iloc[int(o)]
            if avail(r.a) and avail(r.b):
                pairs.append(dict(a=r.a, b=r.b, n=int(r.n), m=int(r.m), TM=float(r.TM), level=str(r.level), bin=lab)); got += 1
            else: skipped.append(dict(a=r.a, b=r.b, bin=lab))
            _runs.clear()
            if got == k: break
        log('구간', lab, got, '쌍')
    reuse = {}
    for p in pairs:
        for d in (p['a'], p['b']): reuse[d] = reuse.get(d, 0) + 1
    jdump(dict(written=time.strftime('%F %T'), seeds=SEED, e7=e7, trunc=tr, pairs=pairs, skipped=skipped,
               domain_reuse={d: c for d, c in reuse.items() if c > 1},
               rule='E7 5 = P2 67 − 7c 바름 점검 5 에서 default_rng(20261020).choice(5); 잘라 낸 사슬 10 = E7 5 + 나머지에서 default_rng(20261022).choice(5); '
                    '30 쌍 = r37 후보 − 281 쌍, 구간마다 (a, b) 정렬 뒤 default_rng(20261021).permutation 순서로 두 도메인 모두 320 K 접힌 실행 ≥ 2 인 쌍을 5·5·5·5·5·3·2 개'), fo, bk=True)
    log('pick 끝', e7, tr[5:], len(pairs))


# ---------- b0 ----------
def part_b0():
    fo = RES + '/r38_s0_b0.json'
    if os.path.exists(fo): log('b0 있음'); return
    pk = json.load(open(OUT + '/pick.json'))['pairs']
    small = sorted(pk, key=lambda p: (max(p['n'], p['m']), p['a'], p['b']))[:10]; rows = []
    for p in small:
        x, y = medoid(p['a']), medoid(p['b']); n, m = len(x), len(y); tol = REL * 0.75 * max(n, m)
        t0 = time.time(); S_h = float(fs.align(x, y).score); th = time.time() - t0
        cert = fs.bnb(x, y, rel_tol=REL, device=DEV)
        t0 = time.time(); Rg = random_rotations(100000, seed=1, dtype=DT, device=DEV); Sg = []
        for c in range(0, len(Rg), 5000):
            K = torch.einsum('nxz,gmxz->gnm', x, torch.einsum('gxy,myz->gmxz', Rg[c:c + 5000], y)); Sg.append(dp_score(K))
        Sg = torch.cat(Sg); top = Rg[torch.topk(Sg, 256).indices]
        S_grid = max(float(Sg.max()), float(fs.align(x, y, starts=top).score)); tg = time.time() - t0
        S_best = max(cert['S_lb'], S_grid)
        r = dict(a=p['a'], b=p['b'], n=n, m=m, TM=p['TM'], tol=tol, S_heur=S_h, S_lb=cert['S_lb'], S_ub=cert['S_ub'], S_grid=S_grid,
                 certified=bool(cert['certified']), bnb_sec=cert['sec'], heur_sec=th, grid_sec=tg, levels=cert['levels'],
                 bnb_ge_heur=bool(cert['S_lb'] >= S_h - 1e-9), grid_le_ub=bool(S_grid <= cert['S_ub'] + 1e-9),
                 grid_vs_bnb_simlong=float(abs(S_grid - cert['S_lb']) / (0.75 * max(n, m))),
                 d_heur=float(np.sqrt(max(float(sq_norm(x) + sq_norm(y)) - 2 * S_h, 0))), d_best=float(np.sqrt(max(float(sq_norm(x) + sq_norm(y)) - 2 * S_best, 0))))
        r['pass'] = bool(r['certified'] and r['bnb_ge_heur'] and r['grid_le_ub'] and r['grid_vs_bnb_simlong'] <= REL)
        rows.append(r); log('b0', p['a'], p['b'], n, m, 'BnB', round(cert['sec'], 1), 's', '통과' if r['pass'] else '실패', round(r['grid_vs_bnb_simlong'], 5))
    jdump(dict(written=time.strftime('%F %T'), rows=rows, n_pass=sum(r['pass'] for r in rows), verdict='통과' if all(r['pass'] for r in rows) else '실패',
               bnb_sec_median=float(np.median([r['bnb_sec'] for r in rows])), rule='prereg 0.0: BnB S ≥ 128 출발 S, 격자 S ≤ BnB 상한, |격자 − BnB| ≤ 0.005 (sim_long), 인증'), fo, bk=True)


# ---------- e7 ----------
def part_e7():
    sel = json.load(open(OUT + '/pick.json'))['e7']; os.makedirs(OUT + '/e7', exist_ok=True)
    for d in sel:
        fo = f'{OUT}/e7/{d}.json'
        if os.path.exists(fo): continue
        t0 = time.time(); rng = seeded('r38_s0_e7', d); _runs.clear()
        ref = torch.as_tensor(np.load(f'{ROOT}/derived/r38/md/{d}.npz')['ref'])
        res = dict(dom=d, cells={}, nruns={}, diag={})
        for W_ in ('full', 'half'):
            R = {T: dom_runs(d, T, W_, TT)[0] for T in TT}; L = dom_runs(d, '320', W_, TT)[1]; ne = L - 1; res.update(L=L, ne=ne)
            for T in TT: res['nruns'][f'{T}|{W_}'] = len(R[T])
            Ts = [T for T in TT if len(R[T]) >= 2]
            if '320' not in Ts: continue
            E = {T: coords(np.concatenate([X for _, X in R[T]])) for T in Ts}
            A = {T: fs.align_to_reference(frames(np.concatenate([X for _, X in R[T]])), ref).numpy().reshape(len(R[T]) * N, -1) for T in Ts}
            pm = {}
            def go(Ta, Tb, rr):
                ms = MeanStart(E[Ta], E[Tb]); o = runpair_means(ms, rr, rng)
                ea = np.sqrt(((A[Ta][o['IA']] - A[Tb][o['IB']]) ** 2).sum(1))
                o['a2000'] = [float(ea[o['GR'] == g].mean()) for g in range(len(rr))]
                o['aall'] = [float(np.sqrt(cost(A[Ta][a * N:(a + 1) * N], A[Tb][b * N:(b + 1) * N])).mean()) for a, b in rr]
                res['diag'][f'{Ta}-{Tb}|{W_}'] = dict(widened=o['widened'], kmax=o['kmax'], n_pairs=len(o['dv']), S_mean=ms.S_mean)
                return {k: {rp: o[k][g] for g, rp in enumerate(rr)} for k in ('d', 'a2000', 'aall')}
            n320 = len(R['320']); pm[('320', '320')] = go('320', '320', list(itertools.combinations(range(n320), 2)))
            for T in Ts[1:]:
                nT = len(R[T]); pm[(T, T)] = go(T, T, list(itertools.combinations(range(nT), 2)))
                pm[('320', T)] = go('320', T, list(itertools.product(range(n320), range(nT))))
                cell = {}
                for v in ('d', 'a2000', 'aall'):
                    cell[{'d': 'D', 'a2000': 'A2000', 'aall': 'Aall'}[v]] = estats(np.mean(list(pm[('320', T)][v].values())), np.mean(list(pm[('320', '320')][v].values())), np.mean(list(pm[(T, T)][v].values())), ne)
                cell.update(nP=n320, nQ=nT); res['cells'][f'MD{T}-MD320|{W_}'] = cell
            if n320 == 5:
                for i, g2 in enumerate(itertools.combinations(range(5), 2)):
                    g3 = [k for k in range(5) if k not in g2]; cell = {}
                    for v in ('d', 'a2000', 'aall'):
                        f = lambda a, b: pm[('320', '320')][v][(min(a, b), max(a, b))]
                        cell[{'d': 'D', 'a2000': 'A2000', 'aall': 'Aall'}[v]] = estats(np.mean([f(a, b) for a in g2 for b in g3]), f(*g2), np.mean([f(a, b) for a, b in itertools.combinations(g3, 2)]), ne)
                    res['cells'][f'null320|{W_}|{i}'] = cell
        res['sec'] = time.time() - t0; jdump(res, fo); log('e7', d, len(res['cells']), '칸', round(res['sec']), 's')
    summ_e7()


def summ_e7():
    rows = []
    for f in sorted(glob.glob(OUT + '/e7/*.json')):
        j = json.load(open(f))
        for k, c in j['cells'].items(): rows.append(dict(dom=j['dom'], cell=k, kind='R0' if k.startswith('null') else 'T', win=k.split('|')[1], c=c))
    def cmp(sub_, ref):
        o = {}
        for q in ('D2', 'wbar', 'ratio'):
            x = np.array([r['c']['D'][q] for r in sub_]); y = np.array([r['c'][ref][q] for r in sub_]); ok = np.isfinite(x) & np.isfinite(y) & (np.abs(y) > 1e-12)
            o[q] = dict(spearman=float(spearmanr(x[ok], y[ok]).correlation) if ok.sum() > 2 else None, median_ratio=float(np.median(x[ok] / y[ok])) if ok.any() else None, n=int(ok.sum()))
        return o
    # 주 판: 온도 칸 + R0 는 도메인·창마다 10 나눔의 중앙값 한 칸
    r0 = {}
    for r in rows:
        if r['kind'] == 'R0': r0.setdefault((r['dom'], r['win']), []).append(r)
    r0rows = [dict(dom=k[0], cell=f'R0|{k[1]}', kind='R0med', win=k[1], c={v: {q: float(np.median([x['c'][v][q] for x in xs])) for q in ('D2', 'wbar', 'ratio')} for v in ('D', 'A2000', 'Aall')}) for k, xs in r0.items()]
    main = [r for r in rows if r['kind'] == 'T'] + r0rows
    out = dict(written=time.strftime('%F %T'), n_dom=len({r['dom'] for r in rows}), n_T=sum(r['kind'] == 'T' for r in rows), n_R0_split=sum(r['kind'] == 'R0' for r in rows),
               main_D_vs_Aall=cmp(main, 'Aall'), main_D_vs_A2000=cmp(main, 'A2000'), T_only_D_vs_Aall=cmp([r for r in rows if r['kind'] == 'T'], 'Aall'),
               all_splits_D_vs_Aall=cmp(rows, 'Aall'),
               rule='prereg 3 절·0.0b: Δ², w̄, 비 셋 모두 Spearman ≥ 0.9, D/맞춘 좌표 비의 중앙값 ∈ [0.8, 1.2]. 주 판 = 온도 칸 + R0(도메인·창마다 10 나눔 중앙값 한 칸), 기준 = 맞춘 좌표 전수(r38_compare 와 같은 통계)')
    m_ = out['main_D_vs_Aall']
    out['verdict'] = '통과' if all(m_[q]['spearman'] is not None and m_[q]['spearman'] >= 0.9 and 0.8 <= m_[q]['median_ratio'] <= 1.2 for q in m_) else '실패'
    jdump(out, RES + '/r38_s0_e7.json', bk=True); log('e7 요약', out['verdict'], {q: (round(m_[q]['spearman'] or 0, 3), round(m_[q]['median_ratio'] or 0, 3)) for q in m_})


# ---------- trunc ----------
def part_trunc():
    sel = json.load(open(OUT + '/pick.json'))['trunc']; os.makedirs(OUT + '/trunc', exist_ok=True)
    for d in sel:
        fo = f'{OUT}/trunc/{d}.json'
        if os.path.exists(fo): continue
        t0 = time.time(); rng = seeded('r38_s0_trunc', d); _runs.clear()
        R, L = dom_runs(d); ne = L - 1; nr = len(R)
        if nr < 2: jdump(dict(dom=d, skip='실행 < 2'), fo); continue
        Ef = coords(np.concatenate([X for _, X in R])); Et = Ef[:, :ne - CUT].contiguous()      # 같은 궤적의 좌표에서 끝 10 변만 뺀다
        cut = (Ef[:, ne - CUT:] ** 2).sum((-1, -2, -3)).cpu().numpy()                          # 프레임마다 0.75 Σs²(자른 변)
        ft = runpair_means(MeanStart(Ef, Et), [(a, b) for a in range(nr) for b in range(nr) if a != b], rng, want_gap=True)
        ff = runpair_means(MeanStart(Ef, Ef), list(itertools.combinations(range(nr), 2)), rng, want_gap=True)
        tt = runpair_means(MeanStart(Et, Et), list(itertools.combinations(range(nr), 2)), rng, want_gap=True)
        gap_diff = float(ft['gv'].mean() - ff['gv'].mean()); expect = float(cut[ft['IA']].mean())
        pair = estats(np.mean(ft['dpair']), np.mean(ff['dpair']), np.mean(tt['dpair']), ne)
        ref = torch.as_tensor(np.load(f'{ROOT}/derived/r38/md/{d}.npz')['ref']); null = []
        if nr == 5:
            Aa = fs.align_to_reference(frames(np.concatenate([X for _, X in R])), ref).numpy().reshape(nr * N, -1)
            pd_ = {(a, b): float(np.sqrt(cost(Aa[a * N:(a + 1) * N], Aa[b * N:(b + 1) * N])).mean()) for a, b in itertools.combinations(range(5), 2)}
            for g2 in itertools.combinations(range(5), 2):
                g3 = [k for k in range(5) if k not in g2]
                null.append(estats(np.mean([pd_[(min(a, b), max(a, b))] for a in g2 for b in g3]), pd_[g2], np.mean([pd_[x] for x in itertools.combinations(g3, 2)]), ne)['d2n'])
        res = dict(dom=d, L=L, ne=ne, nruns=nr, gap_diff=gap_diff, expect_075_sum_s2=expect, rel_err=float((gap_diff - expect) / expect),
                   gap_ff=float(ff['gv'].mean()), gap_ft=float(ft['gv'].mean()), pair=pair, null320_d2n=null,
                   in_domain_null=bool(null and min(null) <= pair['d2n'] <= max(null)) if null else None,
                   widened=dict(ft=ft['widened'], ff=ff['widened'], tt=tt['widened']), sec=time.time() - t0)
        jdump(res, fo); log('trunc', d, '상대 오차', round(res['rel_err'], 4), '짝 몫 d2n', round(pair['d2n'], 5))
    J = [json.load(open(f)) for f in sorted(glob.glob(OUT + '/trunc/*.json'))]; J = [j for j in J if 'skip' not in j]
    med = float(np.median([j['pair']['d2n'] for j in J])) if J else None
    out = dict(written=time.strftime('%F %T'), n=len(J), max_abs_rel_err=float(max(abs(j['rel_err']) for j in J)) if J else None,
               rows=[dict(dom=j['dom'], rel_err=j['rel_err'], pair_d2n=j['pair']['d2n'], in_domain_null=j['in_domain_null']) for j in J],
               pair_d2n_median=med, R0_band=R0_BAND, rule='prereg 0.1: |빈칸 몫 차 − 0.75 Σs²| / 0.75 Σs² < 2% (모든 도메인), 짝 몫 Δ²/√변 의 중앙값이 R0 띠 안(도메인 귀무 범위도 함께)')
    out['verdict_gap'] = '통과' if J and out['max_abs_rel_err'] < 0.02 else '실패'
    out['verdict_pair'] = '통과' if J and R0_BAND[0] <= med <= R0_BAND[1] else '실패'
    jdump(out, RES + '/r38_s0_trunc.json', bk=True); log('trunc 요약', out['verdict_gap'], out['verdict_pair'])


# ---------- ens ----------
def part_ens():
    pk = json.load(open(OUT + '/pick.json'))['pairs']; os.makedirs(OUT + '/ens', exist_ok=True)
    for p in pk:
        fo = f"{OUT}/ens/{p['a']}__{p['b']}.json"
        if os.path.exists(fo): continue
        t0 = time.time(); rng = seeded('r38_s0_ens', p['a'], p['b']); _runs.clear()
        EA, na, La = side(p['a']); EB, nb, Lb = side(p['b'])
        cr = runpair_means(MeanStart(EA, EB), list(itertools.product(range(na), range(nb))), rng, want_gap=True)
        fa = runpair_means(MeanStart(EA, EA), list(itertools.combinations(range(na), 2)), rng, want_gap=True)
        fb = runpair_means(MeanStart(EB, EB), list(itertools.combinations(range(nb), 2)), rng, want_gap=True)
        ne = ((La - 1) + (Lb - 1)) / 2
        res = dict(p, nA=na, nB=nb, ne_mean=ne, D=estats(np.mean(cr['d']), np.mean(fa['d']), np.mean(fb['d']), ne),
                   gap_cross=float(cr['gv'].mean()), d2_cross=float((cr['dv'] ** 2).mean()),
                   widened=dict(cross=cr['widened'], A=fa['widened'], B=fb['widened']), kmax=dict(cross=cr['kmax'], A=fa['kmax'], B=fb['kmax']), sec=time.time() - t0)
        jdump(res, fo); log('ens', p['bin'], p['a'], p['b'], 'Δ²_D', round(res['D']['D2'], 3), round(res['sec']), 's')
    J = [json.load(open(f)) for f in sorted(glob.glob(OUT + '/ens/*.json'))]
    jdump(dict(written=time.strftime('%F %T'), rows=[dict(a=j['a'], b=j['b'], bin=j['bin'], TM=j['TM'], ne_mean=j['ne_mean'], **j['D'], gap_share=j['gap_cross'] / j['d2_cross'] if j['d2_cross'] > 0 else None) for j in J]),
          RES + '/r38_s0_ens.json', bk=True)


# ---------- bias ----------
def part_bias():
    pk = json.load(open(OUT + '/pick.json'))['pairs']; fo = OUT + '/bias.jsonl'
    done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: j = json.loads(l); done.add((j['bin'], j['kind'], j['t']))
            except Exception: pass
    bypin = {lab: [p for p in pk if p['bin'] == lab] for lab, *_ in BINS}; msc = {}
    def ms_of(a, b):
        if (a, b) not in msc:
            if len(msc) > 6: msc.clear()
            msc[(a, b)] = MeanStart(side(a)[0], side(b)[0])
        return msc[(a, b)]
    for t in range(200):
        for lab, *_ in BINS:
            ps = bypin[lab]
            for kind in ('cross', 'floor'):
                if (lab, kind, t) in done or not ps: continue
                rng = seeded('r38_s0_bias', lab, kind, t); p = ps[t % len(ps)]
                if kind == 'cross':
                    a, b = p['a'], p['b']; na, nb = side(a)[1], side(b)[1]; i = int(rng.integers(na)) * N + int(rng.integers(N)); j = int(rng.integers(nb)) * N + int(rng.integers(N))
                else:
                    a = b = (p['a'] if (t // len(ps)) % 2 == 0 else p['b']); nr = side(a)[1]; ra, rb = rng.choice(nr, 2, replace=False)
                    i = int(ra) * N + int(rng.integers(N)); j = int(rb) * N + int(rng.integers(N))
                ms = ms_of(a, b); d_red, _, S_red, ku, _ = red_pairs(ms, [i], [j])
                x, y = ms.EA[i], ms.EB[j]; nx, ny = float(sq_norm(x)), float(sq_norm(y))
                cert = fs.bnb(x, y, rel_tol=REL, device=DEV); S_ref = max(cert['S_lb'], float(S_red[0]))
                d_ref = float(np.sqrt(max(nx + ny - 2 * S_ref, 0)))
                rec = dict(bin=lab, kind=kind, t=t, a=a, b=b, i=i, j=j, n=len(x), m=len(y), d_red=float(d_red[0]), d_ref=d_ref, err=float(d_red[0]) - d_ref,
                           S_red=float(S_red[0]), S_lb=cert['S_lb'], S_ub=cert['S_ub'], S_heur=cert['S_heur'], certified=bool(cert['certified']), tol=cert['tol'], k_used=int(ku[0]), bnb_sec=cert['sec'])
                with open(fo, 'a') as fh: fh.write(json.dumps(rec) + '\n')
                if t % 10 == 0 and kind == 'floor': log('bias', t, lab, 'BnB', round(cert['sec'], 1), 's')
    summ_bias()


def summ_bias():
    pk = json.load(open(OUT + '/pick.json'))['pairs']; J = [json.loads(l) for l in open(OUT + '/bias.jsonl')]
    ens = {(r['a'], r['b']): r for r in json.load(open(RES + '/r38_s0_ens.json'))['rows']}; rows = []
    for lab, *_ in BINS:
        c = [j for j in J if j['bin'] == lab and j['kind'] == 'cross']; f = [j for j in J if j['bin'] == lab and j['kind'] == 'floor']
        ps = [p for p in pk if p['bin'] == lab]
        if not ps: continue                      # 첫 세션 쌍이 없는 구간(0.8–1.01)
        D2 = float(np.median([ens[(p['a'], p['b'])]['D2'] for p in ps if (p['a'], p['b']) in ens])) if ps else float('nan')
        ne = float(np.mean([((p['n'] - 1) + (p['m'] - 1)) / 2 for p in ps])) if ps else float('nan')
        ec = float(np.mean([j['err'] for j in c])) if c else float('nan'); ef = float(np.mean([j['err'] for j in f])) if f else float('nan')
        thr = max(0.05 * abs(D2), 0.25 * (R0_BAND[1] - R0_BAND[0]) * np.sqrt(ne))
        rows.append(dict(bin=lab, n_cross=len(c), n_floor=len(f), err_cross=ec, err_floor=ef, bias=2 * (ec - ef), D2_bin=D2, ne_bin=ne, threshold=thr,
                         passed=bool(2 * (ec - ef) < thr), uncertified=sum(not j['certified'] for j in c + f), bnb_sec_median=float(np.median([j['bnb_sec'] for j in c + f])) if c + f else None))
    jdump(dict(written=time.strftime('%F %T'), rows=rows, verdict='통과' if all(r['passed'] for r in rows) else '표시할 구간 있음',
               rule='prereg 0.2: 2 × (교차 오차 평균 − 바닥 오차 평균) < max(구간 Δ²_D 의 5%, R0 띠 폭의 25%). 오차 = 줄인 계산 d − max(BnB 하한, 줄인 계산) 의 d. '
                    'R0 띠 폭 = 0.0039 (Δ²/√변) × √(구간 평균 변 수), 구간 Δ²_D = 그 구간 쌍들의 ens Δ²_D 중앙값. 첫 세션은 30 쌍의 7 구간만'), RES + '/r38_s0_bias.json', bk=True)


# ---------- ustat ----------
def part_ustat():
    fo = RES + '/r38_s0_ustat.json'
    if os.path.exists(fo): log('ustat 있음'); return
    pk = json.load(open(OUT + '/pick.json'))['pairs'][:10]; rows = []
    for p in pk:
        _runs.clear(); EA, na, _ = side(p['a']); EB, nb, _ = side(p['b']); ms = MeanStart(EA, EB); rng = seeded('r38_s0_ustat', p['a'], p['b'])
        ia, ib = np.meshgrid(np.arange(N), np.arange(N), indexing='ij'); d, *_ = red_pairs(ms, ia.ravel(), ib.ravel())
        sub_ = rng.choice(N * N, NPAIR, replace=False); ms_ = float(d[sub_].mean()); mf = float(d.mean())
        se = float(d[sub_].std(ddof=1) / np.sqrt(NPAIR) * np.sqrt(1 - NPAIR / (N * N)))
        rows.append(dict(a=p['a'], b=p['b'], mean_full=mf, mean_2000=ms_, se=se, z=abs(ms_ - mf) / se)); log('ustat', p['a'], p['b'], round(rows[-1]['z'], 2))
    zm = float(np.median([r['z'] for r in rows]))
    jdump(dict(written=time.strftime('%F %T'), rows=rows, median_z=zm, verdict='통과' if zm < 1 else '실패',
               rule='prereg 0.3: 10 쌍(30 쌍 앞 10, 첫 교차 실행 쌍)에서 |2,000 쌍 평균 − 전수 평균| / 표준오차(유한 모집단 보정)의 중앙값 < 1'), fo, bk=True)


# ---------- unc ----------
def part_unc():
    fo = RES + '/r38_s0_unc.json'
    if os.path.exists(fo): log('unc 있음'); return
    pk = json.load(open(OUT + '/pick.json'))['pairs']; rows = []
    for p in pk:
        _runs.clear(); x, y = medoid(p['a']), medoid(p['b']); n, m = len(x), len(y); nx, ny = (x ** 2).sum((-1, -2)), (y ** 2).sum((-1, -2))
        S, U, M = align_all(x, y); Smax = float(S.max()); near = torch.nonzero(S >= (1 - REL) * Smax).flatten().tolist()
        seen, vals = set(), []
        for k in near:
            key = M[k].cpu().numpy().tobytes()
            if key in seen: continue
            seen.add(key); gap = float((nx * (~M[k].any(1))).sum() + (ny * (~M[k].any(0))).sum()); d2 = float(nx.sum() + ny.sum() - 2 * S[k])
            vals.append(dict(S=float(S[k]), d2=d2, gap=gap, pair=d2 - gap))
        b = max(vals, key=lambda v: v['S']); ed = (n + m) / 2
        sp = lambda q: (max(v[q] for v in vals) - min(v[q] for v in vals))
        r = dict(a=p['a'], b=p['b'], bin=p['bin'], TM=p['TM'], n=n, m=m, n_near=len(near), n_distinct=len(vals),
                 pair_best=b['pair'], gap_best=b['gap'], pair_spread_per_edge=sp('pair') / ed, gap_spread_per_edge=sp('gap') / ed,
                 pair_rel=sp('pair') / b['pair'] if b['pair'] > 0 else None, gap_rel=sp('gap') / b['gap'] if b['gap'] > 0 else None, theta_star_spread='B5 에서')
        r['flag'] = bool((r['pair_rel'] or 0) > 0.2 or (r['gap_rel'] or 0) > 0.2); rows.append(r)
        log('unc', p['bin'], p['a'], p['b'], '서로 다른 대응', len(vals), '표시' if r['flag'] else '')
    jdump(dict(written=time.strftime('%F %T'), rows=rows, n_flag=sum(r['flag'] for r in rows),
               rule='prereg 0.4: 128 출발(항등 + 무작위 127, seed 0) 가운데 S ≥ (1 − 0.005) S_최적 인 서로 다른 대응들 사이 짝 몫·빈칸 몫의 최대 − 최소(변당 = /((n + m)/2)). 그 칸 값의 20% 를 넘으면 표시(서술). θ* 는 B5 에서'), fo, bk=True)


ALL = ['pick', 'b0', 'e7', 'trunc', 'ens', 'bias', 'ustat', 'unc']
for part in (ALL if PARTS == ['all'] else PARTS):
    log('부분 시작', part, DEV)
    try:
        globals()['part_' + part](); log('부분 끝', part)
    except Exception as e:
        traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
log('S0_DONE', ','.join(PARTS))
