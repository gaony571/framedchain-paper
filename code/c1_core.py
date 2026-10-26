
"""EV·E0·사다리 CPU 핵심 (2026-10-18). 입력 derived/mdcath_5T (336 도메인, 5 온도 × 5 반복, 2 프레임 간격, 프레임별 Q).
도메인마다 derived/c1_core/<dom>.json (이어 하기 가능). 판정 기준은 연구 계획 4절 EV (평형 진단, C1-B, C1) 와 E0.
맞춘 좌표 = 320 K 접힌 프레임의 평균 점에 모든 온도의 프레임을 맞춘 것(대응 고정, 같은 서열). 지상 비용 = 맞춘 좌표의 제곱 유클리드(d_fix² 의 상계).
W₂ 는 POT 의 정확한 최적 수송(emd2). 중심 항 = 두 구름에서 각자 (무게) 평균을 뺀 W₂² (M2′)."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, time, glob
ROOT = sys.argv[1]; ONLY = sys.argv[2].split(',') if len(sys.argv) > 2 and sys.argv[2] != 'all' else None
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, torch, framespace as fs, ot
from scipy.spatial.transform import Rotation as Rot
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
torch.set_num_threads(1)
st = json.load(open(ROOT + '/results/mdcath_5T_setup.json')); STEP = st['frame_step']; E0SET, C1SET = set(st['e0_domains']), set(st['c1_domains'])
OUTD = ROOT + '/derived/c1_core'; os.makedirs(OUTD, exist_ok=True)
TEMPS = ['320', '348', '379', '413', '450']; QT = 0.75; Tn = {T: float(T) for T in TEMPS}
DS = ' HBEGITSP-'
def load(d):
    z = np.load(f'{ROOT}/derived/mdcath_5T/{d}.npz'); D = {}
    for T in TEMPS:
        reps = sorted({k.split('_')[1] for k in z.files if k.startswith(T + '_') and k.endswith('_X')}, key=int); Xs, Qs, Rp, Qf = [], [], [], []
        for r in reps:
            X = z[f'{T}_{r}_X']; q = z[f'{T}_{r}_Q'].astype(float); Xs.append(X); Qs.append(q[::STEP][:len(X)]); Rp.append(np.full(len(X), int(r))); Qf.append(q)
        D[T] = dict(X=np.concatenate(Xs), Q=np.concatenate(Qs), rep=np.concatenate(Rp), Qfull=Qf)
    ds = z['320_0_dssp']; dssp = np.array([DS[np.bincount(ds[:, k], minlength=10).argmax()] for k in range(ds.shape[1])])
    return D, dssp
def eqdiag(Qfull):
    out = []
    for q in Qfull:
        f = q >= QT; h = len(q) // 2; ref = 0; low = False
        for v in q:
            if v < 0.3: low = True
            elif v >= QT and low: ref += 1; low = False
        out.append(dict(frac=float(f.mean()), half_diff=float(abs(f[:h].mean() - f[h:].mean())), refold=int(ref)))
    fr = [o['frac'] for o in out]; mid = 0.1 <= np.mean(fr) <= 0.9
    ok = all(o['half_diff'] <= 0.10 for o in out) and (max(fr) - min(fr) <= 0.20) and ((not mid) or sum(o['refold'] for o in out) >= 1)
    return dict(mean_frac=round(float(np.mean(fr)), 4), range=round(float(max(fr) - min(fr)), 4), max_half_diff=round(max(o['half_diff'] for o in out), 4), refold=int(sum(o['refold'] for o in out)), passed=bool(ok))
def frames(X):
    return fs.coordinates_from_backbone(torch.as_tensor(X[:, :, 0]), torch.as_tensor(X[:, :, 1]), torch.as_tensor(X[:, :, 2]))
def cost(A, B):
    return np.maximum((A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2 * A @ B.T, 0)
def w2(A, B, wa=None, wb=None, center=False):
    wa = np.full(len(A), 1 / len(A)) if wa is None else wa / wa.sum(); wb = np.full(len(B), 1 / len(B)) if wb is None else wb / wb.sum()
    if center: A = A - wa @ A; B = B - wb @ B
    return float(ot.emd2(wa, wb, cost(A, B), numItermax=500000))
def dw2(A, B, rng, wa=None, wb=None, center=False):
    """편향을 뺀 W₂² (반쪽 교차): ½[W(A1,B2)+W(A2,B1)] − ½W(A1,A2) − ½W(B1,B2). 좁아진(줄어든) 구름에 벌점."""
    wa = np.ones(len(A)) if wa is None else np.asarray(wa, float); wb = np.ones(len(B)) if wb is None else np.asarray(wb, float)
    ia, ib = rng.permutation(len(A)), rng.permutation(len(B)); a1, a2 = ia[: len(A) // 2], ia[len(A) // 2:]; b1, b2 = ib[: len(B) // 2], ib[len(B) // 2:]
    W = lambda X, wx, Z, wz: w2(X, Z, wx, wz, center)
    return 0.5 * (W(A[a1], wa[a1], B[b2], wb[b2]) + W(A[a2], wa[a2], B[b1], wb[b1])) - 0.5 * W(A[a1], wa[a1], A[a2], wa[a2]) - 0.5 * W(B[b1], wb[b1], B[b2], wb[b2])
def edist(A, B, wa=None, wb=None):
    """에너지 거리 2E|a−b| − E|a−a'| − E|b−b'| (무게 허용, 자기 항은 대각 제외)."""
    wa = np.full(len(A), 1 / len(A)) if wa is None else wa / wa.sum(); wb = np.full(len(B), 1 / len(B)) if wb is None else wb / wb.sum()
    dab = np.sqrt(cost(A, B)); daa = np.sqrt(cost(A, A)); dbb = np.sqrt(cost(B, B))
    def selfterm(dm, w):
        W = np.outer(w, w); np.fill_diagonal(W, 0); return float((W * dm).sum() / W.sum())
    return float(2 * wa @ dab @ wb - selfterm(daa, wa) - selfterm(dbb, wb))
def sub(idx, n, rng):
    return idx if len(idx) <= n else rng.choice(idx, n, replace=False)
def polar(M):
    U, S, Vt = np.linalg.svd(M); d = np.sign(np.linalg.det(U @ Vt)); D = np.eye(3); D[2, 2] = d; return U @ D @ Vt
def lie_scale(Y, Rbar, c):
    """Y (F, ne, 3, 3) = (s/2) R. R̄ᵢ exp(c log(R̄ᵢᵀ Rᵢ)), 길이 비 그대로. c: 스칼라 또는 변마다 (ne,)."""
    F, ne = Y.shape[:2]; half = np.linalg.norm(Y, axis=(2, 3)) / np.sqrt(3); R = Y / half[..., None, None]
    rel = np.einsum('eji,fejk->feik', Rbar, R).reshape(-1, 3, 3)
    w = Rot.from_matrix(rel).as_rotvec().reshape(F, ne, 3); big = np.linalg.norm(w, axis=-1) > 2.8
    cc = np.broadcast_to(np.asarray(c, float), (ne,))[None, :, None] * np.ones((F, ne, 1)); cc[big] = 1.0
    Rc = np.einsum('eij,fejk->feik', Rbar, Rot.from_rotvec((w * cc).reshape(-1, 3)).as_matrix().reshape(F, ne, 3, 3))
    return Rc * half[..., None, None], float(big.mean())
def tot_spread(Yf):
    return float(((Yf - Yf.mean(0)) ** 2).sum(1).mean())
def solve_c(Y, Rbar, target, mask=None):
    lo, hi = 0.3, 4.0
    for _ in range(25):
        c = 0.5 * (lo + hi); cv = np.ones(Y.shape[1]) if mask is None else np.where(mask, c, 1.0)
        Yc, _ = lie_scale(Y, Rbar, c if mask is None else cv); s = tot_spread((Yc if mask is None else Yc[:, mask]).reshape(len(Y), -1))
        if s < target: lo = c
        else: hi = c
    return 0.5 * (lo + hi)
def feats(Ya, ref, cls, pca, Q, CA):
    dev = np.log(((Ya - ref[None]) ** 2).sum((-1, -2)) + 1e-6); K = cls.max() + 1
    cm = np.stack([dev[:, cls == k].mean(1) if (cls == k).any() else np.zeros(len(dev)) for k in range(K)], 1)
    rg = np.sqrt(((CA - CA.mean(1, keepdims=True)) ** 2).sum(-1).mean(1))
    return np.c_[cm, pca.transform(Ya.reshape(len(Ya), -1)), Q, rg]
def logit_cv(X, y, rep, Cs=(0.01, 0.1, 1.0)):
    """반복 단위 교차 적합 로짓 f(x) = log p(y=1)/p(y=0). 정칙화는 첫 겹에서 고름."""
    mu, sd = X.mean(0), X.std(0) + 1e-9; Xs = (X - mu) / sd; f = np.full(len(y), np.nan); best = None
    for r in np.unique(rep):
        tr, te = rep != r, rep == r
        if len(np.unique(y[tr])) < 2 or te.sum() == 0: continue
        if best is None:
            sc = []
            for C in Cs:
                m = LogisticRegression(C=C, max_iter=2000).fit(Xs[tr], y[tr]); p = np.clip(m.predict_proba(Xs[te])[:, 1], 1e-6, 1 - 1e-6)
                sc.append(-np.mean(y[te] * np.log(p) + (1 - y[te]) * np.log(1 - p)))
            best = Cs[int(np.argmin(sc))]
        m = LogisticRegression(C=best, max_iter=2000).fit(Xs[tr], y[tr]); f[te] = m.decision_function(Xs[te])
    return f, best
def multinom(X, Tlab, rep, Tlist, mode, iters=300):
    """교차 적합 다항 로지스틱. mode 'free' / 'beta' (로그 비가 β 에 선형: logit_T = β_T (w·x) + c_T) / 'betag' (+ g(T) 방향)."""
    mu, sd = X.mean(0), X.std(0) + 1e-9; Xs = torch.tensor((X - mu) / sd, dtype=torch.float64); y = torch.tensor(Tlab); K = len(Tlist); dim = Xs.shape[1]
    beta = torch.tensor([1 / Tn[T] for T in Tlist], dtype=torch.float64); beta = (beta - beta.mean()) / beta.std()
    T0 = 320.0; g = torch.tensor([np.log(Tn[T] / T0) - 1 + T0 / Tn[T] for T in Tlist], dtype=torch.float64); g = (g - g.mean()) / (g.std() + 1e-12)
    ll = []
    for r in np.unique(rep):
        tr, te = torch.tensor(rep != r), torch.tensor(rep == r)
        if te.sum() == 0: continue
        if mode == 'free': W = torch.zeros(dim, K, dtype=torch.float64, requires_grad=True); prm = [W]
        else: w1 = torch.zeros(dim, dtype=torch.float64, requires_grad=True); prm = [w1]
        if mode == 'betag': w2_ = torch.zeros(dim, dtype=torch.float64, requires_grad=True); prm.append(w2_)
        c = torch.zeros(K, dtype=torch.float64, requires_grad=True); prm.append(c)
        opt = torch.optim.LBFGS(prm, max_iter=iters, line_search_fn='strong_wolfe')
        def logits(Xb):
            if mode == 'free': return Xb @ W + c
            L = (Xb @ w1)[:, None] * beta[None, :] + c
            return L + (Xb @ w2_)[:, None] * g[None, :] if mode == 'betag' else L
        def closure():
            opt.zero_grad(); L = logits(Xs[tr]); loss = torch.nn.functional.cross_entropy(L, y[tr]) + 0.1 * sum((p_ ** 2).sum() for p_ in prm[:-1]) / dim; loss.backward(); return loss
        opt.step(closure)
        with torch.no_grad(): ll.append(float(torch.nn.functional.cross_entropy(logits(Xs[te]), y[te], reduction='sum')))
    return float(np.sum(ll) / len(Tlab))
def run(d):
    rng = np.random.default_rng(abs(hash(d)) % (2 ** 32)); t0 = time.time()
    D, dssp = load(d); res = dict(dom=d, eq={T: eqdiag(D[T]['Qfull']) for T in TEMPS})
    fold = {T: np.nonzero(D[T]['Q'] >= QT)[0] for T in TEMPS}; res['n_folded'] = {T: int(len(fold[T])) for T in TEMPS}
    if len(fold['320']) < 100: res['skip'] = 'few folded 320'; return res
    E = {T: frames(D[T]['X']) for T in TEMPS}
    ref = fs.mean_point(E['320'][fold['320']][:600]); Y = {T: fs.align_to_reference(E[T], ref).numpy() for T in TEMPS}; refn = ref.numpy(); ne = refn.shape[0]
    flat = {T: Y[T].reshape(len(Y[T]), -1) for T in TEMPS}
    # 클래스 (C3 와 같은 굵은 클래스, 묻힘은 도메인 안 삼분위)
    CA0 = D['320']['X'][fold['320'][0], :, 0]; M = np.linalg.norm(CA0[:, None] - CA0[None], axis=-1); cn = (M < 10).sum(1) - 1; eb = (cn[:-1] + cn[1:]) / 2
    s_ = np.array(['H' if c in 'HGI' else 'E' if c in 'EB' else 'L' for c in dssp]); ess = np.where(s_[:-1] == s_[1:], s_[:-1], 'L'); term = np.zeros(ne, int); term[:3] = 1; term[-3:] = 1
    cls = np.unique([f'{a}{b}{t}' for a, b, t in zip(ess, np.searchsorted(np.quantile(eb, [1 / 3, 2 / 3]), eb), term)], return_inverse=True)[1]
    # ---- 사다리(잡음·온도 행)와 E0 ----
    N = 100; A_, B_ = [0, 1], [2, 3, 4]
    f320A = fold['320'][np.isin(D['320']['rep'][fold['320']], A_)]; f320B = fold['320'][np.isin(D['320']['rep'][fold['320']], B_)]
    lad = {}
    if len(f320A) >= 50 and len(f320B) >= 50:
        a, b = sub(f320A, N, rng), sub(f320B, N, rng); lad['floor_320'] = dict(w2=w2(flat['320'][a], flat['320'][b]), w2c=w2(flat['320'][a], flat['320'][b], center=True), dw2=dw2(flat['320'][a], flat['320'][b], rng), ed=edist(flat['320'][a], flat['320'][b]))
        for T in TEMPS[1:]:
            if len(fold[T]) >= 50:
                c_ = sub(fold[T], N, rng); lad[f'320_{T}'] = dict(w2=w2(flat['320'][a], flat[T][c_]), w2c=w2(flat['320'][a], flat[T][c_], center=True), dw2=dw2(flat['320'][a], flat[T][c_], rng), ed=edist(flat['320'][a], flat[T][c_]), fracfold=round(len(fold[T]) / len(D[T]['Q']), 4))
    res['ladder'] = lad
    if d in E0SET:
        e0 = []
        for n in (25, 50, 100, 200):
            for k in range(5):
                if len(fold['348']) >= n and len(f320A) >= n and len(f320B) >= n:
                    a, b, c_ = sub(f320A, n, rng), sub(f320B, n, rng), sub(fold['348'], n, rng)
                    e0.append(dict(N=n, w2_temp=w2(flat['320'][a], flat['348'][c_]), w2_floor=w2(flat['320'][a], flat['320'][b]), dw2_temp=dw2(flat['320'][a], flat['348'][c_], rng), dw2_floor=dw2(flat['320'][a], flat['320'][b], rng),
                                   mean2_temp=float(((flat['320'][a].mean(0) - flat['348'][c_].mean(0)) ** 2).sum()), mean2_floor=float(((flat['320'][a].mean(0) - flat['320'][b].mean(0)) ** 2).sum())))
        res['e0'] = e0
    if d not in C1SET: res['sec'] = round(time.time() - t0, 1); return res
    # ---- C1-B ----
    Tb = [T for T in ('320', '348', '379', '413') if len(fold[T]) >= 100]; res['c1b_temps'] = Tb
    pca = PCA(10, random_state=0).fit(flat['320'][fold['320']])
    Fx = {T: feats(Y[T][fold[T]], refn, cls, pca, D[T]['Q'][fold[T]], D[T]['X'][fold[T], :, 0]) for T in Tb}; Rp = {T: D[T]['rep'][fold[T]] for T in Tb}
    c1b = {}
    if len(Tb) >= 3:
        def pair(T1, T2):
            X = np.r_[Fx[T1], Fx[T2]]; y = np.r_[np.zeros(len(Fx[T1])), np.ones(len(Fx[T2]))]; rp = np.r_[Rp[T1], Rp[T2]]; f, C = logit_cv(X, y, rp); return X, f
        # 가법성: 세 짝의 로짓을 모든 프레임(세 온도)에서
        T1, T2, T3 = Tb[:3]; Xall = np.r_[Fx[T1], Fx[T2], Fx[T3]]; rall = np.r_[Rp[T1], Rp[T2], Rp[T3]]
        def pair_on(Ta, Tb_):
            Xp = np.r_[Fx[Ta], Fx[Tb_]]; yp = np.r_[np.zeros(len(Fx[Ta])), np.ones(len(Fx[Tb_]))]; mu, sd = Xp.mean(0), Xp.std(0) + 1e-9
            out = np.full(len(Xall), np.nan)
            for r in np.unique(rall):
                trp = np.r_[Rp[Ta], Rp[Tb_]] != r
                m = LogisticRegression(C=0.1, max_iter=2000).fit((Xp[trp] - mu) / sd, yp[trp]); te = rall == r; out[te] = m.decision_function((Xall[te] - mu) / sd)
            return out
        f12, f23, f13 = pair_on(T1, T2), pair_on(T2, T3), pair_on(T1, T3); resid = f13 - (f12 + f23)
        # 바닥: 같은 온도(T1) 반복 묶음 A 대 B
        ya = np.isin(Rp[T1], A_).astype(float); half = np.zeros(len(ya), int)
        for r in np.unique(Rp[T1]):
            ix = np.nonzero(Rp[T1] == r)[0]; half[ix[len(ix) // 2:]] = 1
        ffl, _ = logit_cv(Fx[T1], ya, half)   # 시간 반쪽 교차 적합(반복 하나를 빼면 그 반복의 묶음 이름을 못 봐 정의가 깨진다)
        floor_rms = float(np.sqrt(np.nanmean((ffl - np.nanmean(ffl)) ** 2)))
        c1b['additivity'] = dict(temps=[T1, T2, T3], resid_rms=round(float(np.sqrt(np.nanmean((resid - np.nanmean(resid)) ** 2))), 4), f13_rms=round(float(np.sqrt(np.nanmean((f13 - np.nanmean(f13)) ** 2))), 4),
                                 floor_rms=round(floor_rms, 4), corr_f13_vs_sum=round(float(np.corrcoef(f13[np.isfinite(f13)], (f12 + f23)[np.isfinite(f13)])[0, 1]), 4))
        c1b['additivity']['pass'] = bool(c1b['additivity']['resid_rms'] <= 1.5 * floor_rms)
        # 다항: 자유 대 온도 기저
        Xm = np.concatenate([Fx[T] for T in Tb]); ym = np.concatenate([np.full(len(Fx[T]), k) for k, T in enumerate(Tb)]); rm = np.concatenate([Rp[T] for T in Tb])
        llf, llb, llbg = multinom(Xm, ym, rm, Tb, 'free'), multinom(Xm, ym, rm, Tb, 'beta'), (multinom(Xm, ym, rm, Tb, 'betag') if len(Tb) >= 3 else np.nan)
        # 바닥: 같은 온도 반복 묶음 분류의 로그 손실 이득(우연 ln2 대비)
        p = 1 / (1 + np.exp(-np.clip(ffl, -30, 30))); yy = ya; flg = float(np.log(2) - np.nanmean(-(yy * np.log(p + 1e-9) + (1 - yy) * np.log(1 - p + 1e-9))))
        c1b['multinomial'] = dict(temps=Tb, ll_free=round(llf, 4), ll_beta=round(llb, 4), ll_betag=round(llbg, 4) if np.isfinite(llbg) else None, floor_gain=round(flg, 4),
                                  pass_beta=bool((llb - llf) <= 1.5 * max(flg, 1e-3)), pass_betag=bool(np.isfinite(llbg) and (llbg - llf) <= 1.5 * max(flg, 1e-3)))
    c1b['eq_pass_temps'] = [T for T in Tb if res['eq'][T]['passed']]
    res['c1b'] = c1b
    # ---- C1 경로 비교 (F 판) ----
    c1 = {}
    for (Ta, Tm, Tz) in (('320', '348', '379'), ('348', '379', '413')):
        fa, fm, fz = fold[Ta], fold[Tm], fold[Tz]
        if min(len(fa), len(fz)) < 100 or len(fm) < 60: continue
        tA = np.isin(D[Tm]['rep'][fm], [0, 1, 2]); truth, tB = fm[tA], fm[~tA]
        if len(truth) < 30 or len(tB) < 30: continue
        n = min(150, len(fa), len(fz)); ia, iz = sub(fa, n, rng), sub(fz, n, rng); YA, YZ = flat[Ta][ia], flat[Tz][iz]
        t = (1 / Tn[Ta] - 1 / Tn[Tm]) / (1 / Tn[Ta] - 1 / Tn[Tz]); TR = flat[Tm][sub(truth, 150, rng)]; TB = flat[Tm][sub(tB, 150, rng)]
        P = {}
        G = ot.emd(np.full(n, 1 / n), np.full(n, 1 / n), cost(YA, YZ)); jz = G.argmax(1); P['P1_transport'] = (YA + t * (YZ[jz] - YA), None)
        # P2 볼츠만 재가중, P3 혼합: 양 끝 프레임(반복 교차 적합 로짓)
        FA = feats(Y[Ta][ia], refn, cls, pca, D[Ta]['Q'][ia], D[Ta]['X'][ia, :, 0]); FZ = feats(Y[Tz][iz], refn, cls, pca, D[Tz]['Q'][iz], D[Tz]['X'][iz, :, 0])
        fpool, _ = logit_cv(np.r_[FA, FZ], np.r_[np.zeros(n), np.ones(n)], np.r_[D[Ta]['rep'][ia], D[Tz]['rep'][iz]])
        fpool = np.nan_to_num(fpool); wB = np.exp(t * fpool - np.logaddexp(0, fpool)); pool = np.r_[YA, YZ]
        P['P2_boltzmann'] = (pool, wB); P['P3_mixture'] = (pool, np.r_[np.full(n, 1 - t), np.full(n, t)]); P['P4_identity'] = (YA, None)
        # P6 계열: 아래 끝 구름의 리 대수 배수
        YA4 = Y[Ta][ia]; Rbar = np.stack([polar(YA4[:, e].mean(0)) for e in range(ne)])
        sA, sZ = tot_spread(YA), tot_spread(YZ); target = sA + t * (sZ - sA)
        Yh, bigf = lie_scale(YA4, Rbar, np.sqrt(Tn[Tm] / Tn[Ta])); P['P6_harmonic'] = (Yh.reshape(n, -1), None)
        c6 = solve_c(YA4, Rbar, target); Yp, _ = lie_scale(YA4, Rbar, c6); P['P6p_isotropic'] = (Yp.reshape(n, -1), None)
        cvec = np.ones(ne); YZ4 = Y[Tz][iz]
        for k in np.unique(cls):
            mk = cls == k; sAk = tot_spread(YA4[:, mk].reshape(n, -1)); sZk = tot_spread(YZ4[:, mk].reshape(n, -1)); cvec[mk] = solve_c(YA4, Rbar, sAk + t * (sZk - sAk), mask=mk)
        Ypp, _ = lie_scale(YA4, Rbar, cvec); P['P6pp_class'] = (Ypp.reshape(n, -1), None)
        out = dict(t=round(t, 4), n=int(n), c_isotropic=round(c6, 4), c_harmonic=round(float(np.sqrt(Tn[Tm] / Tn[Ta])), 4), frac_rot_near_pi=round(bigf, 4),
                   floor=dict(w2=w2(TB, TR), w2c=w2(TB, TR, center=True), dw2=dw2(TB, TR, rng), dw2c=dw2(TB, TR, rng, center=True), ed=edist(TB, TR)))
        for k, (Yp_, wp) in P.items(): out[k] = dict(w2=w2(Yp_, TR, wa=wp), w2c=w2(Yp_, TR, wa=wp, center=True), dw2=dw2(Yp_, TR, rng, wa=wp), dw2c=dw2(Yp_, TR, rng, wa=wp, center=True), ed=edist(Yp_, TR, wa=wp),
                                               spread=float(((Yp_ - (np.full(len(Yp_), 1 / len(Yp_)) if wp is None else wp / wp.sum()) @ Yp_) ** 2).sum(1) @ (np.full(len(Yp_), 1 / len(Yp_)) if wp is None else wp / wp.sum())))
        out['truth_spread'] = float(((TR - TR.mean(0)) ** 2).sum(1).mean())
        c1[f'{Ta}_{Tm}_{Tz}'] = out
    res['c1'] = c1; res['sec'] = round(time.time() - t0, 1)
    return res
doms = ONLY if ONLY else st['domains']
for d in doms:
    p = f'{OUTD}/{d}.json'
    if os.path.exists(p) and not ONLY: continue
    try:
        r = run(d); json.dump(r, open(p, 'w'))
        print(time.strftime('%H:%M:%S'), d, r.get('sec'), flush=True)
    except Exception as e:
        import traceback; print('ERR', d, repr(e)[:300], traceback.format_exc()[-800:], flush=True)
if not ONLY: open(ROOT + '/results/logs/c1_core.done', 'w').write(time.strftime('%F %T'))
