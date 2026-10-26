"""(round41 g6, 2026-10-21 UTC; Drive JSON round41_2026-10-21) 일반 보조정리의 단백질 밖 예 (보충, NCS 판에만).
가족 다섯, 세 짝 (x, y, z) 1,000 개씩 (seed 20261087):
  (i)   ℝ 값 열, 군 없음, ℓ2 판: D² = min Σ_짝 (x_i − y_j)² + Σ_짝없음 a²  (c = |a|; 이 구성)          — DP 가 정확
  (ii)  같은 열의 ℓ1 판 = ERP (g = 0): min Σ_짝 |x_i − y_j| + Σ_짝없음 |a|  ("ERP 는 p = 1 판")    — DP 가 정확
  (iii) 평면 곡선(변 벡터 ∈ ℝ², SO(2) 몫, c = ‖a‖) ℓ2 판: D² = ‖x‖² + ‖y‖² − 2 max_θ S(θ), S(θ) = 순서를 지키는 부분 대응의 최대 Σ ⟨x_i, R_θ y_j⟩.
        θ 는 립시츠 분지한정(상수 √(Σ‖x‖² Σ‖y‖²), 상대 허용 1e-7)으로 인증 → D 의 구간 [D_lo, D_hi].
  음성 대조 (iv) (i) 에서 c = 2|a| (1-립시츠 깨짐), (v) DTW (ℓ2, 빈칸 없음).
판정: 삼각 부등식 셋(세 짝마다 세 방향). 정확한 가족은 위반 = d(x,z) > d(x,y) + d(y,z) + 1e-9·(1 + 우변). (iii) 은 확실한 위반 = D_lo(x,z) > D_hi(x,y) + D_hi(y,z),
  구간이 겹쳐 가를 수 없으면 "판정 불가" 로 센다. 짝·빈칸 분해 몫(빈칸 몫 / D²)과 편집 수를 함께 낸다. 점검: 회전한 같은 곡선 D = 0, 영 항을 끼운 열 D = 0.
사용: python r41_synth.py <ROOT>   출력: results/r41_synth.json, results/figs/r41_synth.png"""
import os, sys, json, time, shutil
import numpy as np
ROOT = sys.argv[1]; RES = ROOT + '/results'; FO = RES + '/r41_synth.json'; NT = 1000; SEED = 20261087; TOL = 1e-7


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def jdump(o, p):
    json.dump(o, open(p + '.tmp', 'w'), ensure_ascii=False, indent=1, default=float); os.replace(p + '.tmp', p)
    b = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + os.path.basename(p))


# ---------- 정확한 DP (ℝ 값) ----------
def dp_min(x, y, sub, gap):
    """min Σ_짝 sub(x_i, y_j) + Σ_짝없음 gap(·) 와 그 정렬의 (짝 비용, 빈칸 비용, 짝 수)."""
    n, m = len(x), len(y); F = np.zeros((n + 1, m + 1)); B = np.zeros((n + 1, m + 1), int)
    gx, gy = np.array([gap(a) for a in x]), np.array([gap(b) for b in y])
    F[1:, 0] = np.cumsum(gx); F[0, 1:] = np.cumsum(gy); B[1:, 0] = 1; B[0, 1:] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            c = (F[i - 1, j - 1] + sub(x[i - 1], y[j - 1]), F[i - 1, j] + gx[i - 1], F[i, j - 1] + gy[j - 1]); k = int(np.argmin(c)); F[i, j] = c[k]; B[i, j] = k
    i, j, pc, gc, nm = n, m, 0.0, 0.0, 0
    while i > 0 or j > 0:
        k = B[i, j]
        if k == 0: pc += sub(x[i - 1], y[j - 1]); nm += 1; i, j = i - 1, j - 1
        elif k == 1: gc += gx[i - 1]; i -= 1
        else: gc += gy[j - 1]; j -= 1
    return F[n, m], pc, gc, nm


def dtw(x, y):
    n, m = len(x), len(y); F = np.full((n + 1, m + 1), np.inf); F[0, 0] = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1): F[i, j] = (x[i - 1] - y[j - 1]) ** 2 + min(F[i - 1, j - 1], F[i - 1, j], F[i, j - 1])
    return np.sqrt(F[n, m])


FAM = dict(
    i=lambda x, y: (lambda r: (np.sqrt(r[0]), r))(dp_min(x, y, lambda a, b: (a - b) ** 2, lambda a: a * a)),
    ii=lambda x, y: (lambda r: (r[0], r))(dp_min(x, y, lambda a, b: abs(a - b), lambda a: abs(a))),
    iv=lambda x, y: (lambda r: (np.sqrt(r[0]), r))(dp_min(x, y, lambda a, b: (a - b) ** 2, lambda a: 4 * a * a)),
    v=lambda x, y: (dtw(x, y), None))


# ---------- 평면 곡선 (SO(2) 몫), 인증 ----------
def dp_max(K):
    """K (T, n, m) -> (T,) 순서를 지키는 부분 대응의 최대 Σ K (framespace.distance.dp_score 와 같은 점화식)."""
    F = np.zeros(K.shape[:1] + (K.shape[2] + 1,))
    for i in range(K.shape[1]):
        G = np.maximum(F[:, 1:], F[:, :-1] + K[:, i, :]); F = np.concatenate([F[:, :1], np.maximum.accumulate(G, 1)], 1)
    return F[:, -1]


def S_theta(AB, th):
    A, B = AB; return dp_max(np.cos(th)[:, None, None] * A[None] + np.sin(th)[:, None, None] * B[None])


def max_S(x, y, tol=TOL, n0=256, cap=400000):
    A = x @ y.T; B = x[:, 1:2] * y[:, 0][None] - x[:, 0:1] * y[:, 1][None]; Lc = np.sqrt((x ** 2).sum() * (y ** 2).sum())
    h = 2 * np.pi / n0; c = (np.arange(n0) + 0.5) * h; v = S_theta((A, B), c); LB = v.max(); tb = c[v.argmax()]; hs = np.full(n0, h); nev = n0; capped = False
    while True:
        keep = v + Lc * hs / 2 > LB + tol * max(LB, 1e-12)
        if not keep.any(): break
        c, hs = c[keep], hs[keep]
        if len(c) > cap: capped = True; break
        c = np.concatenate([c - hs / 4, c + hs / 4]); hs = np.concatenate([hs / 2, hs / 2]); v = np.concatenate([S_theta((A, B), c[i:i + 4096]) for i in range(0, len(c), 4096)]); nev += len(c)
        if v.max() > LB: LB = v.max(); tb = c[v.argmax()]
    UB = max(LB + tol * max(LB, 1e-12), float((v + Lc * hs / 2).max()) if capped else 0.0)
    return LB, UB, tb, nev, capped


def plane_D(x, y):
    LB, UB, tb, nev, capped = max_S(x, y); nx, ny = (x ** 2).sum(), (y ** 2).sum()
    lo, hi = np.sqrt(max(nx + ny - 2 * UB, 0)), np.sqrt(max(nx + ny - 2 * LB, 0))
    K = np.cos(tb) * (x @ y.T) + np.sin(tb) * (x[:, 1:2] * y[:, 0][None] - x[:, 0:1] * y[:, 1][None])
    n, m = K.shape; F = np.zeros((n + 1, m + 1)); Bk = np.zeros((n + 1, m + 1), int)       # 그 θ 에서 짝 되짚기
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            c = (F[i - 1, j - 1] + K[i - 1, j - 1], F[i - 1, j], F[i, j - 1]); k = int(np.argmax(c)); F[i, j] = c[k]; Bk[i, j] = k
    i, j, ua, ub = n, m, np.ones(n, bool), np.ones(m, bool)
    while i > 0 and j > 0:
        k = Bk[i, j]
        if k == 0: ua[i - 1] = ub[j - 1] = False; i, j = i - 1, j - 1
        elif k == 1: i -= 1
        else: j -= 1
    gc = float((x[ua] ** 2).sum() + (y[ub] ** 2).sum())
    return lo, hi, dict(gap=gc, D2=hi ** 2, nev=nev, capped=capped)


# ---------- 자료 ----------
def edit(rng, x, newval):
    """기준 열에서 무작위 편집: 지우기·끼우기 확률 U(0, 0.3), 섭동 SD U(0.1, 1.0)."""
    pd, pi, sd = rng.uniform(0, .3), rng.uniform(0, .3), rng.uniform(.1, 1.); out = []; ne = 0
    for a in x:
        if rng.random() < pi: out.append(newval()); ne += 1
        if rng.random() < pd: ne += 1; continue
        out.append(a + sd * rng.normal(size=np.shape(a)))
    if not out: out.append(newval())
    return np.array(out), ne


def series(rng):
    n = rng.integers(15, 41); x = np.zeros(n); x[0] = rng.normal() * 1.4
    for k in range(1, n): x[k] = 0.7 * x[k - 1] + rng.normal()
    return x


def curve(rng):
    n = rng.integers(15, 41); ang = np.cumsum(rng.normal(0, 0.6, n)); ln = 1 + 0.1 * rng.normal(size=n)
    return np.stack([ln * np.cos(ang), ln * np.sin(ang)], 1)


def rot(v, t):
    c, s = np.cos(t), np.sin(t); return v @ np.array([[c, s], [-s, c]])


def main():
    rng = np.random.default_rng(SEED); out = dict(written=None, n_triples=NT, families={}); t0 = time.time()
    trip_s = []
    for k in range(NT):
        x = series(rng); y, ey = edit(rng, x, lambda: rng.normal() * 1.4); src = y if k % 2 else x; z, ez = edit(rng, src, lambda: rng.normal() * 1.4); trip_s.append((x, y, z, ey, ez))
    trip_c = []
    for k in range(NT):
        x = curve(rng); y, ey = edit(rng, x, lambda: rot(np.array([1.0, 0.0]), rng.uniform(0, 2 * np.pi))); src = y if k % 2 else x; z, ez = edit(rng, src, lambda: rot(np.array([1.0, 0.0]), rng.uniform(0, 2 * np.pi)))
        trip_c.append((x, rot(y, rng.uniform(0, 2 * np.pi)), rot(z, rng.uniform(0, 2 * np.pi)), ey, ez))
    for fam in ('i', 'ii', 'iv', 'v'):
        f = FAM[fam]; viol = 0; slack = []; gsh = []; ned = []
        for x, y, z, ey, ez in trip_s:
            dxy, rxy = f(x, y); dyz, _ = f(y, z); dxz, _ = f(x, z)
            for a, b, c in ((dxz, dxy, dyz), (dxy, dxz, dyz), (dyz, dxy, dxz)):
                s = b + c - a; slack.append(s / max(b + c, 1e-12))
                if a > b + c + 1e-9 * (1 + b + c): viol += 1
            if rxy is not None and rxy[0] > 0: gsh.append(rxy[2] / rxy[0]); ned.append(ey)
        o = dict(n_ineq=3 * NT, violations=int(viol), slack_rel_min=float(np.min(slack)), slack_rel_q01=float(np.percentile(slack, 1)))
        if gsh:
            from scipy.stats import spearmanr
            o.update(gap_share_median=float(np.median(gsh)), gap_share_iqr=[float(np.percentile(gsh, 25)), float(np.percentile(gsh, 75))], spearman_gapshare_edits=float(spearmanr(gsh, ned).correlation))
        out['families'][fam] = o; out['_slack_' + fam] = slack; log('가족', fam, o)
    viol = und = 0; slack = []; gsh = []; ned = []; nev = []; capn = 0
    for k, (x, y, z, ey, ez) in enumerate(trip_c):
        Dxy = plane_D(x, y); Dyz = plane_D(y, z); Dxz = plane_D(x, z)
        for (alo, ahi, _), (blo, bhi, _), (clo, chi, _) in ((Dxz, Dxy, Dyz), (Dxy, Dxz, Dyz), (Dyz, Dxy, Dxz)):
            if alo > bhi + chi: viol += 1
            elif ahi > blo + clo: und += 1
            slack.append((blo + clo - ahi) / max(bhi + chi, 1e-12))
        info = Dxy[2]; nev += [Dxy[2]['nev'], Dyz[2]['nev'], Dxz[2]['nev']]; capn += sum(D[2]['capped'] for D in (Dxy, Dyz, Dxz))
        if info['D2'] > 0: gsh.append(info['gap'] / info['D2']); ned.append(ey)
        if k % 200 == 0: log('평면', k, round(time.time() - t0), 's')
    from scipy.stats import spearmanr
    out['families']['iii'] = dict(n_ineq=3 * NT, violations_certain=int(viol), undecided=int(und), slack_rel_min=float(np.min(slack)), slack_rel_q01=float(np.percentile(slack, 1)),
                                  gap_share_median=float(np.median(gsh)), gap_share_iqr=[float(np.percentile(gsh, 25)), float(np.percentile(gsh, 75))], spearman_gapshare_edits=float(spearmanr(gsh, ned).correlation),
                                  bnb_evals_median=float(np.median(nev)), bnb_capped=int(capn), tol_rel=TOL)
    out['_slack_iii'] = slack; log('가족 iii', out['families']['iii'])
    # 점검: 회전한 같은 곡선, 영 항 끼운 열
    x = curve(np.random.default_rng(SEED + 1)); lo, hi, _ = plane_D(x, rot(x, 1.234)); xs = series(np.random.default_rng(SEED + 2)); xz = np.insert(xs, [3, 7, 7], 0.0)
    out['checks'] = dict(plane_rotated_self_D_hi=float(hi), series_zero_inserted_D=float(FAM['i'](xs, xz)[0]), erp_zero_inserted_D=float(FAM['ii'](xs, xz)[0]))
    out['prediction'] = dict(i_ii_iii_zero=bool(out['families']['i']['violations'] == 0 and out['families']['ii']['violations'] == 0 and out['families']['iii']['violations_certain'] == 0),
                             controls_positive={f: bool(out['families'][f]['violations'] > 0) for f in ('iv', 'v')})
    out['written'] = time.strftime('%F %T'); out['sec'] = round(time.time() - t0, 1)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
        for f, lab in (('i', '(i) ℓ2, c=|a|'), ('ii', '(ii) ℓ1 = ERP'), ('iii', '(iii) plane/SO(2)'), ('iv', '(iv) c=2|a|'), ('v', '(v) DTW')):
            s = np.array(out['_slack_' + f]); ax[0].hist(np.clip(s, -0.5, 1), bins=60, histtype='step', label=f"{lab}: {out['families'][f].get('violations', out['families'][f].get('violations_certain'))}/{3 * NT}")
        ax[0].axvline(0, color='k', lw=0.8); ax[0].set_xlabel('triangle slack / (d(x,y)+d(y,z))'); ax[0].set_ylabel('count'); ax[0].legend(fontsize=7)
        ax[1].text(0.02, 0.5, '\n'.join(f"{f}: gap share median {out['families'][f]['gap_share_median']:.2f}, ρ(edits) {out['families'][f]['spearman_gapshare_edits']:.2f}" for f in ('i', 'ii', 'iii')), fontsize=8); ax[1].axis('off')
        os.makedirs(RES + '/figs', exist_ok=True); fig.tight_layout(); fig.savefig(RES + '/figs/r41_synth.png', dpi=150)
    except Exception as e:
        out['fig_err'] = repr(e)[:200]
    for f in ('i', 'ii', 'iii', 'iv', 'v'): out.pop('_slack_' + f, None)
    out['rule'] = 'round41 g6 (Drive JSON round41_2026-10-21): 예측 (i)–(iii) 위반 0/3,000 (iii 는 인증 구간으로 확실한 위반), (iv)·(v) 는 > 0 예상(0 이어도 반증 아님).'
    jdump(out, FO); log('G6_DONE', out['prediction'], out['checks'])


if __name__ == '__main__':
    main()
