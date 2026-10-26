"""탄성 모양 분석(ESA) 의 재구현: 비교 기준선 (연구 기획 E5, E7. 선행: Liu, Srivastava, Zhang 2010/2011; Srivastava, Klassen 2016).

단백질 백본을 원자 N, Cα, C 를 잇는 곡선 β 로 보고 길이를 1 로 맞춘다. 제곱근 속도 함수 q = β' / √|β'| 는 평행이동에 불변이고
‖q‖² = 길이 = 1 이라 단위 구 위에 있다. 두 곡선의 모양 거리는 회전 O 와 재매개화 γ (끝점 고정, 단조) 에 대한
d = min arccos ⟨q1, O (q2 ∘ γ) √γ'⟩ 이다. 회전은 SVD 로, γ 는 격자 위 동적계획법으로 번갈아 푼다.
이 공간과 다른 점: 끝점이 서로 맞붙어야 하고(부분 대응이 아님), 길이를 지워 모양만 보며, 프레임(펩타이드 평면의 방향)이 없다,
전역 최적의 보증이 없다(번갈아 푸는 국소해).
"""
import math

import numpy as np

NBRS = [(i, j) for i in range(1, 7) for j in range(1, 7) if math.gcd(i, j) == 1]


def backbone_curve(N, CA, C):
    """잔기마다의 N, Cα, C (L, 3) -> 곡선 점 (3L, 3) (N1, Cα1, C1, N2, ...)."""
    return np.stack([np.asarray(N, float), np.asarray(CA, float), np.asarray(C, float)], 1).reshape(-1, 3)


def resample(beta, T):
    """호 길이로 고르게 T 점을 다시 뽑고, 다시 뽑은 꺾은선의 길이를 1 로 맞춘다 (중심은 원점)."""
    beta = np.asarray(beta, float)
    seg = np.linalg.norm(np.diff(beta, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, s[-1], T)
    out = np.stack([np.interp(t, s, beta[:, k]) for k in range(3)], 1)
    out /= np.linalg.norm(np.diff(out, axis=0), axis=1).sum()      # 다시 뽑은 꺾은선의 길이를 1 로
    return out - out.mean(0)


def srvf(beta):
    """꺾은선 곡선 (T, 3) (매개 t ∈ [0, 1], 점 사이 간격 1/(T−1)) -> 선분마다의 q (T−1, 3).
    선분 위에서 β' 이 일정하므로 ‖q‖² = Σ |β'| Δt = 곡선의 길이가 정확히 나온다."""
    T = len(beta)
    v = np.diff(beta, axis=0) * (T - 1)
    sp = np.linalg.norm(v, axis=1)
    return v / np.sqrt(np.maximum(sp, 1e-12))[:, None]


def inner(q1, q2):
    """⟨q1, q2⟩ = ∫ q1·q2 dt (선분마다 일정)."""
    return float((q1 * q2).sum() / len(q1))


def best_rotation(q1, q2):
    """argmax_O ⟨q1, O q2⟩ (O ∈ SO(3))."""
    M = np.einsum('ti,tj->ij', q1, q2)                 # Σ q1 q2ᵀ
    U, _, Vt = np.linalg.svd(M)
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    return U @ D @ Vt


def dp_reparam(q1, q2):
    """격자 (T × T) 위의 동적계획법으로 ∫ |q1(t) − √γ'(t) q2(γ(t))|² dt 를 최소로 하는 γ (T,) (끝점 고정, 단조).
    마디 (k, l) → (i, j) 의 기울기는 (j − l)/(i − k), 이웃은 NBRS. 구간 적분은 그 사이의 격자점에서 사다리꼴로."""
    T = len(q1)
    t = np.linspace(0, 1, T)                                    # 선분 값을 이 격자의 표본으로 본다
    E = np.full((T, T), np.inf)
    E[0, 0] = 0.0
    back = np.zeros((T, T, 2), int)
    jj = np.arange(T)
    for i in range(1, T):
        for di, dj in NBRS:
            k = i - di
            if k < 0:
                continue
            l = jj - dj
            ok = l >= 0
            if not ok.any():
                continue
            js, ls = jj[ok], l[ok]
            m = dj / di
            ts = t[k:i + 1]                                     # 구간의 격자점 (di + 1,)
            g = t[ls][:, None] + m * (ts[None] - t[k])          # (nj, di+1)
            q2g = np.stack([np.interp(g, t, q2[:, c]) for c in range(3)], -1)
            r = ((q1[k:i + 1][None] - math.sqrt(m) * q2g) ** 2).sum(-1)
            cost = np.trapezoid(r, dx=1.0 / (T - 1), axis=1)
            cand = E[k, ls] + cost
            better = cand < E[i, js]
            E[i, js[better]] = cand[better]
            back[i, js[better]] = np.stack([np.full(better.sum(), k), ls[better]], 1)
    path = [(T - 1, T - 1)]
    while path[-1] != (0, 0):
        i, j = path[-1]
        path.append(tuple(back[i, j]))
    path = np.array(path[::-1])
    return np.interp(t, t[path[:, 0]], t[path[:, 1]])


def warp(beta, gamma):
    """곡선 (T, 3) 을 γ 로 다시 매개한다: β(γ(t))."""
    t = np.linspace(0, 1, len(beta))
    return np.stack([np.interp(gamma, t, beta[:, k]) for k in range(3)], 1)


def esa_distance(beta1, beta2, T=None, iters=3, reparam=True):
    """두 백본 곡선 (점 열) 의 모양 거리. T: 다시 뽑을 점 수(기본은 두 곡선 점 수의 큰 쪽, 400 이하로).
    반환 dict(d (구 위의 각), l2 (‖q1 − q̃2‖), inner, O, gamma)."""
    T = T or int(min(400, max(len(beta1), len(beta2))))
    b1, b2 = resample(beta1, T), resample(beta2, T)
    q1 = srvf(b1)
    O = np.eye(3)
    gamma = np.linspace(0, 1, T)
    for _ in range(iters):
        b2w = warp(b2, gamma) @ O.T
        q2 = srvf(b2w)
        O = best_rotation(q1, q2) @ O
        if reparam:
            q2 = srvf(b2 @ O.T)
            g = dp_reparam(q1, q2)
            gamma = np.interp(np.linspace(0, 1, T), np.linspace(0, 1, len(g)), g)
    q2 = srvf(warp(b2, gamma) @ O.T)
    O2 = best_rotation(q1, q2)
    q2 = q2 @ O2.T
    ip = inner(q1, q2)
    return dict(d=float(np.arccos(np.clip(ip, -1, 1))), l2=float(math.sqrt(max(inner(q1 - q2, q1 - q2), 0))), inner=ip, O=O2 @ O, gamma=gamma)
