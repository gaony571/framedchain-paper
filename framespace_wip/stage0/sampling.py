"""앙상블 표본 (연구 기획 M1): 프레임 사이 거리로 잰 자기상관과 그 간격으로 뽑기.

한 궤적 x_t 에서 g(τ) = E_t d_fix²(x_t, x_{t+τ}) 는 τ 가 커지면 2·Var 로 간다 (Var 는 평균 점까지의 제곱 거리의 평균).
자기상관 C(τ) = 1 − g(τ) / g_∞ 이고, 적분 자기상관 시간 τ_int = 1/2 + Σ_{τ≥1} C(τ) 를 C 가 처음 0 이하가 되는 곳에서 끊어 낸다.
독립에 가까운 표본의 간격은 2 τ_int 이다. 같은 단백질 안이라 대응은 변 번호 그대로다(d_fix, SVD 한 번으로 정확).
"""
import math

import numpy as np
import torch

import framespace as fs
from framespace.distance import sq_norm


def paired_fixed_sq(x, y):
    """쌍 (x_p, y_p) 마다의 d_fix² (P,). x, y (P, n, 3, 3), 대응은 변 번호 그대로."""
    x, y = fs.as_tensor(x), fs.as_tensor(y).to(fs.as_tensor(x))
    H = torch.einsum('pkcb,pkab->pca', y, x)
    sv = torch.linalg.svdvals(H)
    S = sv[..., 0] + sv[..., 1] + torch.sign(torch.det(H)) * sv[..., 2]
    return (sq_norm(x) + sq_norm(y) - 2 * S).clamp_min(0)


def distance_autocorr(E, max_lag=None, lags=None, n_ref=2000, seed=0):
    """한 궤적 E (T, n, 3, 3) 의 거리 자기상관. 반환 dict(lags, g, g_inf, C, tau_int, spacing).
    g_inf 는 서로 멀리 떨어진 무작위 프레임 쌍(시간 간격 ≥ T/4)의 d_fix² 평균. lags 를 주지 않으면 1..max_lag (기본 T/4)."""
    E = fs.as_tensor(E)
    T = len(E)
    max_lag = max_lag or max(1, T // 4)
    lags = np.arange(1, max_lag + 1) if lags is None else np.asarray(lags)
    g = np.array([float(paired_fixed_sq(E[:T - L], E[L:]).mean()) for L in lags])
    rng = np.random.default_rng(seed)
    i = rng.integers(0, T, n_ref)
    j = rng.integers(0, T, n_ref)
    far = np.abs(i - j) >= T // 4
    g_inf = float(paired_fixed_sq(E[torch.as_tensor(i[far])], E[torch.as_tensor(j[far])]).mean()) if far.any() else float(g[-1])
    C = 1 - g / g_inf
    cut = np.argmax(C <= 0) if (C <= 0).any() else len(C)
    tau = 0.5 + float(C[:cut].sum())
    return dict(lags=lags, g=g, g_inf=g_inf, C=C, tau_int=tau, spacing=int(math.ceil(2 * tau)), truncated=cut == len(C))


def spaced_indices(T, spacing, N=None, skip=0):
    """길이 T 의 궤적에서 skip 프레임을 버리고 spacing 간격으로 고른 번호. N 을 주면 앞에서부터 N 개까지."""
    idx = np.arange(skip, T, max(1, int(spacing)))
    return idx if N is None else idx[:N]


def pooled_sample(E, groups, N, spacing=None, skip=0.0, seed=0):
    """replica 여러 개 (groups) 에서 모두 N 개를 고르게 뽑는다. spacing 을 주면 replica 안에서 그 간격 이상으로.
    채우지 못하면 있는 만큼 돌려주고 유효 표본 수(n_eff)를 함께 낸다."""
    g = np.asarray(groups)
    labs = np.unique(g)
    per = int(math.ceil(N / len(labs)))
    out = []
    for v in labs:
        idx = np.where(g == v)[0]
        idx = idx[int(len(idx) * skip):]
        if spacing:
            idx = idx[::max(1, int(spacing))]
        if len(idx) > per:
            idx = idx[np.linspace(0, len(idx) - 1, per).astype(int)]
        out.append(idx)
    idx = np.concatenate(out)
    if len(idx) > N:
        idx = np.sort(np.random.default_rng(seed).choice(idx, N, replace=False))
    return idx, len(idx)
