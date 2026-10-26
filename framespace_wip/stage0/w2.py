"""앙상블 사이의 W₂ (연구 기획 M2, M3, E0).

앙상블은 표본들의 (무게 붙은) 경험 분포다. 지상 비용은 표본 쌍마다 이 공간의 d² (정규화하지 않은 값)이고,
W₂² = min over 수송 계획 Σ P_kl d²(A_k, B_l). 표본 쌍의 d 를 상계로 계산하므로 낸 W₂ 도 참값의 상계다.
길이가 다른 쌍을 한 표에 놓을 때는 W₂ / √(0.75 (n + m)) 눈금(scaled)을 쓴다.

지상 비용을 내는 방식 셋:
  full: 쌍마다 항등 + 무작위 출발 회전(기본 128개)과 전체 DP. 검증용 기준.
  mean: 평균 점 출발(band.MeanStart). 출발 회전 1–4개, 전체 DP.
  band: 평균 점 출발 + π̄ 를 중심선으로 한 띠 DP. 띠 가장자리에 짝이 닿은 쌍은 반폭을 두 배로 넓혀 다시 푼다.
"""
import math
import time

import numpy as np
import torch

import framespace as fs
from framespace.distance import sq_norm
from framespace.transport import optimal_plan as emd   # 정확한 최적 수송은 패키지(0.2.0)로 옮겼다

from .band import MeanStart, align_band, band_windows


# ---------- 최적 수송 ----------

def w2_from_costs(C, wa=None, wb=None):
    """지상 비용 행렬 (d²) -> dict(W2, W2sq, plan)."""
    C = np.asarray(C, dtype=np.float64)
    P = emd(C, wa, wb)
    w2sq = float(max((P * C).sum(), 0.0))
    return dict(W2=math.sqrt(w2sq), W2sq=w2sq, plan=P)


def scaled(W2, n, m):
    """길이가 다른 쌍의 눈금 W₂ / √(0.75 (n + m)) (n, m 은 영 아닌 항의 수)."""
    return W2 / math.sqrt(0.75 * (n + m))


# ---------- 지상 비용 ----------

def _pairs(I, J, device):
    ii, jj = torch.meshgrid(torch.arange(I, device=device), torch.arange(J, device=device), indexing='ij')
    return ii.flatten(), jj.flatten()


def pair_costs(EA, EB, method='band', n_starts=128, n_local=1, k=5, k_max=None, iters=10, chunk=512, seed=0,
               keep=False, mean_start=None, device=None, dtype=None):
    """모든 표본 쌍의 d² (I, J).
    method: 'full' (출발 n_starts 개, 전체 DP), 'mean' (평균 점 출발 n_local 개, 전체 DP), 'band' (평균 점 출발 + 띠 DP).
    k: 띠의 처음 반폭. k_max: 넓힐 수 있는 최대 반폭(기본은 열 수, 그 반폭이면 전체 DP 와 같다).
    keep: 짝 (I*J, n, m) 과 회전 (I*J, 3, 3) 도 돌려준다 (쌍 번호 p = i * J + j).
    반환 dict(cost, sec, method, n_widened (띠를 넓힌 쌍 수), k_used (쌍마다 마지막 반폭), mean_start, [match, rotation])."""
    t0 = time.time()
    EA = fs.as_tensor(EA, device=device, dtype=dtype)
    EB = fs.as_tensor(EB, device=EA.device, dtype=EA.dtype)
    I, J, n, m = len(EA), len(EB), EA.shape[1], EB.shape[1]
    ii, jj = _pairs(I, J, EA.device)
    S = EA.new_zeros(I * J)
    M = torch.zeros(I * J, n, m, dtype=torch.bool, device=EA.device) if keep else None
    UU = EA.new_zeros(I * J, 3, 3) if keep else None
    kused = np.zeros(I * J, int)
    ms = None
    if method in ('mean', 'band'):
        ms = mean_start or MeanStart(EA, EB, n_mean_starts=n_starts, seed=seed)
    k_max = m if k_max is None else k_max
    widened = 0
    for c0 in range(0, I * J, chunk):
        sl = slice(c0, min(c0 + chunk, I * J))
        a, b = EA[ii[sl]], EB[jj[sl]]
        if method == 'full':
            al = fs.align(a, b, n_random=n_starts - 1, iters=iters, seed=seed + c0)
            s_, U_, M_ = al.score, al.rotation, al.match
            kused[sl] = m
        else:
            st = ms.starts(ii[sl], jj[sl], n_local=n_local, seed=seed + c0)
            if method == 'mean':
                al = fs.align(a, b, starts=st, iters=iters)
                s_, U_, M_ = al.score, al.rotation, al.match
                kused[sl] = m
            else:
                kk = k
                todo = torch.arange(len(a), device=EA.device)
                s_ = a.new_zeros(len(a)); U_ = a.new_zeros(len(a), 3, 3)
                M_ = torch.zeros(len(a), n, m, dtype=torch.bool, device=EA.device)
                kc = np.zeros(len(a), int)
                while len(todo):
                    lo, w = band_windows(ms.center, m, kk)
                    sc, Ur, Mr, edge = align_band(a[todo], b[todo], lo, w, st[todo], iters=iters)
                    s_[todo], U_[todo], M_[todo] = sc, Ur, Mr
                    kc[todo.cpu().numpy()] = kk
                    if w >= m or kk >= k_max:
                        break
                    todo = todo[edge]
                    widened += int(len(todo)) if kk == k else 0
                    kk *= 2
                kused[sl] = kc
        S[sl] = s_
        if keep:
            M[sl], UU[sl] = M_, U_
    C = (sq_norm(EA)[ii] + sq_norm(EB)[jj] - 2 * S).clamp_min(0).view(I, J)
    out = dict(cost=C.double().cpu().numpy(), sec=time.time() - t0, method=method, n_widened=widened, k_used=kused.reshape(I, J), mean_start=ms)
    if keep:
        out.update(match=M, rotation=UU)
    return out


def ensemble_w2(EA, EB, wa=None, wb=None, **kw):
    """앙상블 사이 W₂. 반환 dict(W2, W2sq, plan, cost, sec, ...). kw 는 pair_costs 로."""
    pc = pair_costs(EA, EB, **kw)
    res = w2_from_costs(pc['cost'], wa, wb)
    res.update({k: v for k, v in pc.items() if k != 'cost'}, cost=pc['cost'])
    return res


def soft_correspondence(plan, match, I, J):
    """방식 B (M2″): 수송 계획의 무게로 쌍마다의 짝을 평균한 확률 대응 P (n, m), A 의 변마다 짝지어질 확률 p (n,),
    행마다의 엔트로피 (n,) (짝이 미끄러지는 정도)."""
    w = torch.tensor(np.asarray(plan).ravel(), dtype=torch.float64, device=match.device)
    Pc = torch.einsum('q,qnm->nm', w, match.double()).cpu().numpy()
    p = Pc.sum(1)
    q = np.clip(Pc / np.maximum(p[:, None], 1e-12), 1e-300, 1)
    H = -(np.where(Pc > 0, q * np.log(q), 0.0)).sum(1)
    return Pc, p, H


# ---------- 바닥값, 편향, 안정성 ----------

def split_by_groups(groups, first):
    g = np.asarray(groups)
    return np.isin(g, first), ~np.isin(g, first)


def replica_floors(E, groups, n_per=None, seed=0, **kw):
    """독립 replica 쌍마다의 W₂ (M3). groups (F,) 는 replica 번호. n_per 를 주면 replica 마다 고르게 n_per 개.
    반환: [(g1, g2, W2)]."""
    E = fs.as_tensor(E)
    g = np.asarray(groups)
    labs = np.unique(g)
    rng = np.random.default_rng(seed)
    pick = {}
    for v in labs:
        idx = np.where(g == v)[0]
        if n_per is not None and len(idx) > n_per:
            idx = idx[np.linspace(0, len(idx) - 1, n_per).astype(int)]
        pick[v] = idx
    out = []
    for a in range(len(labs)):
        for b in range(a + 1, len(labs)):
            r = ensemble_w2(E[torch.as_tensor(pick[labs[a]])], E[torch.as_tensor(pick[labs[b]])], **kw)
            out.append((labs[a], labs[b], r['W2']))
    return out


def half_split_floor(E, **kw):
    """replica 가 없을 때의 바닥값: 궤적 앞 반 대 뒤 반 (표시해서 쓴다)."""
    E = fs.as_tensor(E)
    h = len(E) // 2
    return ensemble_w2(E[:h], E[h:2 * h], **kw)['W2']


def fit_bias(Ns, w2sq):
    """W₂² = a + b / N 을 최소제곱으로 맞춘다. 반환 (a, b)."""
    Ns, y = np.asarray(Ns, dtype=np.float64), np.asarray(w2sq, dtype=np.float64)
    X = np.stack([np.ones_like(Ns), 1.0 / Ns], 1)
    (a, b), *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(a), float(b)


def bootstrap_w2(C, B=200, seed=0):
    """지상 비용 행렬에서 표본을 다시 뽑아(행과 열을 복원 추출) W₂ 를 다시 낸다. 반환 dict(W2s, mean, sd, cv)."""
    C = np.asarray(C)
    I, J = C.shape
    rng = np.random.default_rng(seed)
    v = []
    for _ in range(B):
        r, c = rng.integers(0, I, I), rng.integers(0, J, J)
        v.append(w2_from_costs(C[np.ix_(r, c)])['W2'])
    v = np.array(v)
    return dict(W2s=v, mean=float(v.mean()), sd=float(v.std(ddof=1)), cv=float(v.std(ddof=1) / v.mean()))


def triangle_violation(w_ab, w_bc, w_ac):
    """세 앙상블의 W₂ 에서 삼각 부등식 위반의 크기 (세 방향 중 가장 큰 것, 위반이 없으면 0)."""
    return max(0.0, w_ac - w_ab - w_bc, w_ab - w_ac - w_bc, w_bc - w_ab - w_ac)
