"""평균과 퍼짐으로 나누기 (연구 기획 M2′, M2″ 방식 A, M4).

맞춘 좌표: 두 평균 점 사이의 대응 π̄ 와 회전 U0 를 하나 고정하고, 각 표본을 자기 평균에 회전으로 맞춘 뒤 (B 쪽은 U0 를 더 건다)
π̄ 의 열(짝지은 변, A 에만 있는 변, B 에만 있는 변)로 늘어놓는다. 없는 쪽은 영행렬이다. 이 좌표는 유클리드 공간 R^{L×9} 의 점이고,
그 위의 W₂ (W2_aligned) 에서는 W₂² = ‖μ_A − μ_B‖² + W₂²(중심화) 가 정확하며 최적 계획도 같다(교차항의 주변 합이 0).
열마다 나누면 c_col = (평균 몫) + (퍼짐 몫) 이고 Σ c_col = W₂² 이다.
표본 쌍마다 대응과 회전을 다시 최적화한 W₂ (w2.ensemble_w2) 는 W2_aligned 이하다. 둘을 함께 적는다.

짝 없는 변 i* 의 열은 다른 쪽이 늘 영이라 비용이 E‖x_i*‖² = 0.75 E[s_i*²] 로 정해져 있다(정적 비용). 그 평균 몫과 퍼짐 몫의 비율은
그 변 자신의 흩어짐이 정한다. 자리 찾기 판정에서는 짝 없는 열을 빼고 짝지은 열(동역학의 차이)만 쓴다.
평균 점이 수렴하면 맞춘 표본의 평균은 평균 점과 같다(info 의 converged 로 확인).
"""
import numpy as np
import torch

import framespace as fs

from .band import MeanStart
from .w2 import emd


def _cross(P, X, Y):
    """열마다 Σ_ij P_ij ⟨X_i, Y_j⟩ (L,). X (I, L, 9), Y (J, L, 9)."""
    T = P @ Y.reshape(len(Y), -1)
    return (X.reshape(len(X), -1) * T).sum(0).reshape(X.shape[1], -1).sum(1)


def aligned_coordinates(ms):
    """MeanStart -> 맞춘 좌표 XA (FA, L, 3, 3), XB (FB, L, 3, 3), 열 cols [(k 또는 −1, l 또는 −1)]."""
    A = torch.einsum('ixy,ikyz->ikxz', ms.RA, ms.EA)
    B = torch.einsum('xy,iyw,ikwz->ikxz', ms.U0, ms.RB, ms.EB)
    pairs = sorted(zip(*[t.tolist() for t in torch.nonzero(ms.match, as_tuple=True)]))
    n, m = A.shape[1], B.shape[1]
    cols, i, j = [], 0, 0
    for k, l in pairs:
        cols += [(q, -1) for q in range(i, k)] + [(-1, q) for q in range(j, l)] + [(k, l)]
        i, j = k + 1, l + 1
    cols += [(q, -1) for q in range(i, n)] + [(-1, q) for q in range(j, m)]
    L = len(cols)
    XA = A.new_zeros(len(A), L, 3, 3)
    XB = B.new_zeros(len(B), L, 3, 3)
    ka = [c for c, (k, _) in enumerate(cols) if k >= 0]
    kb = [c for c, (_, l) in enumerate(cols) if l >= 0]
    XA[:, ka] = A[:, [cols[c][0] for c in ka]]
    XB[:, kb] = B[:, [cols[c][1] for c in kb]]
    return XA, XB, cols


def decompose(EA, EB, wa=None, wb=None, mean_start=None, n_mean_starts=128, plan=None, seed=0):
    """M2′ 와 M2″ 방식 A. 반환 dict:
      W2sq_aligned, mean2 (= ‖μ_A − μ_B‖², 평균 점 사이의 d_fix²), spread2 (= W₂²(중심화)), identity_gap (W2sq − mean2 − spread2),
      varA, varB, var_ratio (각 앙상블의 퍼짐: 맞춘 좌표에서 자기 평균까지의 제곱 거리의 평균),
      cols, col_mean (L,), col_spread (L,), matched (L,) bool,
      static (짝 없는 열의 총비용), static_expected (0.75 Σ E[s²]), dynamic (짝지은 열의 총비용), dynamic_mean, dynamic_spread,
      edgeA (A 의 변마다의 몫 (n,)), edgeB (B 의 변마다 (m,)), plan, mean_start, d2_means (평균 점 사이의 d² 상계).
    plan 을 주면 그 계획으로 나눈다(기본은 맞춘 좌표의 최적 계획). 무게 wa, wb 는 수송의 주변 분포와 평균·퍼짐에 쓴다
    (평균 점 자체는 무게 없이 낸다: 무게 붙은 평균 점은 E6 에서 정한다)."""
    ms = mean_start or MeanStart(EA, EB, n_mean_starts=n_mean_starts, seed=seed)
    XA, XB, cols = aligned_coordinates(ms)
    FA, FB, L = len(XA), len(XB), len(cols)
    a = np.full(FA, 1.0 / FA) if wa is None else np.asarray(wa, float) / np.sum(wa)
    b = np.full(FB, 1.0 / FB) if wb is None else np.asarray(wb, float) / np.sum(wb)
    A2 = XA.reshape(FA, L, 9).double().cpu().numpy()
    B2 = XB.reshape(FB, L, 9).double().cpu().numpy()
    muA = np.einsum('i,ilc->lc', a, A2)
    muB = np.einsum('j,jlc->lc', b, B2)
    Ac, Bc = A2 - muA, B2 - muB
    if plan is None:
        C = (A2 ** 2).sum((1, 2))[:, None] + (B2 ** 2).sum((1, 2))[None] - 2 * np.einsum('ilc,jlc->ij', A2, B2)
        plan = emd(np.maximum(C, 0), None if wa is None else a, None if wb is None else b)
    P = np.asarray(plan)
    # 열마다: 평균 몫, 퍼짐 몫 Σ_ij P_ij ‖ã_i − b̃_j‖² = Σ a_i ‖ã_i‖² + Σ b_j ‖b̃_j‖² − 2 Σ P_ij ⟨ã_i, b̃_j⟩
    col_mean = ((muA - muB) ** 2).sum(1)
    col_spread = (np.einsum('i,ilc->l', a, Ac ** 2) + np.einsum('j,jlc->l', b, Bc ** 2)
                  - 2 * _cross(P, Ac, Bc))
    col_total = (np.einsum('i,ilc->l', a, A2 ** 2) + np.einsum('j,jlc->l', b, B2 ** 2)
                 - 2 * _cross(P, A2, B2))      # 열마다의 비용 Σ P_ij ‖a_i − b_j‖²
    W2sq = float(col_total.sum())
    mean2, spread2 = float(col_mean.sum()), float(col_spread.sum())
    matched = np.array([k >= 0 and l >= 0 for k, l in cols])
    varA = float(np.einsum('i,ilc->', a, Ac ** 2))
    varB = float(np.einsum('j,jlc->', b, Bc ** 2))
    edgeA = np.zeros(ms.EA.shape[1])
    edgeB = np.zeros(ms.EB.shape[1])
    for c, (k, l) in enumerate(cols):
        if k >= 0 and l >= 0:                                     # 짝지은 열은 반씩
            edgeA[k] += col_total[c] / 2
            edgeB[l] += col_total[c] / 2
        elif k >= 0:
            edgeA[k] += col_total[c]
        else:
            edgeB[l] += col_total[c]
    sA = (fs.as_tensor(ms.EA)[:, :, :, 0].pow(2).sum(-1) * 4).double().cpu().numpy()   # s² (FA, n)
    sB = (fs.as_tensor(ms.EB)[:, :, :, 0].pow(2).sum(-1) * 4).double().cpu().numpy()
    static_exp = 0.0
    for k, l in cols:
        if l < 0:
            static_exp += 0.75 * float(a @ sA[:, k])
        elif k < 0:
            static_exp += 0.75 * float(b @ sB[:, l])
    d2_means = float(fs.sq_norm(ms.mA) + fs.sq_norm(ms.mB) - 2 * ms.S_mean)
    return dict(W2sq_aligned=W2sq, W2_aligned=float(np.sqrt(max(W2sq, 0))), mean2=mean2, spread2=spread2,
                identity_gap=W2sq - mean2 - spread2, d2_means=d2_means,
                varA=varA, varB=varB, var_ratio=varA / varB if varB > 0 else float('inf'),
                cols=cols, col_mean=col_mean, col_spread=col_spread, col_total=col_total, matched=matched,
                static=float(col_total[~matched].sum()), static_expected=static_exp,
                dynamic=float(col_total[matched].sum()), dynamic_mean=float(col_mean[matched].sum()),
                dynamic_spread=float(col_spread[matched].sum()),
                edgeA=edgeA, edgeB=edgeB, plan=P, mean_start=ms,
                converged=(ms.infoA['converged'], ms.infoB['converged']))
