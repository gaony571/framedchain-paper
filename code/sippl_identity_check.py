"""(2026-11-02 추가, 기록 없음) S2 Text 8.3 의 비교 점검: 선형 꼴 반례의 (a, c) 값을 항등 대응의 값과 견준다.
반례 스크립트 framespace_wip/experiments/sippl_linear_counterexample.py (2026-10-17, 실행본 그대로)의 머리 주석
"항등 대응의 값과 0.3% 안" 은 항등 대응에 항등 겹침(U = I, t = 0; 그 스크립트의 tm_cost(P, Q, I, E, z, d0))을 놓은 값과의 비교다.
S2 Text 8.3 은 항등 대응에 최적 겹침을 놓은 값과 견준다. 이 파일은 같은 사슬(seed 0, 같은 난수 순서)과 같은 함수(stage0.rulers)로
최적화 값, 항등 대응·항등 겹침 값, 항등 대응·최적 겹침 값을 함께 낸다.
최적 겹침: 항등 대응(잔기 i ↔ i)을 고정하고 가중 최소제곱(무게 1/(1 + r²/d₀²)², rulers._kabsch)을 200 번 되풀이한다.
실행 (저장소 뿌리에서):
  uv run --offline --no-project --python 3.12 --with numpy --with numba --with torch --with scipy python code/sippl_identity_check.py
출력 (2026-11-02, Python 3.12, numpy 2.5.4, numba 0.68.1):
k=1 r=0.5: ac opt 0.4966 | ident pose 0.5000 (+0.69%) | ident corr + opt pose 0.4966 (-0.000%)
k=5 r=0.5: ac opt 2.4888 | ident pose 2.5000 (+0.45%) | ident corr + opt pose 2.4888 (+0.000%)
k=10 r=0.5: ac opt 4.9831 | ident pose 5.0000 (+0.34%) | ident corr + opt pose 4.9831 (+0.000%)
k=10 r=0.3: ac opt 2.6268 | ident pose 2.6471 (+0.77%) | ident corr + opt pose 2.6268 (+0.000%)
k=10 r=0.7: ac opt 6.6052 | ident pose 6.6216 (+0.25%) | ident corr + opt pose 6.6052 (+0.000%)
k=20 r=0.5: ac opt 9.9469 | ident pose 10.0000 (+0.53%) | ident corr + opt pose 9.9469 (-0.000%)
읽기: 최적화 값은 항등 대응·최적 겹침 값과 0.001% 안에서 같고, 항등 겹침 값보다는 0.25–0.77% 낮다.
"""
import os
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(REPO, 'framespace_pkg'), os.path.join(REPO, 'framespace_wip')]
from stage0.rulers import tm_align, tm_cost, _kabsch  # noqa: E402

# 반례 스크립트와 같은 사슬: 나선 비슷한 60 잔기, seed 0
rng = np.random.default_rng(0)
n, d0 = 60, 4.0
t = np.arange(n)
A = np.c_[2.3 * np.cos(1.75 * t), 2.3 * np.sin(1.75 * t), 1.5 * t] + rng.normal(0, 0.3, (n, 3))
I, E, z = np.eye(n, dtype=bool), np.eye(3), np.zeros(3)
ii = np.arange(n)


def opt_pose(P, Q):
    """항등 대응을 고정한 최적 겹침: 무게 1/(1 + r²/d₀²)² 의 가중 최소제곱을 200 번 되풀이(값은 줄지 않음)."""
    U, tt = E.copy(), z.copy()
    for _ in range(200):
        Y = Q @ U.T + tt
        r2 = ((P - Y) ** 2).sum(1) / d0 ** 2
        U, tt = _kabsch(P, Q, ii, ii, 1.0 / (1.0 + r2) ** 2)
    return U, tt


for k, r in [(1, 0.5), (5, 0.5), (10, 0.5), (10, 0.3), (10, 0.7), (20, 0.5)]:
    # c 는 반례 스크립트와 같은 난수 순서로 만든다(방향 u 를 같은 순서로 뽑음; b 는 난수를 쓰지 않고 여기서는 필요 없음)
    idx = np.linspace(3, n - 4, k).round().astype(int)
    u = rng.normal(size=(k, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    C = A.copy(); C[idx] += 2 * r * d0 * u
    res = tm_align(A, C, d0)['D2']                    # 최적화 값 (반례 표의 ac)
    idp = tm_cost(A, C, I, E, z, d0)                  # 항등 대응 + 항등 겹침 (반례 스크립트의 ident)
    U, tt = opt_pose(A, C)
    ido = tm_cost(A, C, I, U, tt, d0)                 # 항등 대응 + 최적 겹침 (S2 Text 8.3)
    print(f"k={k} r={r}: ac opt {res:.4f} | ident pose {idp:.4f} ({(idp-res)/res:+.2%}) | ident corr + opt pose {ido:.4f} ({(ido-res)/res:+.3%})")
