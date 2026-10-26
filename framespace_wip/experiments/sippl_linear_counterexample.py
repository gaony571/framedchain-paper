"""Sippl 꼴(선형) 반례 (2026-10-17, 결과 문서 7.4A).

Sippl 2008 의 거리 D = L_a + L_b − 2S 에 S_TM 을 넣으면 2·D_TM² 이다. 이 꼴이 정확한 최적에서도 삼각 부등식을 어기는지 본다.
a: 나선 비슷한 60 잔기 사슬. b: a 의 k 곳을 r·d₀ 만큼 옮김. c: 같은 방향으로 2r·d₀ 만큼 옮김.
짝마다 δ²(2r) > 2δ²(r) (r² < 1/2) 이므로 선형 꼴은 깨지고, δ 는 거리라 제곱근 꼴은 선다.
ac 의 값은 최적화의 상계다. 항등 대응의 값과 0.3% 안이라 참 최적이 20% 낮을 여지는 보이지 않는다.

로컬 실행 (atlas 에서):
PYTHONPATH="framespace;framespace_wip" uv run --no-project --python 3.12 --with numpy --with numba --with torch --with scipy python framespace_wip/experiments/sippl_linear_counterexample.py

결과 (seed 0, d₀ = 4 Å):
k= 1 r=0.5: ab 0.195 bc 0.195 ac 0.497 | 선형 위반 +0.107 (+21.6%) | 제곱근 여유 −0.177
k= 5 r=0.5: ab 0.982 bc 0.982 ac 2.489 | +0.524 (+21.1%) | −0.405
k=10 r=0.5: ab 1.974 bc 1.975 ac 4.983 | +1.034 (+20.8%) | −0.578
k=10 r=0.3: ab 0.814 bc 0.814 ac 2.627 | +0.998 (+38.0%) | −0.184
k=10 r=0.7: ab 3.227 bc 3.227 ac 6.605 | +0.151 (+2.3%)  | −1.023
k=20 r=0.5: ab 3.929 bc 3.935 ac 9.947 | +2.083 (+20.9%) | −0.812
"""
import numpy as np

from stage0.rulers import tm_align, tm_cost

rng = np.random.default_rng(0)
n, d0 = 60, 4.0
t = np.arange(n)
A = np.c_[2.3 * np.cos(1.75 * t), 2.3 * np.sin(1.75 * t), 1.5 * t] + rng.normal(0, 0.3, (n, 3))
I, E, z = np.eye(n, dtype=bool), np.eye(3), np.zeros(3)
for k, r in [(1, 0.5), (5, 0.5), (10, 0.5), (10, 0.3), (10, 0.7), (20, 0.5)]:
    idx = np.linspace(3, n - 4, k).round().astype(int)
    u = rng.normal(size=(k, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    B = A.copy(); B[idx] += r * d0 * u
    C = A.copy(); C[idx] += 2 * r * d0 * u
    res = {nm: tm_align(P, Q, d0)['D2'] for nm, (P, Q) in dict(ab=(A, B), bc=(B, C), ac=(A, C)).items()}
    ident = {nm: tm_cost(P, Q, I, E, z, d0) for nm, (P, Q) in dict(ab=(A, B), bc=(B, C), ac=(A, C)).items()}
    lin = res['ac'] - (res['ab'] + res['bc'])
    sq = np.sqrt(res['ac']) - (np.sqrt(res['ab']) + np.sqrt(res['bc']))
    print(f"k={k:2d} r={r}: D2 ab={res['ab']:.3f} ({ident['ab']:.3f}) bc={res['bc']:.3f} ({ident['bc']:.3f}) "
          f"ac={res['ac']:.3f} ({ident['ac']:.3f}) | lin {lin:+.3f} ({lin / res['ac']:+.1%}) | sqrt {sq:+.3f}")
