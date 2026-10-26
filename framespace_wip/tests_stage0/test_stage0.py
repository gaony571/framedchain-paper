"""0단계 코드 시험. 합성 프레임 열 앙상블로 잰다 (pytest)."""
import math

import numpy as np
import pytest
import torch

import framespace as fs
import stage0 as s0


def chain(n, seed, noise=0.15):
    """이웃 변이 비슷한 회전을 이루는 실제 구조 (n, 3, 3) = s R / 2."""
    g = torch.Generator().manual_seed(seed)
    step = fs.expmap(torch.tensor([0.3, 1.2, 0.4], dtype=torch.float64))
    R = [torch.eye(3, dtype=torch.float64)]
    for _ in range(n - 1):
        R.append(R[-1] @ step @ fs.expmap(noise * torch.randn(3, generator=g, dtype=torch.float64)))
    s = 1 + 0.016 * torch.randn(n, generator=g, dtype=torch.float64)
    return s[:, None, None] * torch.stack(R) / 2


def ensemble(base, F, seed, sigma=0.08, rotate=True):
    """base (n, 3, 3) 의 변마다 작은 회전 요동을 준 표본 F 개. rotate 면 표본마다 전체 회전도 건다(같은 점)."""
    g = torch.Generator().manual_seed(seed)
    n = len(base)
    W = fs.expmap(sigma * torch.randn(F, n, 3, generator=g, dtype=torch.float64))
    E = torch.einsum('nxy,fnyz->fnxz', base, W)
    if rotate:
        E = torch.einsum('fxy,fnyz->fnxz', fs.random_rotations(F, seed=seed + 7), E)
    return E


def insert_edge(base, p, seed=0):
    """변 p 앞에 변 하나를 끼운 구조 (n+1, 3, 3) (끼운 변은 이웃 프레임을 조금 돌린 것)."""
    g = torch.Generator().manual_seed(seed)
    new = base[p] @ fs.expmap(0.5 * torch.randn(3, generator=g, dtype=torch.float64))
    return torch.cat([base[:p], new[None], base[p:]], 0)


# ---------- 띠 DP ----------

def test_band_equals_masked_full_dp():
    """띠 DP 는 띠 밖의 짝을 막은 전체 DP 와 정확히 같다 (중심선이 단조일 때)."""
    rng = np.random.default_rng(0)
    for trial in range(20):
        n, m, k = rng.integers(5, 25), rng.integers(5, 25), int(rng.integers(1, 4))
        K = torch.tensor(rng.normal(size=(3, n, m)))
        center = np.sort(rng.integers(0, m, n))
        lo, w = s0.band_windows(center, m, k)
        col = lo[:, None] + np.arange(w)[None]
        Kb = K[:, np.arange(n)[:, None], col]
        mask = torch.full((n, m), -1e9, dtype=torch.float64)
        mask[np.arange(n)[:, None], col] = 0
        ref = fs.dp_score(K + mask)
        assert torch.allclose(s0.dp_score_band(Kb, lo), ref, atol=1e-9), trial


def test_band_full_width_equals_dp():
    rng = np.random.default_rng(1)
    K = torch.tensor(rng.normal(size=(4, 12, 9)))
    lo, w = s0.band_windows(np.zeros(12, int), 9, 20)
    assert w == 9 and torch.allclose(s0.dp_score_band(K, lo), fs.dp_score(K), atol=1e-12)


def test_align_band_full_width_matches_align():
    a, b = chain(20, 1), chain(23, 2)
    st = fs.random_rotations(8, seed=3)[None]
    full = fs.align(a[None], b[None], starts=st)
    lo, w = s0.band_windows(np.zeros(20, int), 23, 30)
    S, U, M, edge = s0.align_band(a[None], b[None], lo, w, st)
    assert abs(float(S[0] - full.score[0])) < 1e-9 and bool((M[0] == full.match[0]).all()) and not bool(edge[0])


# ---------- 지상 비용과 W₂ ----------

@pytest.fixture(scope='module')
def two_ensembles():
    base = chain(30, 11)
    A = ensemble(base, 12, 1)
    B = ensemble(insert_edge(base, 15, 4), 12, 2)
    return base, A, B


def test_pair_cost_methods_agree(two_ensembles):
    """줄인 계산(평균 점 출발, 띠)이 시작점 128개의 전체 계산과 같은 값을 낸다 (쉬운 합성 자료)."""
    _, A, B = two_ensembles
    full = s0.pair_costs(A, B, method='full')['cost']
    mean = s0.pair_costs(A, B, method='mean', n_local=1)['cost']
    band = s0.pair_costs(A, B, method='band', n_local=1, k=3)
    # 한 출발점의 국소해라 드물게 조금 나쁠 수 있다 (2026-10-16: 144쌍 중 1쌍이 상대 2.5e-6 높음). 크기는 E0 에서 잰다.
    assert np.abs(mean - full).max() < 1e-4 * full.max()
    assert np.abs(band['cost'] - full).max() < 1e-4 * full.max()
    assert (band['cost'] >= full - 1e-9).all()        # 띠는 더 좁힌 문제라 같거나 크다(상계)


def test_band_widens_when_needed():
    """중심선이 반폭만큼 어긋나면 최적 경로가 띠 가장자리에 닿고, 넓혀서 바로잡는다."""
    base = chain(25, 5)
    A = ensemble(base, 4, 1)
    B = ensemble(base, 4, 2)
    ms = s0.MeanStart(A, B)
    ms.center = np.clip(ms.center + 2, 0, 24)          # 반폭 k = 2 만큼 밀어 둔다: 참 경로가 띠의 왼쪽 끝
    band = s0.pair_costs(A, B, method='band', k=2, mean_start=ms)
    full = s0.pair_costs(A, B, method='full')['cost']
    assert band['n_widened'] > 0 and np.abs(band['cost'] - full).max() < 1e-4 * full.max()


def test_band_edge_rule_is_not_a_guarantee():
    """중심선이 반폭보다 훨씬 어긋나면 띠 안의 짝 점수가 모두 나빠 짝을 잡지 않고, 가장자리에 닿지 않은 채 틀린 값이 남는다.
    그래서 가장자리 규칙은 점검일 뿐이고 정확성은 전체 DP·분기 한정과 견주어 잰다(연구 기획 E0). 값은 여전히 상계다."""
    base = chain(25, 5)
    A = ensemble(base, 4, 1)
    B = ensemble(base, 4, 2)
    ms = s0.MeanStart(A, B)
    ms.center = np.clip(ms.center + 6, 0, 24)
    band = s0.pair_costs(A, B, method='band', k=2, mean_start=ms)
    full = s0.pair_costs(A, B, method='full')['cost']
    assert (band['cost'] >= full - 1e-9).all() and (band['cost'] - full).max() > 0.1 * full.max()


def test_emd_assignment_and_lp_agree():
    rng = np.random.default_rng(2)
    C = rng.random((7, 7))
    P1 = s0.emd(C)
    P2 = s0.emd(C, np.ones(7), np.ones(7))            # 무게를 주면 일반 수송으로 푼다
    assert abs((P1 * C).sum() - (P2 * C).sum()) < 1e-9
    a, b = rng.random(5), rng.random(8)
    P = s0.emd(rng.random((5, 8)), a, b)
    assert np.allclose(P.sum(1), a / a.sum(), atol=1e-9) and np.allclose(P.sum(0), b / b.sum(), atol=1e-9)


def test_w2_basic_properties(two_ensembles):
    _, A, B = two_ensembles
    same = s0.ensemble_w2(A, A, method='mean')['W2']
    ab = s0.ensemble_w2(A, B, method='mean')['W2']
    ba = s0.ensemble_w2(B, A, method='mean')['W2']
    assert same < 1e-6 and abs(ab - ba) < 1e-6 * ab and ab > 0
    assert abs(s0.scaled(ab, 30, 31) - ab / math.sqrt(0.75 * 61)) < 1e-12


# ---------- 분해 ----------

def test_decomposition_identity_and_static_cost(two_ensembles):
    """맞춘 좌표에서 W₂² = 평균 몫 + 퍼짐 몫이 정확하고, 짝 없는 변의 비용은 0.75 s² 다."""
    _, A, B = two_ensembles
    r = s0.decompose(A, B)
    assert all(r['converged'])
    assert abs(r['identity_gap']) < 1e-9 * r['W2sq_aligned']
    assert abs(r['col_total'].sum() - r['W2sq_aligned']) < 1e-9 * r['W2sq_aligned']
    assert np.allclose(r['col_mean'] + r['col_spread'], r['col_total'], atol=1e-10)
    assert (~r['matched']).sum() == 1                  # 끼운 변 하나
    assert abs(r['static'] - r['static_expected']) < 1e-9
    assert abs(r['mean2'] - r['d2_means']) < 1e-6 * max(r['mean2'], 1e-12)   # 평균 점 사이 d² 와 같다
    # 표본 쌍마다 다시 최적화한 W₂ 는 맞춘 좌표의 W₂ 이하
    w = s0.ensemble_w2(A, B, method='mean')['W2sq']
    assert w <= r['W2sq_aligned'] + 1e-9
    # 변마다의 몫의 합
    assert abs(r['edgeA'].sum() + r['edgeB'].sum() - r['W2sq_aligned']) < 1e-9 * r['W2sq_aligned']


def test_decomposition_separates_shift_and_spread():
    """평균만 옮긴 쌍은 평균 몫이, 퍼짐만 넓힌 쌍은 퍼짐 몫이 커진다."""
    base = chain(20, 3)
    A = ensemble(base, 30, 1, sigma=0.05)
    shifted = base.clone()
    shifted[8:12] = torch.einsum('nxy,yz->nxz', shifted[8:12], fs.expmap(torch.tensor([0.0, 0.0, 0.4], dtype=torch.float64)))
    Bs = ensemble(shifted, 30, 2, sigma=0.05)
    Bw = ensemble(base, 30, 3, sigma=0.15)
    rs, rw = s0.decompose(A, Bs), s0.decompose(A, Bw)
    assert rs['mean2'] > rs['spread2'] and rw['spread2'] > rw['mean2']
    assert rw['var_ratio'] < 0.5                        # A 가 B 보다 좁다


def test_soft_correspondence(two_ensembles):
    _, A, B = two_ensembles
    pc = s0.pair_costs(A, B, method='mean', keep=True)
    plan = s0.w2_from_costs(pc['cost'])['plan']
    Pc, p, H = s0.soft_correspondence(plan, pc['match'], len(A), len(B))
    assert np.all(p <= 1 + 1e-9) and np.all(H >= -1e-12) and abs(Pc.sum() - p.sum()) < 1e-9


# ---------- 바닥값, 편향, 안정성 ----------

def test_fit_bias_recovers():
    Ns = np.array([25, 50, 100, 200, 400])
    a, b = s0.fit_bias(Ns, 0.3 + 2.0 / Ns)
    assert abs(a - 0.3) < 1e-10 and abs(b - 2.0) < 1e-8


def test_floor_smaller_than_between(two_ensembles):
    base, A, B = two_ensembles
    E = ensemble(base, 24, 9)
    groups = np.repeat([0, 1, 2], 8)
    fl = s0.replica_floors(E, groups, method='mean')
    between = s0.ensemble_w2(E[:8], B[:8], method='mean')['W2']
    assert len(fl) == 3 and max(f[2] for f in fl) < between


def test_bootstrap_and_triangle(two_ensembles):
    _, A, B = two_ensembles
    C = s0.pair_costs(A, B, method='mean')['cost']
    bt = s0.bootstrap_w2(C, B=30)
    assert 0 < bt['cv'] < 1
    assert s0.triangle_violation(1.0, 1.0, 1.5) == 0 and abs(s0.triangle_violation(1.0, 0.2, 1.5) - 0.3) < 1e-12


# ---------- 표본 ----------

def test_autocorr_ar1():
    """변마다 AR(1) 요동 (계수 φ) 이면 거리 자기상관이 φ^τ 이고 τ_int = 1/2 + φ/(1 − φ)."""
    phi, T, n = 0.8, 4000, 15
    rng = np.random.default_rng(0)
    z = np.zeros((T, n, 3))
    for t in range(1, T):
        z[t] = phi * z[t - 1] + math.sqrt(1 - phi ** 2) * rng.normal(size=(n, 3))
    base = chain(n, 1)
    E = torch.einsum('nxy,tnyz->tnxz', base, fs.expmap(torch.tensor(0.03 * z)))
    ac = s0.distance_autocorr(E, max_lag=60)
    expect = 0.5 + phi / (1 - phi)
    assert abs(ac['tau_int'] - expect) < 0.25 * expect
    assert ac['spacing'] == math.ceil(2 * ac['tau_int'])


def test_pooled_sample():
    g = np.repeat([0, 1, 2, 3], 100)
    idx, neff = s0.pooled_sample(np.zeros(400), g, 40, spacing=5)
    assert neff == 40 and len(set(g[idx])) == 4


# ---------- ESA ----------

def helix_atoms(L, seed=0, wobble=0.0):
    rng = np.random.default_rng(seed)
    t = np.arange(3 * L) * (2 * np.pi / 10.8)
    pts = np.stack([2.3 * np.cos(t), 2.3 * np.sin(t), 0.5 * np.arange(3 * L)], 1)
    pts += wobble * rng.normal(size=pts.shape)
    return pts


def test_esa_invariances():
    b1 = helix_atoms(20)
    R = fs.random_rotations(1, seed=4)[0].numpy()
    same = s0.esa_distance(b1, b1 @ R.T + 5.0, T=120)
    assert same['d'] < 1e-3
    other = s0.esa_distance(b1, helix_atoms(20, seed=1, wobble=1.5), T=120)
    assert other['d'] > 0.1
    # 점 밀도를 바꿔도(재매개화) 같은 곡선이면 가깝다
    dense = np.concatenate([np.linspace(b1[i], b1[i + 1], 4, endpoint=False) for i in range(len(b1) - 1)] + [b1[-1:]])
    assert s0.esa_distance(b1, dense, T=120)['d'] < 0.02
    rev = s0.esa_distance(helix_atoms(20, seed=1, wobble=1.5), b1, T=120)
    assert abs(rev['d'] - other['d']) < 0.05
