"""Self-contained demonstration on synthetic backbones (no downloads, CPU, about a minute).

It uses only the as-run libraries of this repository (framespace 0.2.0 in framespace_pkg/ and the research module
framespace_wip/stage0) and applies the paper's computations to small synthetic ensembles:

  1. invariances of the distance D (global rotation, inserted zero terms) and the triangle inequality on a triple;
  2. the reduced computation used in all ensemble statistics (mean-point starts + banded dynamic programming;
     the function `reduced_d` below is `red_pairs` of code/r38_s0.py with the data-loading parts removed);
  3. the energy-distance statistic between ensembles, Δ² = 2C − B_P − B_Q with one mean per pair of runs
     (`estats` of code/r38_s0.py), for a null comparison (same distribution), a shift (larger fluctuations),
     and two ensembles of different length (an insertion of 4 edges).

Usage (from the repository root):
    python examples/demo_synthetic.py
The printed values should match examples/demo_expected_output.txt to the printed precision (the run is
deterministic: fixed seeds, float64, CPU).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)
sys.path[:0] = [os.path.join(REPO, 'framespace_pkg'), os.path.join(REPO, 'framespace_wip')]
import numpy as np, torch
import framespace as fs
from framespace.distance import sq_norm
from framespace.geometry import expmap
from stage0.band import MeanStart, align_band, band_windows

torch.set_num_threads(1)
DT = torch.float64
N_EDGES, RUNS, FRAMES, NPAIR, K0 = 40, 3, 30, 200, 5


def chain(n, seed):
    """A smooth synthetic backbone in the space: frames from a random walk of small rotations, trans edges (3.83 Å)."""
    g = torch.Generator().manual_seed(seed)
    R = [torch.eye(3, dtype=DT)]
    for _ in range(n - 1):
        R.append(expmap(0.6 * torch.randn(3, generator=g, dtype=DT)) @ R[-1])
    return torch.as_tensor(fs.to_coordinates(torch.stack(R).numpy(), np.full(n, fs.L0)), dtype=DT)


def ensemble(x, sigma, runs, frames, seed):
    """runs × frames structures: each edge frame rotated by an independent small rotation (std sigma rad per axis),
    and a random global rotation per structure (the distance must not see it)."""
    g = torch.Generator().manual_seed(seed)
    n = len(x)
    out = []
    for _ in range(runs * frames):
        P = expmap(sigma * torch.randn(n, 3, generator=g, dtype=DT))
        G = expmap(torch.randn(3, generator=g, dtype=DT))
        out.append(torch.einsum('xy,nyz,nzw->nxw', G, P, x))
    return torch.stack(out)


def full_d(x, y):
    a = fs.align(x, y, n_random=127, seed=0)
    return float((sq_norm(x) + sq_norm(y) - 2 * a.score).clamp_min(0).sqrt())


def reduced_d(ms, ia, ib, k0=K0, iters=10):
    """Reduced computation for the pairs (ia[p], ib[p]); identical to red_pairs in code/r38_s0.py
    (memory batching removed). Returns d (P,) and the number of pairs whose band was widened."""
    EA, EB = ms.EA, ms.EB; m = EB.shape[1]
    ia = torch.as_tensor(np.asarray(ia)); ib = torch.as_tensor(np.asarray(ib))
    S = EA.new_zeros(len(ia)); todo = torch.arange(len(ia)); kk = k0; widened = 0
    while len(todo):
        lo, w = band_windows(ms.center, m, kk)
        sc, _, M, edge = align_band(EA[ia[todo]], EB[ib[todo]], lo, w, ms.starts(ia[todo], ib[todo], n_local=1), iters=iters)
        S[todo] = sc
        if w >= m: break
        todo = todo[edge]
        if kk == k0: widened = int(len(todo))
        kk *= 2
    d2 = (sq_norm(EA)[ia] + sq_norm(EB)[ib] - 2 * S).clamp_min(0)
    return d2.sqrt().numpy(), widened


def upairs(a, b, rng, k=NPAIR, nn=FRAMES):
    f = rng.choice(nn * nn, k, replace=False); return a * nn + f // nn, b * nn + f % nn


def estats(C, BP, BQ, ne):
    D2 = 2 * C - BP - BQ; w = (BP + BQ) / 2
    return dict(C=C, B_P=BP, B_Q=BQ, D2=D2, wbar=w, ratio=float(np.sqrt(max((C / w) ** 2 - 1, 0))), d2n=D2 / np.sqrt(ne))


def compare(EA, EB, seed):
    """Δ² with one mean per pair of runs: cross = runs of A × runs of B, floors = pairs of different runs on one side."""
    rng = np.random.default_rng(seed)
    def mean_d(ms, pairs):
        out = []
        for a, b in pairs:
            ia, ib = upairs(a, b, rng); out.append(float(reduced_d(ms, ia, ib)[0].mean()))
        return float(np.mean(out))
    within = [(a, b) for a in range(RUNS) for b in range(RUNS) if a < b]
    C = mean_d(MeanStart(EA, EB), [(a, b) for a in range(RUNS) for b in range(RUNS)])
    BP = mean_d(MeanStart(EA, EA), within); BQ = mean_d(MeanStart(EB, EB), within)
    return estats(C, BP, BQ, ne=(EA.shape[1] + EB.shape[1]) / 2)


def main():
    x, y, z = chain(N_EDGES, 1), chain(N_EDGES, 2), chain(N_EDGES - 6, 3)
    print('1. distance D (full search, 128 starts)')
    G = expmap(torch.tensor([0.3, -1.1, 2.0], dtype=DT))
    print(f'   D(x, Gx)                    = {full_d(x, torch.einsum("xy,nyz->nxz", G, x)):.6f}   (global rotation)')
    xz = torch.cat([x[:20], torch.zeros(3, 3, 3, dtype=DT), x[20:]])
    print(f'   D(x, x with 3 zero terms)   = {full_d(x, xz):.6f}   (zero terms are the same point)')
    dxy, dyz, dxz = full_d(x, y), full_d(y, z), full_d(x, z)
    print(f'   D(x,y) D(y,z) D(x,z)        = {dxy:.4f} {dyz:.4f} {dxz:.4f};  D(x,z) <= D(x,y) + D(y,z): {dxz <= dxy + dyz + 1e-9}')
    print(f'   normalised D(x,y) in [0,1]  = {dxy / float((sq_norm(x) + sq_norm(y)).sqrt()):.4f}')

    P1 = ensemble(x, 0.10, RUNS, FRAMES, 11); P2 = ensemble(x, 0.10, RUNS, FRAMES, 12)
    Q = ensemble(x, 0.20, RUNS, FRAMES, 13)
    xi = torch.cat([x[:20], torch.einsum('xy,nyz->nxz', 2 * x[19], chain(4, 4)), x[20:]])   # insertion of 4 edges
    I = ensemble(xi, 0.10, RUNS, FRAMES, 14)

    print('2. reduced computation vs full search on 10 frame pairs (upper bound; relative difference)')
    ms = MeanStart(P1, I); ia, ib = np.arange(10), np.arange(10) + FRAMES
    dr, wid = reduced_d(ms, ia, ib)
    df = np.array([full_d(P1[i], I[j]) for i, j in zip(ia, ib)])
    print(f'   median (d_reduced - d_full)/d_full = {np.median((dr - df) / df):+.2e}; max = {np.max((dr - df) / df):+.2e}; widened bands: {wid}')

    print('3. energy-distance statistic (one mean per run pair, 200 frame pairs each)')
    for name, (A, B) in [('null  P vs P (same distribution)', (P1, P2)), ('shift P vs Q (sigma 0.10 -> 0.20)', (P1, Q)),
                         ('insertion P vs P+4 edges', (P1, I))]:
        s = compare(A, B, 20261026)
        print(f'   {name:34s} C={s["C"]:.4f} B_P={s["B_P"]:.4f} B_Q={s["B_Q"]:.4f} '
              f'D2={s["D2"]:+.4f} ratio={s["ratio"]:.3f} D2/sqrt(n)={s["d2n"]:+.5f}')


if __name__ == '__main__':
    main()
