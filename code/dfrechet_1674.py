
import sys, os, time
ROOT = sys.argv[1]
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, pandas as pd, warnings
from numba import njit
warnings.filterwarnings('ignore')
import framespace as fs, stage0 as s0
from tmtools import tm_align as tma
g = pd.read_csv(ROOT + '/results/v_residue_unit_graded.tsv', sep='\t')[['A', 'B']]
z = np.load(ROOT + '/derived/xprot_structs.npz', allow_pickle=True)
def ca(name):
    if name + '|CA' in z.files:
        v = z[name + '|CA']; return v[0] if v.ndim == 3 else v
    return fs.read_pdb_backbone(ROOT + f'/derived/atlas_pdb/{name}.pdb')['CA'][0]
@njit(cache=True)
def dfd(P, Y):
    n, m = len(P), len(Y); F = np.empty((n, m))
    for i in range(n):
        for j in range(m):
            d = np.sqrt(((P[i] - Y[j]) ** 2).sum())
            if i == 0 and j == 0: F[i, j] = d
            elif i == 0: F[i, j] = max(d, F[i, j - 1])
            elif j == 0: F[i, j] = max(d, F[i - 1, j])
            else: F[i, j] = max(d, min(F[i - 1, j], F[i, j - 1], F[i - 1, j - 1]))
    # 결합 경로 되짚기
    i, j = n - 1, m - 1; pi = [i]; pj = [j]
    while i > 0 or j > 0:
        if i == 0: j -= 1
        elif j == 0: i -= 1
        else:
            a, b, c = F[i - 1, j - 1], F[i - 1, j], F[i, j - 1]
            if a <= b and a <= c: i -= 1; j -= 1
            elif b <= c: i -= 1
            else: j -= 1
        pi.append(i); pj.append(j)
    return F[n - 1, m - 1], np.array(pi[::-1]), np.array(pj[::-1])
def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); H = (Q - qc).T @ (P - pc); U, S, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T; return R, pc - R @ qc
OUT = ROOT + '/results/dfrechet_1674.tsv'
done = set()
if os.path.exists(OUT): d0 = pd.read_csv(OUT, sep='\t'); done = set(zip(d0.A, d0.B))
else: open(OUT, 'w').write('A\tB\tdf\tdf_tmpose\tsec\n')
buf, t0, k = [], time.time(), 0
for A, B in g.itertuples(index=False):
    if (A, B) in done: continue
    t1 = time.time(); P, Q = ca(A), ca(B)
    res = tma(P, Q, 'A' * len(P), 'A' * len(Q)); U0, t0_ = res.u.T, -res.u.T @ res.t
    poses = list(s0.start_poses(P, Q, n_random=0)) + [(U0, t0_)]
    vals = [dfd(P, Q @ U.T + t)[0] for U, t in poses]
    order = np.argsort(vals)[:8]; best = min(vals); df_tm = vals[-1]
    for o in order:
        U, t = poses[o]
        for _ in range(15):
            v, pi, pj = dfd(P, Q @ U.T + t)
            best = min(best, v)
            U2, t2 = kabsch(P[pi], Q[pj])
            v2 = dfd(P, Q @ U2.T + t2)[0]
            if v2 >= v - 1e-9: break
            U, t = U2, t2; best = min(best, v2)
    buf.append(f'{A}\t{B}\t{best}\t{df_tm}\t{time.time() - t1:.2f}'); k += 1
    if len(buf) >= 50:
        open(OUT, 'a').write('\n'.join(buf) + '\n'); buf = []; print(time.strftime('%H:%M:%S'), k, round((time.time() - t0) / k, 2), flush=True)
if buf: open(OUT, 'a').write('\n'.join(buf) + '\n')
print('FINISHED', k, flush=True)
