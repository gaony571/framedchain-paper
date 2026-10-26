
import sys, os, time, json
ROOT = sys.argv[1]
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import framespace as fs, stage0 as s0
from tmtools import tm_align as tma
g = pd.read_csv(ROOT + '/results/v_residue_unit_graded.tsv', sep='\t')[['A', 'B', 'src', 'TM']]
z = np.load(ROOT + '/derived/xprot_structs.npz', allow_pickle=True)
def bb(name):
    k = name + '|CA'
    if k in z.files:
        f = lambda a: (z[name + '|' + a][0] if z[name + '|' + a].ndim == 3 else z[name + '|' + a])
        return f('N'), f('CA'), f('C'), f('O')
    b = fs.read_pdb_backbone(ROOT + f'/derived/atlas_pdb/{name}.pdb'); return b['N'][0], b['CA'][0], b['C'][0], b['O'][0]
def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); H = (Q - qc).T @ (P - pc); U, S, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T; return R, pc - R @ qc
def gdt_ts(P, Q, ii, jj, R0, t0, L):
    out = []
    for c in (1.0, 2.0, 4.0, 8.0):
        best = 0; seeds = [(R0, t0)]
        k = len(ii)
        for Lf in (k, max(4, k // 2), max(4, k // 4), 4):
            for s in np.unique(np.linspace(0, max(k - Lf, 0), 6).round().astype(int)):
                if k >= 3: seeds.append(kabsch(P[ii[s:s + Lf]], Q[jj[s:s + Lf]]))
        for R, t in seeds:
            for _ in range(10):
                d = np.linalg.norm(P[ii] - (Q[jj] @ R.T + t), axis=1); m = d < c
                best = max(best, int(m.sum()))
                if m.sum() < 3: break
                Rn, tn = kabsch(P[ii[m]], Q[jj[m]])
                if np.allclose(Rn, R) and np.allclose(tn, t): break
                R, t = Rn, tn
        out.append(best / L)
    return float(np.mean(out))
def lddt_ca(P, Q, ii, jj):
    """기준 = Q (B). 기준 안 15 Å 안의 쌍마다, 대응한 P 쌍의 거리 차가 0.5/1/2/4 Å 안이면 보존. 짝 없는 잔기는 보존 안 됨."""
    m = len(Q); DQ = np.linalg.norm(Q[:, None] - Q[None], axis=-1); inc = (DQ < 15) & ~np.eye(m, dtype=bool)
    mp = -np.ones(m, int); mp[jj] = ii
    a, b = np.nonzero(inc); ok = (mp[a] >= 0) & (mp[b] >= 0)
    dp = np.full(len(a), np.inf); dp[ok] = np.linalg.norm(P[mp[a[ok]]] - P[mp[b[ok]]], axis=1)
    diff = np.abs(dp - DQ[a, b]); sc = np.mean([(diff < t) for t in (0.5, 1, 2, 4)], axis=0)
    per = np.bincount(a, weights=sc, minlength=m) / np.maximum(np.bincount(a, minlength=m), 1)
    return float(per.mean())
OUT = ROOT + '/results/hu1_rulers_1674.tsv'
COLS = ['A', 'B', 'src', 'n', 'm', 'tm_a', 'tm_b', 'rmsd_tma', 'n_ali', 'gdt_ts_B', 'lddt_B', 'S_tm4', 'D2_l4', 'D2_l8', 'D2_l16', 'D2_inf', 'esa', 'sec']
done = set()
if os.path.exists(OUT):
    d0 = pd.read_csv(OUT, sep='\t'); done = set(zip(d0.A, d0.B))
else:
    open(OUT, 'w').write('\t'.join(COLS) + '\n')
buf, t_start, k = [], time.time(), 0
for A, B, src, _ in g.itertuples(index=False):
    if (A, B) in done: continue
    t0 = time.time()
    try:
        NA, PA, CA_, OA = bb(A); NB, PB, CB, OB = bb(B)
        res = tma(PA, PB, 'A' * len(PA), 'A' * len(PB))
        ii, jj, i, j = [], [], 0, 0
        for x, y in zip(res.seqxA, res.seqyA):
            if x != '-' and y != '-': ii.append(i); jj.append(j)
            i += x != '-'; j += y != '-'
        ii, jj = np.array(ii), np.array(jj); U, t = res.u.T, -res.u.T @ res.t
        gdt = gdt_ts(PA, PB, ii, jj, U, t, len(PB)); ld = lddt_ca(PA, PB, ii, jj)
        stm = s0.tm_align(PA, PB, 4.0, extra_poses=[(U, t)], keep=16, search=True)['S']
        xa = fs.coordinates_from_backbone(PA, CA_, OA); xb = fs.coordinates_from_backbone(PB, CB, OB)
        dl = [s0.wl_align(xa, xb, ell=l, keep=8)['D2'] for l in (4.0, 8.0, 16.0)]
        dinf = s0.wl_align(xa, xb, ell=float('inf'))['D2']
        try:
            e = float(s0.esa_distance(s0.backbone_curve(NA, PA, CA_), s0.backbone_curve(NB, PB, CB), T=100, iters=2)['d'])
        except Exception:
            e = float('nan')
        row = [A, B, src, len(PA), len(PB), res.tm_norm_chain1, res.tm_norm_chain2, res.rmsd, len(ii), gdt, ld, stm] + dl + [dinf, e, time.time() - t0]
        buf.append('\t'.join(str(v) for v in row)); k += 1
    except Exception as ex:
        print('ERR', A, B, repr(ex)[:200], flush=True)
    if len(buf) >= 20:
        open(OUT, 'a').write('\n'.join(buf) + '\n'); buf = []
        print(time.strftime('%H:%M:%S'), k, round((time.time() - t_start) / max(k, 1), 2), flush=True)
if buf: open(OUT, 'a').write('\n'.join(buf) + '\n')
print('FINISHED', k, flush=True)
