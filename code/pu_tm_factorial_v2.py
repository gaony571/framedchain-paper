
import sys, os, time, json, itertools
ROOT = sys.argv[1]
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import stage0 as s0
from tmtools import tm_align as tma
setup = json.load(open(ROOT + '/' + sys.argv[2]))
blk = setup['block']
zc = np.load(ROOT + '/' + setup['ca_npz'])
CA = {n: zc[n.replace('|', '__')] for n in blk}
PU = ROOT + '/' + sys.argv[3]
D0 = setup['d0_fixed']
COLS = ['a', 'b', 'n', 'm', 'tm_a', 'tm_b', 'S_fix_k1', 'S_fix_k4', 'S_fix_k16', 'S_fix_hi', 'S_fix_max',
        'S_short_k4', 'S_short_hi', 'S_long_k4', 'S_long_hi', 'S_avg_k4', 'S_avg_hi', 'sec']
done = set()
if os.path.exists(PU) and os.path.getsize(PU) > 0:
    d = pd.read_csv(PU, sep='\t'); done = set(zip(d.a, d.b))
else:
    open(PU, 'w').write('\t'.join(COLS) + '\n')
buf, t_start, k = [], time.time(), 0
for a, b in itertools.combinations(blk, 2):
    if (a, b) in done: continue
    t0 = time.time(); P, Q = CA[a], CA[b]; n, m = len(P), len(Q)
    res = tma(P, Q, 'A' * n, 'A' * m)
    ext = [(res.u.T, -res.u.T @ res.t)]
    big = s0.start_poses(P, Q, per_side=20, n_thread=48, n_random=32)
    out = [a, b, n, m, res.tm_norm_chain1, res.tm_norm_chain2]
    out.append(s0.tm_align(P, Q, D0, keep=1, search=False)['S'])
    out.append(s0.tm_align(P, Q, D0, keep=4, search=False)['S'])
    out.append(s0.tm_align(P, Q, D0, keep=16, search=False)['S'])
    out.append(s0.tm_align(P, Q, D0, extra_poses=ext, keep=16, search=True)['S'])
    out.append(s0.tm_align(P, Q, D0, poses=big, extra_poses=ext, keep=48, search=True)['S'])
    for L in (min(n, m), max(n, m), 0.5 * (n + m)):
        d0 = s0.d0_tm(L)
        out.append(s0.tm_align(P, Q, d0, keep=4, search=False)['S'])
        out.append(s0.tm_align(P, Q, d0, extra_poses=ext, keep=16, search=True)['S'])
    out.append(time.time() - t0); buf.append('\t'.join(str(v) for v in out)); k += 1
    if len(buf) >= 25:
        with open(PU, 'a') as f: f.write('\n'.join(buf) + '\n')
        buf = []
        print(time.strftime('%H:%M:%S'), k, round((time.time() - t_start) / k, 2), flush=True)
if buf:
    with open(PU, 'a') as f: f.write('\n'.join(buf) + '\n')
print('FINISHED', time.strftime('%H:%M:%S'), k, flush=True)
