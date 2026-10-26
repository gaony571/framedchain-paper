
import sys, os, time, json, itertools
ROOT = sys.argv[1]
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import stage0 as s0
from tmtools import tm_align as tma
setup = json.load(open(ROOT + '/results/pu_fold208_setup.json'))
blk = setup['block']; zc = np.load(ROOT + '/' + setup['ca_npz']); CA = {n: zc[n] for n in blk}
PU = ROOT + '/results/pu_fold208_pairs.tsv'
COLS = ['a', 'b', 'n', 'm', 'tm_a', 'tm_b', 'S_fix_hi', 'S_avg_hi', 'sec']
done = set()
if os.path.exists(PU) and os.path.getsize(PU) > 0:
    d = pd.read_csv(PU, sep='\t'); done = set(zip(d.a, d.b))
else:
    open(PU, 'w').write('\t'.join(COLS) + '\n')
buf, t_start, k = [], time.time(), 0
for a, b in itertools.combinations(blk, 2):
    if (a, b) in done: continue
    t0 = time.time(); P, Q = CA[a], CA[b]; n, m = len(P), len(Q)
    res = tma(P, Q, 'A' * n, 'A' * m); ext = [(res.u.T, -res.u.T @ res.t)]
    out = [a, b, n, m, res.tm_norm_chain1, res.tm_norm_chain2,
           s0.tm_align(P, Q, setup['d0_fixed'], extra_poses=ext, keep=16, search=True)['S'],
           s0.tm_align(P, Q, s0.d0_tm(0.5 * (n + m)), extra_poses=ext, keep=16, search=True)['S'], 0]
    out[-1] = time.time() - t0; buf.append('\t'.join(str(v) for v in out)); k += 1
    if len(buf) >= 50:
        with open(PU, 'a') as f: f.write('\n'.join(buf) + '\n')
        buf = []
        print(time.strftime('%H:%M:%S'), k, round((time.time() - t_start) / k, 2), flush=True)
if buf:
    with open(PU, 'a') as f: f.write('\n'.join(buf) + '\n')
print('FINISHED', time.strftime('%H:%M:%S'), k, flush=True)
