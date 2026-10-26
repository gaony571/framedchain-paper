
import sys, os, time, json, itertools
ROOT = sys.argv[1]
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import stage0 as s0
from tmtools import tm_align as tma
z = np.load(ROOT + '/derived/xprot_structs.npz', allow_pickle=True)
setup = json.load(open(ROOT + '/results/pu_tm_factorial_setup.json'))
blk = setup['block']
CA = {n: (z[n + '|CA'][0] if z[n + '|CA'].ndim == 3 else z[n + '|CA']) for n in blk}
PU = ROOT + '/results/pu_tm_factorial_pairs.tsv'
COLS = ['a', 'b', 'n', 'm', 'tm_a', 'tm_b', 'S_fix_low', 'S_fix_high', 'S_var_low', 'S_var_high', 'sec']
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
    dv = s0.d0_tm(0.5 * (n + m)); out = [a, b, n, m, res.tm_norm_chain1, res.tm_norm_chain2]
    for d0 in (setup['d0_fixed'], dv):
        out.append(s0.tm_align(P, Q, d0, keep=4, search=False)['S'])
        out.append(s0.tm_align(P, Q, d0, extra_poses=ext, keep=16, search=True)['S'])
    out.append(time.time() - t0); buf.append('\t'.join(str(v) for v in out)); k += 1
    if len(buf) >= 50:
        with open(PU, 'a') as f: f.write('\n'.join(buf) + '\n')
        buf = []
        print(time.strftime('%H:%M:%S'), k, round((time.time() - t_start) / k, 2), flush=True)
if buf:
    with open(PU, 'a') as f: f.write('\n'.join(buf) + '\n')
print('FINISHED', time.strftime('%H:%M:%S'), k, flush=True)
