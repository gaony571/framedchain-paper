
import sys, os, json, time
ROOT = sys.argv[1]; RAW = ROOT + '/raw'
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, pandas as pd, h5py, warnings
warnings.filterwarnings('ignore')
import stage0 as s0
from tmtools import tm_align as tma
SPD = ROOT + '/derived/mdcath_spread'
pairs = pd.read_csv(ROOT + '/results/pilot_mdcath_homolog_pairs.tsv', sep='\t')[['a', 'b']]
doms = sorted(set(pairs.a) | set(pairs.b))
def load_ca():
    out = {}
    for f in os.listdir(SPD):
        if f.endswith('.npz'):
            try:
                z = np.load(SPD + '/' + f)
                if 'ca0' in z.files: out[f[:-4]] = z['ca0']
            except Exception:
                pass
    return out
def events(ii, jj):
    ev = []
    for k in range(1, len(ii)):
        da, db = ii[k] - ii[k - 1] - 1, jj[k] - jj[k - 1] - 1
        if da > 0: ev.append(['a', int(ii[k - 1] + 1), int(jj[k - 1] + 1), int(da)])
        if db > 0: ev.append(['b', int(ii[k - 1] + 1), int(jj[k - 1] + 1), int(db)])
    return ev
OUT = ROOT + '/results/p2_gap_consensus_raw.tsv'
done = set()
if os.path.exists(OUT):
    d0 = pd.read_csv(OUT, sep='\t'); done = set(zip(d0.a, d0.b))
else:
    open(OUT, 'w').write('a\tb\tn\tm\ttm_a\ttm_b\tS_ours\tk_tma\tk_ours\tev_tma\tev_ours\n')
buf, t0, n = [], time.time(), 0
while True:
  CA = load_ca(); todo = [(a, b) for a, b in pairs.itertuples(index=False) if (a, b) not in done and a in CA and b in CA]
  print(time.strftime('%H:%M:%S'), 'available pairs', len(todo), 'CA', len(CA), flush=True)
  for a, b in todo:
      P, Q = CA[a], CA[b]
      res = tma(P, Q, 'A' * len(P), 'A' * len(Q))
      ti, tj, i, j = [], [], 0, 0
      for x, y in zip(res.seqxA, res.seqyA):
          if x != '-' and y != '-': ti.append(i); tj.append(j)
          i += x != '-'; j += y != '-'
      r = s0.tm_align(P, Q, 4.0, extra_poses=[(res.u.T, -res.u.T @ res.t)], keep=16, search=True)
      oi, oj = np.nonzero(r['match'])
      buf.append('\t'.join(map(str, [a, b, len(P), len(Q), res.tm_norm_chain1, res.tm_norm_chain2, r['S'], len(ti), len(oi), json.dumps(events(ti, tj)), json.dumps(events(list(oi), list(oj)))])))
      n += 1
      if len(buf) >= 100:
          open(OUT, 'a').write('\n'.join(buf) + '\n'); buf = []
          print(time.strftime('%H:%M:%S'), n, round((time.time() - t0) / n, 2), flush=True)
  if buf: open(OUT, 'a').write('\n'.join(buf) + '\n'); buf = []
  for a, b in todo: done.add((a, b))
  if os.path.exists(ROOT + '/results/logs/mdcath_spread.done') and not todo: break
  time.sleep(600)
print('FINISHED', n, flush=True)
