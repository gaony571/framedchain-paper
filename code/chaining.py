
"""다중 도메인 사슬 잇기 (prereg round17). results/chaining.json, derived/chaining.tsv."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, time, gzip, urllib.request
ROOT = sys.argv[1]; NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 150
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, pandas as pd
from tmtools import tm_align as tma
import stage0 as s0
RAW = ROOT + '/raw/chaining'; os.makedirs(RAW + '/pdb', exist_ok=True)
cath = {}
for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
    if l.startswith('#') or not l.strip(): continue
    f = l.split(); cath[f[0]] = f[1:5]
cand = []
for l in open(RAW + '/cath-domain-boundaries.txt'):
    if l.startswith('#'): continue
    t = l.split()
    if len(t) < 4 or t[1] != 'D02' or t[2] != 'F00': continue
    i = 3; doms = []
    try:
        for d in range(2):
            ns = int(t[i]); i += 1; segs = []
            for _ in range(ns): segs.append((t[i], int(t[i + 1]), t[i + 2], int(t[i + 4]), t[i + 5])); i += 6
            doms.append(segs)
    except Exception: continue
    if any(len(s) != 1 for s in doms): continue
    (a, b) = doms; ch = t[0]
    if a[0][1] > b[0][1]: a, b = b, a
    la, lb = a[0][3] - a[0][1] + 1, b[0][3] - b[0][1] + 1
    if not (50 <= la <= 200 and 50 <= lb <= 200) or a[0][3] >= b[0][1]: continue
    da, db = ch + '01', ch + '02'
    if da not in cath or db not in cath: continue
    Ta, Tb = '.'.join(cath[da][:3]), '.'.join(cath[db][:3])
    if Ta == Tb: continue
    cand.append(dict(chain=ch, a=a[0], b=b[0], Ta=Ta, Tb=Tb, key=tuple(sorted((Ta, Tb)))))
rng = np.random.default_rng(0); rng.shuffle(cand); seen, sel = set(), []
for c in cand:
    if c['key'] in seen: continue
    seen.add(c['key']); sel.append(c)
print('candidates', len(cand), 'topology pairs', len(sel), flush=True)
def ca_of(pdb, chain):
    p = f'{RAW}/pdb/{pdb}.pdb.gz'
    if not os.path.exists(p): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb}.pdb.gz', p)
    out = {}
    for l in gzip.open(p, 'rt'):
        if l.startswith('ENDMDL'): break
        if l.startswith('ATOM') and l[12:16].strip() == 'CA' and l[21] == chain and l[16] in ' A':
            key = (int(l[22:26]), l[26].strip()); out.setdefault(key, (np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l[17:20]))
    return out
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V', MSE='M')
def piece(ca, lo, hi):
    ks = [k for k in ca if lo <= k[0] <= hi]; ks.sort()
    return np.array([ca[k][0] for k in ks]), ''.join(AA.get(ca[k][1], 'X') for k in ks)
rows = []; t0 = time.time()
for c in sel:
    if len(rows) >= NMAX: break
    try:
        ca = ca_of(c['chain'][:4].lower(), c['chain'][4])
        X, sx = piece(ca, c['a'][1], c['a'][3]); Y, sy = piece(ca, c['b'][1], c['b'][3]); XY, sxy = piece(ca, c['a'][1], c['b'][3])
        if min(len(X), len(Y)) < 40 or len(XY) < len(X) + len(Y): continue
        r = dict(chain=c['chain'], Ta=c['Ta'], Tb=c['Tb'], nX=len(X), nY=len(Y), nXY=len(XY))
        P = dict(X=(X, sx), Y=(Y, sy), XY=(XY, sxy)); d0f = s0.d0_tm((len(X) + len(Y)) / 2)
        for u, v in (('X', 'Y'), ('X', 'XY'), ('Y', 'XY')):
            (A, sa), (B, sb) = P[u], P[v]; t = tma(A, B, sa, sb); tm1, tm2 = t.tm_norm_chain1, t.tm_norm_chain2
            short = tm1 if len(A) <= len(B) else tm2; long_ = tm2 if len(A) <= len(B) else tm1
            Lb = (len(A) + len(B)) / 2; dn = s0.tm_align(A, B, d0=s0.d0_tm(Lb)); df = s0.tm_align(A, B, d0=d0f)
            r.update({f'{u}{v}_tms': 1 - short, f'{u}{v}_tml': 1 - long_, f'{u}{v}_dn': float(np.sqrt(max(dn['D2'], 0) / Lb)), f'{u}{v}_df': float(np.sqrt(max(df['D2'], 0)))})
        rows.append(r); print(time.strftime('%H:%M:%S'), len(rows), c['chain'], len(XY), round(time.time() - t0, 1), flush=True)
    except Exception as e: print('ERR', c['chain'], repr(e)[:150], flush=True)
D = pd.DataFrame(rows); D.to_csv(ROOT + '/derived/chaining.tsv', sep='\t', index=False)
out = dict(n=len(D), n_candidates=len(cand), n_topology_pairs=len(sel))
TAU = dict(tms=0.5, tml=0.5, dn=np.sqrt(0.5))
def bci(v):
    v = np.asarray(v, float); b = [v[rng.integers(0, len(v), len(v))].mean() for _ in range(2000)]; return [round(float(v.mean()), 3)] + [round(float(x), 3) for x in np.percentile(b, [2.5, 97.5])]
for m in ('tms', 'tml', 'dn', 'df'):
    xy, xx, yy = D[f'XY_{m}'], D[f'XXY_{m}'], D[f'YXY_{m}']
    o = dict(median_X_Y=round(float(xy.median()), 3), median_X_XY=round(float(xx.median()), 3), median_Y_XY=round(float(yy.median()), 3),
             triangle_violation=bci(xy > xx + yy + 1e-6), triangle_violation_rel1pct=bci(xy > (xx + yy) * 1.01))
    if m in TAU:
        tau = TAU[m]; o['single_link_merge'] = bci((np.maximum(xx, yy) < tau) | (xy < tau)); o['X_Y_same_fold_direct'] = bci(xy < tau)
        # 평균 연결 (세 점): 가장 가까운 짝을 먼저 합치고 남은 점은 평균 거리
        mg = []
        for a, b_, c_ in zip(xy, xx, yy):
            d = {('X', 'Y'): a, ('X', 'XY'): b_, ('Y', 'XY'): c_}; (p1, p2), dmin = min(d.items(), key=lambda kv: kv[1])
            if dmin >= tau: mg.append(False); continue
            rest = ({'X', 'Y', 'XY'} - {p1, p2}).pop(); dm = np.mean([d[tuple(sorted((rest, q), key=['X', 'Y', 'XY'].index))] for q in (p1, p2)])
            mg.append(bool({'X', 'Y'} <= {p1, p2} or dm < tau))
        o['average_link_merge'] = bci(mg)
    out[m] = o
out['note'] = 'tms = 1 − TM_short, tml = 1 − TM_long (TM-align), dn = D_TM/√L̄ (d₀(L̄)), df = D_TM 고정 d₀ (세 짝 공통 d₀(평균 도메인 길이), 정리가 서는 판). τ: TM 0.5 에 맞춤.'
out['written'] = time.strftime('%F %T'); json.dump(out, open(ROOT + '/results/chaining.json', 'w'), indent=1); print(json.dumps(out, indent=1))
open(ROOT + '/results/logs/chaining.done', 'w').write(time.strftime('%F %T'))
