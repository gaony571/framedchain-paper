
"""EV·E0·사다리 공통 mdCATH 좌표 뽑기 (results/mdcath_5T_setup.json). h5 는 도메인마다 한 번만 연다.
저장(도메인마다 derived/mdcath_5T/<dom>.npz, 압축): 온도·반복마다 Cα·C·O 좌표(frame_step 간격, float32), 모든 프레임의 Q(float16), 같은 간격의 Rg 와 DSSP 코드, 서열, 원래 프레임 수."""
import sys, os, json, time
ROOT = sys.argv[1]; RAW = ROOT + '/raw'; OUTD = ROOT + '/derived/mdcath_5T'; os.makedirs(OUTD, exist_ok=True)
import numpy as np, h5py
st = json.load(open(ROOT + '/results/mdcath_5T_setup.json')); STEP = st['frame_step']
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', HSD='H', HSE='H', HSP='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
DS = {c: k for k, c in enumerate(' HBEGITSP-')}
def parse(h):
    txt = h['pdbProteinAtoms'][()]; txt = txt.decode() if isinstance(txt, bytes) else str(txt)
    res, order, nm, k = {}, [], {}, 0
    for line in txt.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            key = (line[21], line[22:27]); an = line[12:16].strip()
            if key not in res: res[key] = {}; order.append(key); nm[key] = line[17:21].strip()
            if an in ('N', 'CA', 'C', 'O') and an not in res[key]: res[key][an] = k
            k += 1
    keep = [r for r in order if all(a in res[r] for a in ('N', 'CA', 'C', 'O'))]
    return np.array([[res[r][a] for a in ('CA', 'C', 'O')] for r in keep]), ''.join(AA.get(nm[r], 'X') for r in keep)
t0 = time.time(); n = 0
for d in st['domains']:
    out = f'{OUTD}/{d}.npz'
    if os.path.exists(out): continue
    try:
        with h5py.File(RAW + f'/mdcath/mdcath_dataset_{d}.h5', 'r') as f:
            h = f[d]; idx, seq = parse(h); ca = idx[:, 0]
            X0 = h['320']['0']['coords'][0][ca]; D0 = np.linalg.norm(X0[:, None] - X0[None], axis=-1); L = len(ca)
            ii, jj = np.triu_indices(L, 4); m = D0[ii, jj] < 8.0; ii, jj, dn = ii[m], jj[m], D0[ii[m], jj[m]]
            res = dict(seq=np.array(seq), native_i=ii, native_j=jj, native_d=dn.astype(np.float32))
            for T in st['temps']:
                for r in sorted(h[T].keys(), key=int):
                    C = h[T][r]['coords'][:]; k = f'{T}_{r}'
                    Xc = C[:, ca]; dd = np.linalg.norm(Xc[:, ii] - Xc[:, jj], axis=-1)
                    res[k + '_Q'] = (dd < 1.2 * dn).mean(1).astype(np.float16); res[k + '_nframes'] = np.int32(len(C))
                    res[k + '_X'] = C[::STEP][:, idx.ravel()].reshape(-1, L, 3, 3).astype(np.float32)
                    g = h[T][r]
                    if 'gyrationRadius' in g: res[k + '_rg'] = np.asarray(g['gyrationRadius'][::STEP], np.float32)
                    if 'dssp' in g:
                        ds = np.asarray(g['dssp'][::STEP]); ds = np.vectorize(lambda x: x.decode() if isinstance(x, bytes) else str(x))(ds)
                        res[k + '_dssp'] = np.vectorize(lambda c: DS.get(c, 0))(ds).astype(np.uint8)
        np.savez_compressed(out, **res); n += 1
        print(time.strftime('%H:%M:%S'), n, d, round((time.time() - t0) / n, 1), 's/dom', flush=True)
    except Exception as e:
        print('ERR', d, repr(e)[:200], flush=True)
print('FINISHED', n, flush=True)
open(ROOT + '/results/logs/mdcath_5T.done', 'w').write(time.strftime('%Y-%m-%d %H:%M:%S'))
