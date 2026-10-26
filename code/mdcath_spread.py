
import sys, os, time, json
ROOT = sys.argv[1]; RAW = ROOT + '/raw'
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, h5py, torch
import framespace as fs
OUTD = ROOT + '/derived/mdcath_spread' + (sys.argv[3] if len(sys.argv) > 3 else ''); os.makedirs(OUTD, exist_ok=True)
temps = sys.argv[2].split(',') if len(sys.argv) > 2 else ['320', '348']
doms = sorted(f.replace('mdcath_dataset_', '').replace('.h5', '') for f in os.listdir(RAW + '/mdcath') if f.endswith('.h5'))
PR = sys.argv[4] if len(sys.argv) > 4 else ROOT + '/results/mdcath_spread_priority.json'
if os.path.exists(PR):
    pr = [d for d in json.load(open(PR)) if d in set(doms)]; doms = pr + [d for d in doms if d not in set(pr)]
def atom_idx(h):
    txt = h['pdbProteinAtoms'][()]; txt = txt.decode() if isinstance(txt, bytes) else str(txt)
    res = {}; order = []; k = 0
    for line in txt.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            key = (line[21], line[22:27]); an = line[12:16].strip()
            if key not in res: res[key] = {}; order.append(key)
            if an in ('N', 'CA', 'C', 'O') and an not in res[key]: res[key][an] = k
            k += 1
    keep = [r for r in order if all(a in res[r] for a in ('N', 'CA', 'C', 'O'))]
    return np.array([[res[r][a] for a in ('N', 'CA', 'C', 'O')] for r in keep])
def dssp_mode(h, T):
    D = np.concatenate([np.asarray(h[T][r]['dssp'][()]) for r in sorted(h[T].keys(), key=int)])
    D = np.vectorize(lambda x: x.decode() if isinstance(x, bytes) else str(x))(D)
    out = []
    for j in range(D.shape[1]):
        v, c = np.unique(D[:, j], return_counts=True); out.append(v[np.argmax(c)])
    return np.array(out)
t0 = time.time(); n = 0
for d in doms:
    out = f'{OUTD}/{d}.npz'
    if os.path.exists(out):
        try:
            zz = np.load(out)
            if 'ca0' in zz.files and 'dssp0' in zz.files: continue
            old = {k: zz[k] for k in zz.files}
            with h5py.File(RAW + f'/mdcath/mdcath_dataset_{d}.h5', 'r') as f:
                h = f[d]; idx = atom_idx(h)
                if 'ca0' not in old: old['ca0'] = h[temps[0]]['0']['coords'][0][idx[:, 1]]
                old['dssp0'] = dssp_mode(h, temps[0])
            np.savez(out, **old); continue
        except Exception as e:
            print('BACKFILL ERR', d, repr(e)[:150], flush=True); continue
    try:
        with h5py.File(RAW + f'/mdcath/mdcath_dataset_{d}.h5', 'r') as f:
            h = f[d]; idx = atom_idx(h); res = {}
            for T in temps:
                Xs, reps = [], []
                for r in sorted(h[T].keys(), key=int):
                    C = h[T][r]['coords'][::3][:, idx.ravel()].reshape(-1, len(idx), 4, 3)
                    X = fs.coordinates_from_backbone(C[:, :, 1], C[:, :, 2], C[:, :, 3]); Xs.append(X); reps.append(len(X))
                XA = np.concatenate(Xs); E = torch.as_tensor(XA)
                xbar = fs.mean_point(E[:600]); Ea = fs.align_to_reference(E, xbar)
                dev = ((Ea - xbar[None]) ** 2).sum((-1, -2)).numpy()
                cut = np.cumsum([0] + reps)
                res[T + '_spread'] = dev.mean(0); res[T + '_spread_rep'] = np.array([dev[cut[i]:cut[i + 1]].mean(0) for i in range(len(reps))])
                res[T + '_xbar'] = xbar.numpy(); res[T + '_nframes'] = len(XA)
                if 'ca0' not in res: res['ca0'] = h[T]['0']['coords'][0][idx[:, 1]]; res['dssp0'] = dssp_mode(h, T)
        np.savez(out, **res); n += 1
        if n % 25 == 0: print(time.strftime('%H:%M:%S'), n, round((time.time() - t0) / n, 1), 's/dom', flush=True)
    except Exception as e:
        print('ERR', d, repr(e)[:200], flush=True)
print('FINISHED', n, flush=True)
open(ROOT + f'/results/logs/mdcath_spread{sys.argv[3] if len(sys.argv) > 3 else ""}.done', 'w').write(time.strftime('%Y-%m-%d %H:%M:%S'))
