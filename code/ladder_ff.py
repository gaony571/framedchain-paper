
"""사다리 역장 행 C5 (2026-10-18; 예측 prereg round10_ff). SM 대 BioEmu CATH MD, 같은 CATH 도메인. 도메인마다 derived/ladder_ff/<dom>.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, time, glob, zipfile, difflib, shutil
ROOT = sys.argv[1]; ONLY = sys.argv[2].split(',') if len(sys.argv) > 2 else None
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, torch, mdtraj as md, framespace as fs, ot
torch.set_num_threads(1)
src = open(ROOT + '/code/ladder_model_error.py').read()
exec(src[src.index('AA = dict('):src.index('chains = ONLY or')])
OUTD = ROOT + '/derived/ladder_ff'; os.makedirs(OUTD, exist_ok=True); TMP = '/content/tmp_ff'; os.makedirs(TMP, exist_ok=True)
ZB = {k: zipfile.ZipFile(f'{ROOT}/raw/bioemu_cath/{k}.zip') for k in ('MSR_cath2', 'ONE_cath1')}; bio = {}
for k, Z in ZB.items():
    for n in Z.namelist():
        p = n.split('/')
        if len(p) >= 2 and p[1]: bio.setdefault(p[1].split('_', 1)[1], (k, p[1]))
smd = sorted(os.path.basename(p)[:-4] for p in glob.glob(ROOT + '/raw/mfm_sm/*.npz')); doms = ONLY or [d for d in smd if d in bio]
def sm_backbone(d, every=20):
    z = np.load(f'{ROOT}/raw/mfm_sm/{d}.npz', allow_pickle=True); meta = json.loads(str(z['meta'])); L = z['atom_lines']
    res, order, seq = {}, [], []
    for i, l in enumerate(L):
        rn, an, key = l[17:20].strip(), l[12:16].strip(), l[21:27]
        if rn in ('ACE', 'NME', 'NHE'): continue
        if key not in res: res[key] = {}; order.append((key, rn))
        if an in ('CA', 'C', 'O') and an not in res[key]: res[key][an] = i
    keep = [(k, rn) for k, rn in order if all(a in res[k] for a in ('CA', 'C', 'O'))]
    idx = np.array([[res[k][a] for a in ('CA', 'C', 'O')] for k, _ in keep]); seq = ''.join(AA.get(rn, 'X') for _, rn in keep)
    P = z['positions']; F = len(P); per = F // meta['sims']; run = np.arange(F) // per
    sel = np.arange(0, F, every); X = P[sel][:, idx.ravel()].reshape(len(sel), len(idx), 3, 3).astype(np.float64)
    return X, seq, run[sel], meta
for d in doms:
    out = f'{OUTD}/{d}.json'
    if os.path.exists(out) and not ONLY: continue
    t0 = time.time(); rng = np.random.default_rng(abs(hash(d)) % (2 ** 32)); tmp = f'{TMP}/{d}'; os.makedirs(tmp, exist_ok=True)
    try:
        Xs, seq_s, run_s, meta = sm_backbone(d)
        k, sysn = bio[d]; Z = ZB[k]; pre = f'{k}/{sysn}/'
        Z.extract(pre + 'topology.pdb', tmp); top = md.load(f'{tmp}/{pre}topology.pdb')
        Xb, run_b = [], []
        for i, n in enumerate(sorted(n for n in Z.namelist() if n.startswith(pre) and n.endswith('.xtc'))):
            Z.extract(n, tmp); tr = md.load(f'{tmp}/{n}', top=top); X, seq_b = backbone(tr); Xb.append(X.astype(np.float64)); run_b.append(np.full(len(X), i))
        Xb, run_b = np.concatenate(Xb), np.concatenate(run_b)
        sm_ = difflib.SequenceMatcher(None, seq_s, seq_b, autojunk=False); ia, ib = [], []
        for blk in sm_.get_matching_blocks(): ia += range(blk.a, blk.a + blk.size); ib += range(blk.b, blk.b + blk.size)
        ia, ib = np.array(ia), np.array(ib)
        if len(ia) < 20: raise RuntimeError('few matched residues')
        Es, Eb = frames(Xs[:, ia]), frames(Xb[:, ib]); ref = fs.mean_point(Es[::max(1, len(Es) // 600)][:600])
        F = lambda E: np.asarray(fs.align_to_reference(E, ref)).reshape(len(E), -1); S, B = F(Es), F(Eb)
        N = 250; pick = lambda A, m: A[rng.choice(len(A), min(m, len(A)), replace=False)]
        rs, rb = np.unique(run_s), np.unique(run_b); hs, hb = rng.permutation(rs), rng.permutation(rb)
        S1, S2 = S[np.isin(run_s, hs[: len(hs) // 2])], S[np.isin(run_s, hs[len(hs) // 2:])]; B1, B2 = B[np.isin(run_b, hb[: len(hb) // 2])], B[np.isin(run_b, hb[len(hb) // 2:])]
        res = dict(dom=d, n_res_matched=int(len(ia)), len_sm=len(seq_s), len_bio=len(seq_b), frames_sm=int(len(S)), frames_bio=int(len(B)), runs_sm=int(len(rs)), runs_bio=int(len(rb)), bio_source=k,
                   floor_sm=stats(pick(S1, N), pick(S2, N), rng), floor_bio=stats(pick(B1, N), pick(B2, N), rng),
                   cross=stats(pick(S1, N), pick(B1, N), rng), cross2=stats(pick(S2, N), pick(B2, N), rng))
        res['sec'] = round(time.time() - t0, 1); json.dump(res, open(out, 'w')); print(time.strftime('%H:%M:%S'), d, res['sec'], flush=True)
    except Exception as e:
        import traceback; print('ERR', d, repr(e)[:200], traceback.format_exc()[-400:], flush=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
if not ONLY: open(ROOT + '/results/logs/ladder_ff.done', 'w').write(time.strftime('%F %T'))
