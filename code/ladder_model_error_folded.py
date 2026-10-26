
"""사다리 모델 오차 행 (2026-10-18; 예측 prereg round9_model_error). 사슬마다 derived/ladder_model_error/<chain>.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, time, glob, zipfile, difflib, shutil, tempfile
ROOT = sys.argv[1]; ONLY = sys.argv[2].split(',') if len(sys.argv) > 2 else None
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, torch, mdtraj as md, framespace as fs, ot
torch.set_num_threads(1)
OUTD = ROOT + '/derived/ladder_model_error_folded'; os.makedirs(OUTD, exist_ok=True); TMP = '/content/tmp_me'; os.makedirs(TMP, exist_ok=True)
MODELS = ['alphaflow_md_base_202402', 'alphaflow_md_templates_base_202402', 'esmflow_md_base_202402', 'esmflow_md_templates_base_202402']
ZM = {m: zipfile.ZipFile(f'/content/me_raw/{m}.zip') for m in MODELS}
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', HID='H', HIE='H', HIP='H', HSD='H', HSE='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
def backbone(traj):
    res, seq = [], []
    for r in traj.topology.residues:
        if not r.is_protein: continue
        a = {at.name: at.index for at in r.atoms}
        if all(k in a for k in ('CA', 'C', 'O')): res.append([a['CA'], a['C'], a['O']]); seq.append(AA.get(r.name, 'X'))
    idx = np.array(res); X = traj.xyz[:, idx.ravel()].reshape(traj.n_frames, len(idx), 3, 3) * 10.0
    return X, ''.join(seq)
def frames(X):
    E = fs.coordinates_from_backbone(torch.as_tensor(np.ascontiguousarray(X[:, :, 0], np.float64)), torch.as_tensor(np.ascontiguousarray(X[:, :, 1], np.float64)), torch.as_tensor(np.ascontiguousarray(X[:, :, 2], np.float64)))
    return torch.as_tensor(np.asarray(E, np.float64)) if not torch.is_tensor(E) else E.double()
def cost(A, B): return np.maximum((A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2 * A @ B.T, 0)
def w2(A, B, center=False):
    if center: A = A - A.mean(0); B = B - B.mean(0)
    return float(ot.emd2(np.full(len(A), 1 / len(A)), np.full(len(B), 1 / len(B)), cost(A, B), numItermax=500000))
def dw2(A, B, rng, center=False):
    ia, ib = rng.permutation(len(A)), rng.permutation(len(B)); a1, a2 = ia[: len(A) // 2], ia[len(A) // 2:]; b1, b2 = ib[: len(B) // 2], ib[len(B) // 2:]
    return 0.5 * (w2(A[a1], B[b2], center) + w2(A[a2], B[b1], center)) - 0.5 * w2(A[a1], A[a2], center) - 0.5 * w2(B[b1], B[b2], center)
def edist(A, B):
    dab, daa, dbb = np.sqrt(cost(A, B)), np.sqrt(cost(A, A)), np.sqrt(cost(B, B)); n, m = len(A), len(B)
    return float(2 * dab.mean() - daa.sum() / (n * (n - 1)) - dbb.sum() / (m * (m - 1)))
def stats(A, B, rng):
    return dict(dw2=dw2(A, B, rng), dw2c=dw2(A, B, rng, True), ed=edist(A, B), mean2=float(((A.mean(0) - B.mean(0)) ** 2).sum()),
                spreadA=float(((A - A.mean(0)) ** 2).sum(1).mean()), spreadB=float(((B - B.mean(0)) ** 2).sum(1).mean()))
chains = ONLY or sorted(os.path.basename(p).replace('_analysis.zip', '') for p in glob.glob('/content/me_raw/af82/*.zip'))
NPER = 250
for c in chains:
    out = f'{OUTD}/{c}.json'
    if os.path.exists(out) and not ONLY: continue
    t0 = time.time(); rng = np.random.default_rng(abs(hash(c)) % (2 ** 32)); d = f'{TMP}/{c}'; os.makedirs(d, exist_ok=True)
    try:
        Za = zipfile.ZipFile(f'/content/me_raw/af82/{c}_analysis.zip')
        for n in Za.namelist():
            if n.endswith(('.pdb', '.xtc')): Za.extract(n, d)
        top = md.load(f'{d}/{c}.pdb'); Xmd, seq_md = [], None
        for r in (1, 2, 3):
            tr = md.load(f'{d}/{c}_R{r}.xtc', top=top); X, seq_md = backbone(tr); Xmd.append(X)
        res = dict(chain=c, n_res=len(seq_md))
        mods = {}
        for m in MODELS:
            nm = f'{m}/{c}.pdb'
            if nm not in ZM[m].namelist(): continue
            p = f'{d}/{m}.pdb'; open(p, 'wb').write(ZM[m].read(nm)); trm = md.load(p); X, sq = backbone(trm); mods[m] = (X, sq)
        # 대응: 서열 맞춤(같은 서열이면 항등)
        def match(sa, sb):
            sm = difflib.SequenceMatcher(None, sa, sb, autojunk=False); ia, ib = [], []
            for blk in sm.get_matching_blocks():
                ia += range(blk.a, blk.a + blk.size); ib += range(blk.b, blk.b + blk.size)
            return np.array(ia), np.array(ib)
        common = np.arange(len(seq_md))
        for m, (X, sq) in mods.items():
            ia, ib = match(seq_md, sq); common = np.intersect1d(common, ia)
        res['n_common_res'] = int(len(common))
        if len(common) < 20: raise RuntimeError('few common residues')
        CA0 = backbone(top)[0][0, common, 0]; Dn = np.linalg.norm(CA0[:, None] - CA0[None], axis=-1); qi, qj = np.triu_indices(len(common), 4); mq = Dn[qi, qj] < 8.0; qi, qj, qd = qi[mq], qj[mq], Dn[qi, qj][mq]
        Qf = lambda C: (np.linalg.norm(C[:, qi] - C[:, qj], axis=-1) < 1.2 * qd).mean(1)
        ff = []
        for i_, X in enumerate(Xmd):
            q = Qf(X[:, common, 0]); ff.append(float((q >= 0.75).mean())); keep = np.nonzero(q >= 0.75)[0]
            if len(keep) < 50: raise RuntimeError('few folded MD frames')
            Xmd[i_] = X[keep[np.linspace(0, len(keep) - 1, NPER).round().astype(int)]]
        res['fracfold_md'] = ff
        # 이어진 잔기만(변 정의를 위해 모든 원천에서 같은 잔기 집합)
        Emd = [frames(X[:, common]) for X in Xmd]; ref = fs.mean_point(torch.cat(Emd)[::3][:600])
        F = lambda E: np.asarray(fs.align_to_reference(E, ref)).reshape(len(E), -1)
        MD = [F(E) for E in Emd]; MDall = np.concatenate(MD)
        fl = [stats(MD[i], MD[j], rng) for i, j in ((0, 1), (0, 2), (1, 2))]; res['floor'] = {k: float(np.mean([f[k] for f in fl])) for k in fl[0]}
        for m, (X, sq) in mods.items():
            ia, ib = match(seq_md, sq); mp = dict(zip(ia, ib)); sel = np.array([mp[i] for i in common])
            Xm = X[:, sel]; qm = Qf(Xm[:, :, 0]); fm_ = float((qm >= 0.75).mean()); Xm = Xm[qm >= 0.75]
            if len(Xm) < 30: res[m] = dict(skip='few folded samples', fracfold=fm_); continue
            Xm = Xm[rng.choice(len(Xm), min(len(Xm), NPER), replace=False)]
            Fm = F(frames(Xm)); per = [stats(Fm, MD[i], rng) for i in range(3)]   # 바닥(반복 대 반복)과 같은 짜임: 모델 대 반복 하나씩, 평균
            res[m] = dict(n_samples=int(len(X)), fracfold=fm_, **{k: float(np.mean([q[k] for q in per])) for k in per[0]}, pooled=stats(Fm, MDall[rng.choice(len(MDall), NPER, replace=False)], rng))
        res['sec'] = round(time.time() - t0, 1); json.dump(res, open(out, 'w'))
        print(time.strftime('%H:%M:%S'), c, res['sec'], flush=True)
    except Exception as e:
        import traceback; print('ERR', c, repr(e)[:200], traceback.format_exc()[-500:], flush=True)
    finally:
        shutil.rmtree(d, ignore_errors=True)
if not ONLY: open(ROOT + '/results/logs/ladder_model_error_folded.done', 'w').write(time.strftime('%F %T'))
