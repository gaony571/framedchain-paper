
"""C4 치환 축 (prereg round14; 계획 4절 C4). MegaSim 변이체 대 야생형(ff99SB-disp 궤적만), Tsuboyama ddG.
derived/c4_sub/<WT>.json (이어 하기), 끝에 results/c4_substitution.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, time, zipfile, glob, shutil
ROOT = sys.argv[1]; ONLY = sys.argv[2].split(',') if len(sys.argv) > 2 and sys.argv[2] != 'all' else None
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, pandas as pd, torch, mdtraj as md, framespace as fs
from scipy.stats import spearmanr
torch.set_num_threads(1)
QT = 0.75; OUTD = ROOT + '/derived/c4_sub'; os.makedirs(OUTD, exist_ok=True); TMP = '/content/tmp_c4'; os.makedirs(TMP, exist_ok=True)
def zpath(local, drive): return local if os.path.exists(local) else drive
ZM = zipfile.ZipFile(zpath('/content/megasim_merge.zip', ROOT + '/raw/megasim/MSR_megasim_merge.zip'))
ZU = zipfile.ZipFile(zpath('/content/megasim_mut.zip', ROOT + '/raw/megasim/MSR_megasim_mutants_disp_allatom.zip'))
NM, NU = set(ZM.namelist()), ZU.namelist()
MUTS = {}
for n in NU:
    if n.endswith('trj_mutant_folded.xtc'): s = n.split('/')[1]; MUTS.setdefault(s.split('__')[0], []).append(s)
T = pd.read_csv(zipfile.ZipFile(ROOT + '/raw/tsuboyama/Processed_K50_dG_datasets.zip').open('Processed_K50_dG_datasets/Tsuboyama2023_Dataset2_Dataset3_20230416.csv'), usecols=['aa_seq', 'mut_type', 'WT_name', 'dG_ML', 'ddG_ML'], low_memory=False)
T['ddg'] = pd.to_numeric(T.ddG_ML, errors='coerce'); T['dg'] = pd.to_numeric(T.dG_ML, errors='coerce'); T.loc[T.dg.isna(), 'ddg'] = np.nan
WTSEQ = T[T.mut_type == 'wt'].drop_duplicates('WT_name').set_index('WT_name').aa_seq.to_dict()
DDG = T[T.mut_type != 'wt'].groupby(['WT_name', 'mut_type']).ddg.mean().to_dict()
def tsu_name(w):
    return w.split('_', 1)[0] + '.pdb_' + w.split('_', 1)[1] if '_' in w else w + '.pdb'
def ext(Z, member):
    p = TMP + '/' + member.replace('/', '__')
    with Z.open(member) as f, open(p, 'wb') as g: shutil.copyfileobj(f, g)
    return p
def bb_idx(t):
    out = []
    for r in t.residues:
        if not r.is_protein: continue
        a = {at.name: at.index for at in r.atoms}; out.append([a.get('N', -1), a.get('CA', -1), a.get('C', -1), a.get('O', -1)])
    return np.array(out)
def seq_of(t): return ''.join((r.code or 'X') for r in t.residues if r.is_protein)
def cost(A, B): return np.maximum((A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2 * A @ B.T, 0)
def edist(A, B):
    dab = np.sqrt(cost(A, B)); daa = np.sqrt(cost(A, A)); dbb = np.sqrt(cost(B, B)); n, m = len(A), len(B)
    return float(2 * dab.mean() - daa.sum() / (n * (n - 1)) - dbb.sum() / (m * (m - 1)))
def frames(xyz):
    return fs.coordinates_from_backbone(torch.as_tensor(xyz[:, :, 1], dtype=torch.float64), torch.as_tensor(xyz[:, :, 2], dtype=torch.float64), torch.as_tensor(xyz[:, :, 3], dtype=torch.float64))
def even(n, k): return np.arange(n) if n <= k else np.linspace(0, n - 1, k).round().astype(int)
def resid(y, x):
    X = np.c_[np.ones(len(x)), x]; return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
def run(w):
    t0 = time.time(); base = f'MSR_megasim_merge/{w}'; dflt = json.loads(ZM.read(base + '/dataset.json'))['force_field']
    ref = md.load(ext(ZM, base + '/reference.pdb')); topp = ext(ZM, base + '/topology.pdb') if base + '/topology.pdb' in NM else None
    top = md.load(topp).topology if topp else ref.topology
    bref, btop = bb_idx(ref.topology), bb_idx(top); assert len(bref) == len(btop)
    ok = (bref >= 0).all(1) & (btop >= 0).all(1); CAr = ref.xyz[0][bref[ok, 1]] * 10; nres = int(ok.sum()); ne = nres - 1
    Dn = np.linalg.norm(CAr[:, None] - CAr[None], axis=-1); ii, jj = np.triu_indices(nres, 4); m_ = Dn[ii, jj] < 8.0; ii, jj, dn = ii[m_], jj[m_], Dn[ii, jj][m_]
    Qof = lambda xyz: (np.linalg.norm(xyz[:, ii, 1] - xyz[:, jj, 1], axis=-1) < 1.2 * dn).mean(1)
    tr_ = {}
    for x in sorted(n for n in NM if n.startswith(base + '/trajs/') and n.endswith('.xtc')):
        jn = x.replace('.cmprsd.xtc', '.json'); ff = json.loads(ZM.read(jn))['force_field'] if jn in NM else dflt
        if 'disp' not in ff: continue
        p = ext(ZM, x); t = md.load(p, top=top, atom_indices=btop[ok].ravel()); os.remove(p)
        xyz = t.xyz.reshape(t.n_frames, -1, 4, 3) * 10; tr_[x.split('/')[-1]] = (xyz, Qof(xyz))
    F = [k for k, (xyz, q) in tr_.items() if (q >= QT).sum() >= 20]
    res = dict(wt=w, n_disp_traj=len(tr_), n_folded_traj=len(F), nres=nres)
    if len(F) < 3: res['skip'] = 'folded trajectories < 3'; return res
    G1, G3 = F[0::3], F[2::3]; G2 = [k for k in tr_ if k not in G1 and k not in G3]
    fold = lambda k: tr_[k][0][tr_[k][1] >= QT]
    X1 = np.concatenate([fold(k) for k in G1]); rng = np.random.default_rng(0)
    xbar = fs.mean_point(frames(X1[rng.choice(len(X1), min(600, len(X1)), replace=False)]))
    def Yof(xyz): return fs.align_to_reference(frames(xyz), xbar).numpy()
    xb = xbar.numpy()
    def spread(Y4): return ((Y4 - xb[None]) ** 2).sum((-1, -2)).mean(0)
    Y1 = Yof(X1); s1 = spread(Y1); C1 = Y1.reshape(len(Y1), -1)[even(len(Y1), 300)]
    if len(G1) >= 2:
        h = [np.concatenate([fold(k) for k in G1[0::2]]), np.concatenate([fold(k) for k in G1[1::2]])]; s1h = [spread(Yof(z)) for z in h]
    else:
        n1 = len(Y1) // 2; s1h = [spread(Y1[:n1]), spread(Y1[n1:])]
    level = np.log(spread(Yof(np.concatenate([fold(k) for k in G3]))))
    G2f = np.concatenate([fold(k) for k in G2 if (tr_[k][1] >= QT).sum() > 0]) if any((tr_[k][1] >= QT).sum() > 0 for k in G2) else None
    G2m = [tr_[k][0][(tr_[k][1] >= 0.5) & (tr_[k][1] < QT)] for k in G2]; G2m = np.concatenate(G2m) if G2m else np.zeros((0,))
    u = None
    if G2f is not None and len(G2f) >= 20 and len(G2m) >= 20:
        u = np.log(spread(Yof(G2m[even(len(G2m), 600)]))) - np.log(spread(Yof(G2f[even(len(G2f), 600)]))); ur = resid(u, level); res['u_kind'] = 'mid_Q'
    elif G2f is not None and len(G2f) >= 80:
        q2 = np.concatenate([tr_[k][1][tr_[k][1] >= QT] for k in G2 if (tr_[k][1] >= QT).sum() > 0]); lo = q2 <= np.quantile(q2, 0.25)
        u = np.log(spread(Yof(G2f[lo]))) - np.log(spread(Yof(G2f[~lo][even(int((~lo).sum()), 600)]))); ur = resid(u, level); res['u_kind'] = 'lowQ_quartile'
    shifts = rng.integers(5, max(6, ne - 5), 40) if ne > 12 else np.arange(1, 41)
    fl = []
    for k in F:
        if k in G1: continue
        Yk = Yof(fold(k)[even(len(fold(k)), 150)]); fl.append(edist(Yk.reshape(len(Yk), -1), C1))
    res.update(G=[G1, G2, G3], floor=float(np.median(fl)), floor_all=fl, n_mid=int(len(G2m)), has_dir=u is not None, n_g1_frames=int(len(Y1)))
    tn = tsu_name(w); wseq = WTSEQ.get(tn); tseq = seq_of(top); okseq = ''.join(np.array(list(tseq))[ok]) if len(tseq) == len(ok) else None
    off = okseq.find(wseq) if (wseq and okseq) else -1; res.update(tsu=tn, offset=int(off))
    mid = (CAr[:-1] + CAr[1:]) / 2; cn = (Dn < 10).sum(1) - 1
    M = []
    for s in sorted(MUTS.get(w, [])):
        mt = s.split('__')[1]; key = mt.replace('_', ':'); r = dict(m=s, ddg=DDG.get((tn, key), None), single=':' not in key)
        r['ddg'] = None if r['ddg'] is None or not np.isfinite(r['ddg']) else float(r['ddg'])
        try:
            tp = ext(ZU, f'MSR_megasim_mutants_disp_allatom/{s}/topology.pdb'); mtop = md.load(tp).topology; bm = bb_idx(mtop); os.remove(tp)
            if len(bm) != len(ok) or not (bm[ok] >= 0).all(): r['err'] = 'residues'; M.append(r); continue
            xp = ext(ZU, f'MSR_megasim_mutants_disp_allatom/{s}/trajs/trj_mutant_folded.xtc'); t = md.load(xp, top=mtop, atom_indices=bm[ok].ravel()); os.remove(xp)
            xyz = t.xyz.reshape(t.n_frames, -1, 4, 3) * 10; q = Qof(xyz); fm = q >= QT; r.update(n_frames=int(len(q)), n_folded=int(fm.sum()), frac_folded=float(fm.mean()))
            sites = []
            for tok in mt.split('_'):
                a, pos = tok[0], int(tok[1:-1]); ix = off + pos - 1 if off >= 0 else -1
                if 0 <= ix < nres: sites.append(ix); r.setdefault('site_match', []).append(bool(okseq[ix] == a))
            r['sites'] = sites
            if fm.sum() < 20: M.append(r); continue
            Ym = Yof(xyz[fm][even(int(fm.sum()), 150)]); r['ed'] = edist(Ym.reshape(len(Ym), -1), C1)
            sm = spread(Ym); D = np.log(sm) - np.log(s1); r['dnorm'] = float((D ** 2).mean())
            nh = len(Ym) // 2; D1 = np.log(spread(Ym[:nh])) - np.log(s1h[0]); D2 = np.log(spread(Ym[nh:])) - np.log(s1h[1]); r['cross'] = float((D1 * D2).mean())
            if u is not None:
                Dr = resid(D, level); r['dir'] = float(np.corrcoef(Dr, ur)[0, 1]); r['dir_null'] = [float(np.corrcoef(Dr, np.roll(ur, int(k)))[0, 1]) for k in shifts]
            if sites:
                dist = np.min(np.linalg.norm(mid[:, None] - CAr[sites][None], axis=-1), 1); r['decay'] = float(spearmanr(np.abs(D), dist)[0]); r['burial'] = float(np.mean(cn[sites]))
        except Exception as e:
            r['err'] = repr(e)[:150]
        M.append(r)
    res['mutants'] = M; res['sec'] = round(time.time() - t0, 1); return res
todo = ONLY if ONLY else sorted(MUTS)
for w in todo:
    p = f'{OUTD}/{w}.json'
    if os.path.exists(p) and not ONLY: continue
    try:
        r = run(w); json.dump(r, open(p, 'w')); print(time.strftime('%H:%M:%S'), w, len(r.get('mutants', [])), r.get('sec'), r.get('skip', ''), flush=True)
    except Exception as e:
        import traceback; print('ERR', w, repr(e)[:200], traceback.format_exc()[-600:], flush=True)
if ONLY: sys.exit()
# ---- 요약 ----
J = [json.load(open(q)) for q in glob.glob(OUTD + '/*.json')]; J = [j for j in J if 'mutants' in j]; rng = np.random.default_rng(1)
rows = [dict(wt=j['wt'], floor=j['floor'], **{k: v for k, v in m.items() if k not in ('dir_null', 'sites', 'site_match')}) for j in J for m in j['mutants']]
Dm = pd.DataFrame(rows); Dm['absd'] = Dm.ed / Dm.floor; S = dict(n_wt=len(J), n_mut=len(Dm), n_with_ddg=int(Dm.ddg.notna().sum()), n_folded20=int((Dm.n_folded >= 20).sum()) if 'n_folded' in Dm else 0)
ok_ = Dm[Dm.single & Dm.ddg.notna() & Dm.ed.notna()]
def wt_rho(df, a, b, minn=10):
    out = {}
    for w, g in df.groupby('wt'):
        if len(g) >= minn and g[a].nunique() > 2: out[w] = spearmanr(g[a], g[b])[0]
    return pd.Series(out).dropna()
def boot(s):
    v = s.values; b = [np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(2000)]; return dict(n_wt=len(v), median=round(float(np.median(v)), 4), CI95=[round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])], frac_pos=round(float((v > 0).mean()), 3))
ok_ = ok_.assign(neg_ddg=-ok_.ddg)
S['test1_ed'] = boot(wt_rho(ok_, 'absd', 'neg_ddg')); S['test1_cross'] = boot(wt_rho(ok_[ok_.cross.notna()], 'cross', 'neg_ddg'))
rr = []
for w, g in ok_.groupby('wt'):
    if len(g) < 10 or g.burial.isna().any(): continue
    Z = np.c_[g.n_folded.rank(), g.burial.rank()]; rr.append((w, spearmanr(resid(g.absd.rank().values, Z), resid(g.neg_ddg.rank().values, Z))[0]))
S['test1_ed_resid_nfolded_burial'] = boot(pd.Series(dict(rr)).dropna())
S['extra_fracfolded_vs_ddg'] = boot(wt_rho(Dm[Dm.single & Dm.ddg.notna() & Dm.frac_folded.notna()], 'frac_folded', 'ddg'))
# test2
obs, nulls = [], []
for j in J:
    ms = [m for m in j['mutants'] if m.get('dir') is not None and m.get('ddg') is not None and m['ddg'] < -1 and m.get('single')]
    if len(ms) < 5: continue
    obs.append(np.mean([m['dir'] for m in ms])); nulls.append(np.mean([m['dir_null'] for m in ms], 0))
if obs:
    obs, nulls = np.array(obs), np.array(nulls); nm = nulls.mean(0); bo = [obs[rng.integers(0, len(obs), len(obs))].mean() for _ in range(2000)]
    S['test2_direction'] = dict(n_wt=len(obs), mean_obs=round(float(obs.mean()), 4), CI95=[round(float(x), 4) for x in np.percentile(bo, [2.5, 97.5])], null_mean=round(float(nm.mean()), 4), null_q99=round(float(np.quantile(nm, 0.99)), 4), p_shift=float((nm >= obs.mean()).mean()), frac_wt_pos=round(float((obs > 0).mean()), 3))
dc = Dm[Dm.decay.notna()] if 'decay' in Dm else Dm.iloc[:0]
S['test3_decay'] = dict(n=int(len(dc)), median=round(float(dc.decay.median()), 4), frac_neg=round(float((dc.decay < 0).mean()), 3), wt_median=boot(dc.groupby('wt').decay.median()))
S['site_match_rate'] = float(np.mean([x for j in J for m in j['mutants'] for x in m.get('site_match', [])])) if J else None
S['written'] = time.strftime('%F %T'); json.dump(S, open(ROOT + '/results/c4_substitution.json', 'w'), indent=1); print(json.dumps(S, indent=1))
open(ROOT + '/results/logs/c4_sub.done', 'w').write(time.strftime('%F %T'))
