"""(prereg round37 R7) 접힘 대 풀림.
주: Majewski 빠른 접힘(CHARMM22*/TIP3P 350 K, 적응 표본 추출 — MSM 가중 없음, 사유는 사후 필드).
보조: mdCATH 413·450 K 같은 실행의 두 무리.
바탕: 맞춘 좌표 유클리드 + d_fix(대응 고정, 쌍마다 회전 닫힌 해).
천연 접촉 = 자료 안에서 형성 빈도 ≥ 0.6 인 Cα 쌍(|i−j| ≥ 4, < 8 Å; 천연 거리 = 형성 프레임의 중앙값). Q = 1.2 × 천연 거리 안 비율.
derived/r37_unfold/*.json, results/r37_unfold.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, tarfile, itertools
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
import mdtraj as md
OUT = ROOT + '/derived/r37_unfold'; os.makedirs(OUT, exist_ok=True); MJ = ROOT + '/raw/majewski'; TMP = '/content/mj_tmp'; os.makedirs(TMP, exist_ok=True)

def dmat_fix(X, Y):
    """d_fix 쌍 거리. X: (n, ne, 3, 3), Y: (m, ne, 3, 3). d² = ‖x‖² + ‖y‖² − 2 max_R tr(R Σ_e y_e x_eᵀ)."""
    M = np.einsum('yeij,xekj->xyik', Y, X); U, S_, Vt = np.linalg.svd(M); dt = np.sign(np.linalg.det(U @ Vt)); S_[..., -1] *= dt
    nx = (X ** 2).sum((1, 2, 3)); ny = (Y ** 2).sum((1, 2, 3)); return np.sqrt(np.maximum(nx[:, None] + ny[None] - 2 * S_.sum(-1), 0))
def edist_fix(A4, B4):
    dab, daa, dbb = dmat_fix(A4, B4), dmat_fix(A4, A4), dmat_fix(B4, B4)
    off = lambda Dm: Dm[~np.eye(len(Dm), dtype=bool)].mean()
    return float(2 * dab.mean() - off(daa) - off(dbb))
def native_from(CA, stride=20):
    S_ = CA[::stride]; ii, jj = np.triu_indices(CA.shape[1], 4); dd = np.linalg.norm(S_[:, ii] - S_[:, jj], axis=-1)
    fr = (dd < 8.0).mean(0); keep = fr >= 0.6; nd = np.array([np.median(dd[dd[:, k] < 8.0, k]) for k in np.nonzero(keep)[0]])
    return ii[keep], jj[keep], nd
def qval(CA, ii, jj, nd):
    return (np.linalg.norm(CA[:, ii] - CA[:, jj], axis=-1) < 1.2 * nd).mean(1)
def delta2(G, ne, rng, N=100, rep=10):
    """G: {(상태, 무리): (평평한 맞춘 좌표, 4차원 맞춘 좌표)}. 교차 = (f,g) 대 (u,1−g), 바닥 = 같은 상태의 무리 0 대 1."""
    pk = lambda A: rng.choice(len(A), min(N, len(A)), replace=False)
    out = {}
    for nm, fn in (('aligned', lambda a, b: edist(a[0], b[0])), ('dfix', lambda a, b: edist_fix(a[1], b[1]))):
        vals = []
        for _ in range(rep):
            sel = {}
            for k, v in G.items():
                i = pk(v[0]); sel[k] = (v[0][i], v[1][i])
            cr = np.mean([fn(sel[('f', g)], sel[('u', 1 - g)]) for g in (0, 1)]); ff = fn(sel[('f', 0)], sel[('f', 1)]); uu = fn(sel[('u', 0)], sel[('u', 1)])
            vals.append((cr - (ff + uu) / 2) / np.sqrt(ne))
        out[nm] = dict(d2n=float(np.mean(vals)), sd=float(np.std(vals)))
    return out

# ---- 주: Majewski ----
for prot in ('chignolin', 'trpcage', 'villin', 'bba', 'homeodomain'):
    if os.path.exists(f'{OUT}/mj_{prot}.json'): continue
    t0 = time.time(); tgz = f'{MJ}/{prot}_trajectories.tar.gz'; top = None; trajs = []; pend = []
    try:
        with tarfile.open(tgz, 'r|gz') as tf:
            for m in tf:
                base = os.path.basename(m.name)
                if not m.isfile() or base.startswith('.'): continue
                if base == 'filtered.pdb' and top is None:
                    p = f'{TMP}/{prot}_top.pdb'; open(p, 'wb').write(tf.extractfile(m).read()); top = md.load(p).topology; trajs += pend; pend = []
                elif base.endswith('.xtc'):
                    p = f'{TMP}/{prot}_{len(trajs) + len(pend)}.xtc'; open(p, 'wb').write(tf.extractfile(m).read())
                    (pend if top is None else trajs).append(p)
                    if len(trajs) >= 400: break
        sel = [a.index for a in top.atoms if a.name in ('CA', 'C', 'O') and a.residue.is_protein]
        X, gid = [], []; rr = np.random.default_rng(7)
        for p in trajs:
            try: t = md.load(p, top=top, atom_indices=sel, stride=10)
            except Exception: continue
            finally:
                if os.path.exists(p): os.remove(p)
            X.append(t.xyz.reshape(t.n_frames, -1, 3, 3) * 10); gid.append(np.full(t.n_frames, rr.integers(0, 2)))
        X = np.concatenate(X); gid = np.concatenate(gid); CA = X[:, :, 0]
        ii, jj, nd = native_from(CA); Q = qval(CA, ii, jj, nd)
        E = frames(X); ref = fs.mean_point(E[np.nonzero(Q >= 0.75)[0][:600]]); Ya = fs.align_to_reference(E, ref).numpy(); ne = Ya.shape[1]; F = Ya.reshape(len(Ya), -1)
        G = {}
        for s, msk in (('f', Q >= 0.75), ('u', Q < 0.4)):
            for g in (0, 1):
                ix = np.nonzero(msk & (gid == g))[0]
                if len(ix) < 40: raise RuntimeError(f'too few {s}{g}: {len(ix)}')
                G[(s, g)] = (F[ix], Ya[ix])
        res = dict(protein=prot, n_traj=len(trajs), n_frames=int(len(X)), ne=int(ne), n_native=int(len(ii)), frac_folded=float((Q >= 0.75).mean()), frac_unfolded=float((Q < 0.4).mean()),
                   counts={f'{s}{g}': int(len(G[(s, g)][0])) for s, g in G}, sec=round(time.time() - t0, 1), note='MSM 가중 없음(적응 표본 추출), 궤적을 무작위 두 무리로, 400 궤적까지, 10 프레임 간격')
        res.update(delta2(G, ne, np.random.default_rng(37)))
        json.dump(res, open(f'{OUT}/mj_{prot}.json', 'w')); print(time.strftime('%H:%M:%S'), prot, json.dumps(res)[:400], flush=True)
    except Exception as e:
        import traceback; print('ERR', prot, repr(e)[:200], traceback.format_exc()[-500:], flush=True)
    for f in glob.glob(f'{TMP}/{prot}_*'): os.remove(f)

# ---- 보조: mdCATH 413·450 K ----
for d in st['domains']:
    if os.path.exists(f'{OUT}/md_{d}.json'): continue
    try:
        D, _ = load(d); res = dict(dom=d)
        for T in ('413', '450'):
            Qd = D[T]['Q']; rp = D[T]['rep']; fold = Qd >= QT; unf = Qd < 0.4
            reps_f = [r for r in np.unique(rp) if (fold & (rp == r)).sum() >= 40]; reps_u = [r for r in np.unique(rp) if (unf & (rp == r)).sum() >= 40]
            if len(reps_f) < 2 or len(reps_u) < 2: continue
            E = frames(D[T]['X']); ref = fs.mean_point(E[np.nonzero(fold)[0][:600]]); Ya = fs.align_to_reference(E, ref).numpy(); ne = Ya.shape[1]; F = Ya.reshape(len(Ya), -1)
            rng = np.random.default_rng(5); pk = lambda ix: rng.choice(ix, min(100, len(ix)), replace=False)
            fi = {r: np.nonzero(fold & (rp == r))[0] for r in reps_f}; ui = {r: np.nonzero(unf & (rp == r))[0] for r in reps_u}
            prs = [(i, j) for i in reps_f for j in reps_u if i != j][:10]
            out = {}
            for nm, fn in (('aligned', lambda a, b: edist(F[a], F[b])), ('dfix', lambda a, b: edist_fix(Ya[a], Ya[b]))):
                cr = np.mean([fn(pk(fi[i]), pk(ui[j])) for i, j in prs])
                ff = np.mean([fn(pk(fi[i]), pk(fi[j])) for i, j in itertools.combinations(reps_f, 2)][:10]); uu = np.mean([fn(pk(ui[i]), pk(ui[j])) for i, j in itertools.combinations(reps_u, 2)][:10])
                out[nm] = float((cr - (ff + uu) / 2) / np.sqrt(ne))
            res[T] = dict(ne=int(ne), n_reps_f=len(reps_f), n_reps_u=len(reps_u), **out)
        if len(res) > 1: json.dump(res, open(f'{OUT}/md_{d}.json', 'w')); print(time.strftime('%H:%M:%S'), d, flush=True)
    except Exception as e: print('ERR', d, repr(e)[:150], flush=True)

S = dict(majewski={}, mdcath={})
for q in glob.glob(OUT + '/mj_*.json'):
    j = json.load(open(q)); S['majewski'][j['protein']] = {k: j[k] for k in ('ne', 'n_traj', 'n_frames', 'frac_folded', 'frac_unfolded', 'aligned', 'dfix')}
M_ = [json.load(open(q)) for q in glob.glob(OUT + '/md_*.json')]
for T in ('413', '450'):
    sb = [m[T] for m in M_ if T in m]
    if sb: S['mdcath'][T] = dict(n=len(sb), aligned_median=round(float(np.median([x['aligned'] for x in sb])), 4), dfix_median=round(float(np.median([x['dfix'] for x in sb])), 4),
                               aligned_iqr=[round(float(x), 4) for x in np.percentile([x['aligned'] for x in sb], [25, 75])])
S['compare'] = {'+130K_full': 0.066, '+130K_second_half': 0.0663}
S['written'] = time.strftime('%F %T'); json.dump(S, open(ROOT + '/results/r37_unfold.json', 'w'), indent=1, ensure_ascii=False); print(json.dumps(S, indent=1, ensure_ascii=False))
open(ROOT + '/results/logs/r37_unfold.done', 'w').write(time.strftime('%F %T'))
