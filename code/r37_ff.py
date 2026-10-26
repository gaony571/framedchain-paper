"""(prereg round37 R3′) 같은 단백질의 역장 + 규약 칸: ATLAS (CHARMM36m, 300 K, 반복 3) 대 mdCATH 320 K (CHARMM22*, 반복 5).
mdCATH ∩ ATLAS 44 사슬, 공통 잔기(서열 정렬의 가장 긴 연속 일치 구간), 접힌 프레임(mdCATH 320 K 반복 0 첫 프레임의 천연 접촉, Q ≥ 0.75).
같은 도메인의 mdCATH 온도 곡선(320 대 348/379/413/450 K; 전체·뒤 절반)도 같은 통계로 내서 온도 환산.
derived/r37_ff/<chain>.json, results/r37_ff.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools, zipfile, shutil, difflib
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
import h5py, mdtraj as md, pandas as pd
OUT = ROOT + '/derived/r37_ff'; os.makedirs(OUT, exist_ok=True); TMP = '/content/ff_tmp'; os.makedirs(TMP, exist_ok=True)
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', HSD='H', HSE='H', HSP='H', HIE='H', HID='H', HIP='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
TT = ('320', '348', '379', '413', '450')
def md_atoms(pdbtxt):
    res = {}; order = []
    for i, l in enumerate([l for l in pdbtxt.splitlines() if l.startswith(('ATOM', 'HETATM'))]):
        nm, rn, rs = l[12:16].strip(), l[17:21].strip(), (l[21], l[22:27].strip())
        if rs not in res: res[rs] = {'rn': rn}; order.append(rs)
        res[rs][nm] = i
    keep = [r for r in order if all(a in res[r] for a in ('CA', 'C', 'O'))]
    return ''.join(AA.get(res[r]['rn'], 'X') for r in keep), np.array([[res[r]['CA'], res[r]['C'], res[r]['O']] for r in keep])
def at_atoms(top):
    keep = [r for r in top.residues if r.is_protein and all(any(a.name == n for a in r.atoms) for n in ('CA', 'C', 'O'))]
    idx = np.array([[next(a.index for a in r.atoms if a.name == n) for n in ('CA', 'C', 'O')] for r in keep])
    return ''.join(AA.get(r.name, 'X') for r in keep), idx
def block(s1, s2):
    m = max(difflib.SequenceMatcher(None, s1, s2, autojunk=False).get_matching_blocks(), key=lambda b: b.size); return m.a, m.b, m.size
def stats(runsA, runsB, ne, rng, N=100):
    """runsA/B: 실행별 평평한 맞춘 좌표 목록. Δ²/√변 (교차 − 바닥 평균/2), 평균 이동 Δμ²/변, 로그 총퍼짐 비(B/A)."""
    pk = lambda X: X[rng.choice(len(X), min(N, len(X)), replace=False)]
    prs = [(i, j) for i in range(len(runsA)) for j in range(len(runsB))][:10]
    cr = np.mean([edist(pk(runsA[i]), pk(runsB[j])) for i, j in prs])
    fA = np.mean([edist(pk(runsA[i]), pk(runsA[j])) for i, j in itertools.combinations(range(len(runsA)), 2)][:10])
    fB = np.mean([edist(pk(runsB[i]), pk(runsB[j])) for i, j in itertools.combinations(range(len(runsB)), 2)][:10])
    mA = [r.mean(0) for r in runsA]; mB = [r.mean(0) for r in runsB]
    m2 = lambda a, b: float(((a - b) ** 2).sum())
    cm = np.mean([m2(mA[i], mB[j]) for i, j in prs]); fmA = np.mean([m2(mA[i], mA[j]) for i, j in itertools.combinations(range(len(mA)), 2)]); fmB = np.mean([m2(mB[i], mB[j]) for i, j in itertools.combinations(range(len(mB)), 2)])
    sp = lambda rs: np.mean([((r - r.mean(0)) ** 2).sum(1).mean() for r in rs])
    return dict(d2n=float((cr - (fA + fB) / 2) / np.sqrt(ne)), dmu=float((cm - (fmA + fmB) / 2) / ne), logspread=float(np.log(sp(runsB) / sp(runsA))))
for x in [x for x in json.load(open(ROOT + '/results/overlap_mdcath_atlas.json'))['list_cover90'] if x.get('atlas_local')]:
    c, d = x['atlas'], x['mdcath']
    if os.path.exists(f'{OUT}/{c}.json'): continue
    try:
        t0 = time.time(); rng = np.random.default_rng(37)
        h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode())
        z = zipfile.ZipFile(f'{ROOT}/raw/atlas/{c}_protein.zip'); names = z.namelist(); pdbn = [n for n in names if n.endswith('.pdb')][0]; z.extract(pdbn, TMP); top = md.load(f'{TMP}/{pdbn}').topology; sat, iat = at_atoms(top)
        a0, b0, L = block(smd, sat)
        if L < 0.9 * min(len(smd), len(sat)): raise RuntimeError(f'짧은 일치 {L}/{len(smd)},{len(sat)}')
        imd, iat = imd[a0:a0 + L], iat[b0:b0 + L]
        MD = {T: [g[T][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g[T].keys(), key=int)] for T in TT if T in g}
        h.close()
        AT = []
        for xt in sorted(n for n in names if n.endswith('.xtc')):
            z.extract(xt, TMP); t = md.load(f'{TMP}/{xt}', top=top, atom_indices=iat.ravel(), stride=20); os.remove(f'{TMP}/{xt}')
            sub_ = t.xyz.reshape(t.n_frames, L, 3, 3) * 10; AT.append(sub_)
        ca0 = MD['320'][0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
        Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
        foldMD = {T: [X[Qf(X) >= QT] for X in MD[T]] for T in MD}; foldAT = [X[Qf(X) >= QT] for X in AT]
        allf = np.concatenate([X for X in foldMD['320'] if len(X)])[:600]; ref = fs.mean_point(frames(allf))
        A_ = lambda X: fs.align_to_reference(frames(X), ref).numpy().reshape(len(X), -1)
        ok = lambda rs: [A_(X) for X in rs if len(X) >= 40]
        R320 = ok(foldMD['320']); RAT = ok(foldAT); ne = L - 1
        if len(R320) < 2 or len(RAT) < 2: raise RuntimeError('접힌 실행 부족')
        res = dict(chain=c, dom=d, L=int(L), ne=int(ne), n_runs_md=len(R320), n_runs_at=len(RAT), ff=stats(R320, RAT, ne, rng), temp={}, temp_half={})
        for T in TT[1:]:
            if T not in foldMD: continue
            RT = ok(foldMD[T])
            if len(RT) >= 2: res['temp'][T] = stats(R320, RT, ne, rng)
            H320 = [X[len(X) // 2:] for X in R320 if len(X) // 2 >= 40]; HT = [X[len(X) // 2:] for X in RT if len(X) // 2 >= 40]
            if len(H320) >= 2 and len(HT) >= 2: res['temp_half'][T] = stats(H320, HT, ne, rng)
        res['sec'] = round(time.time() - t0, 1); json.dump(res, open(f'{OUT}/{c}.json', 'w')); print(time.strftime('%H:%M:%S'), c, res['ff'], res['sec'], flush=True)
    except Exception as e:
        import traceback; print('ERR', c, repr(e)[:200], traceback.format_exc()[-300:], flush=True)
    for f in glob.glob(TMP + '/*'): os.remove(f) if os.path.isfile(f) else shutil.rmtree(f, ignore_errors=True)
J = [json.load(open(q)) for q in glob.glob(OUT + '/*.json')]; rng = np.random.default_rng(3)
def mc(v, n=2000):
    v = np.asarray(v, float); b = [np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(n)]; return [round(float(np.median(v)), 4)] + [round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])]
S = dict(n=len(J), ff={k: mc([j['ff'][k] for j in J]) for k in ('d2n', 'dmu', 'logspread')})
dT = {'348': 28, '379': 59, '413': 93, '450': 130}
for lab in ('temp', 'temp_half'):
    S[lab] = {T: dict(n=sum(T in j[lab] for j in J), **{k: mc([j[lab][T][k] for j in J if T in j[lab]]) for k in ('d2n', 'dmu', 'logspread')}) for T in dT if any(T in j[lab] for j in J)}
    xs = [0] + [dT[T] for T in S[lab]]; ys = [0.0] + [S[lab][T]['dmu'][0] for T in S[lab]]
    v = S['ff']['dmu'][0]; eq = float(np.interp(v, np.maximum.accumulate(ys), xs, right=np.nan)) if v <= max(ys) else None
    S[f'{lab}_equiv_K'] = dict(ff_dmu=v, equiv_dT=eq, minus20=(eq - 20) if eq is not None else None, note='평균 이동 Δμ²/변 의 같은 도메인 온도 곡선(중앙값, 단조로 만든 것)에 선형 보간. 곡선 밖이면 None(>+130 K).')
S['written'] = time.strftime('%F %T'); json.dump(S, open(ROOT + '/results/r37_ff.json', 'w'), indent=1, ensure_ascii=False); print(json.dumps(S, indent=1, ensure_ascii=False))
open(ROOT + '/results/logs/r37_ff.done', 'w').write(time.strftime('%F %T'))
