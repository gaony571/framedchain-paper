"""(prereg round37 R3″) 규약만 다른 칸: ONE_cath1 대 MSR_cath2 (둘 다 amber ff99sb-ildn, 300 K, 10 ns 저장; ONE 은 시스템당 누적 > 100 µs).
겹치는 16 도메인. 실행 = 반복. 공통 잔기(서열 정렬의 가장 긴 연속 일치), 접힌 프레임(MSR 첫 실행 첫 프레임의 천연 접촉, Q ≥ 0.75).
derived/r37_protocol/<dom>.json, results/r37_protocol.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools, zipfile, shutil, difflib
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
import mdtraj as md
OUT = ROOT + '/derived/r37_protocol'; os.makedirs(OUT, exist_ok=True); TMP = '/content/pr_tmp'; os.makedirs(TMP, exist_ok=True)
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', HSD='H', HSE='H', HIE='H', HID='H', HIP='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')
def bb(top):
    keep = [r for r in top.residues if r.is_protein and all(any(a.name == n for a in r.atoms) for n in ('CA', 'C', 'O'))]
    return ''.join(AA.get(r.name, 'X') for r in keep), np.array([[next(a.index for a in r.atoms if a.name == n) for n in ('CA', 'C', 'O')] for r in keep])
def stats(runsA, runsB, ne, rng, N=100):
    pk = lambda X: X[rng.choice(len(X), min(N, len(X)), replace=False)]
    prs = [(i, j) for i in range(len(runsA)) for j in range(len(runsB))]; rng.shuffle(prs); prs = prs[:10]
    ca = list(itertools.combinations(range(len(runsA)), 2)); cb = list(itertools.combinations(range(len(runsB)), 2)); rng.shuffle(ca); rng.shuffle(cb)
    cr = np.mean([edist(pk(runsA[i]), pk(runsB[j])) for i, j in prs]); fA = np.mean([edist(pk(runsA[i]), pk(runsA[j])) for i, j in ca[:10]]); fB = np.mean([edist(pk(runsB[i]), pk(runsB[j])) for i, j in cb[:10]])
    mA = [r.mean(0) for r in runsA]; mB = [r.mean(0) for r in runsB]; m2 = lambda a, b: float(((a - b) ** 2).sum())
    cm = np.mean([m2(mA[i], mB[j]) for i, j in prs]); fmA = np.mean([m2(mA[i], mA[j]) for i, j in ca[:10]]); fmB = np.mean([m2(mB[i], mB[j]) for i, j in cb[:10]])
    sp = lambda rs: np.mean([((r - r.mean(0)) ** 2).sum(1).mean() for r in rs])
    return dict(d2n=float((cr - (fA + fB) / 2) / np.sqrt(ne)), dmu=float((cm - (fmA + fmB) / 2) / ne), logspread=float(np.log(sp(runsB) / sp(runsA))))
ZB = {k: zipfile.ZipFile(f'{ROOT}/raw/bioemu_cath/{k}.zip') for k in ('MSR_cath2', 'ONE_cath1')}
NL = {k: z.namelist() for k, z in ZB.items()}
one = {n.split('/')[1].replace('cath1_', '') for n in NL['ONE_cath1'] if n.count('/') >= 2}; msr = {n.split('/')[1].replace('cath2_', '') for n in NL['MSR_cath2'] if n.count('/') >= 2}
D16 = sorted(one & msr)
def runs(key, pref, d, top_fallback=None):
    pre = f'{key}/{pref}{d}/'; tn = pre + 'topology.pdb'
    if tn in NL[key]: ZB[key].extract(tn, TMP); top = md.load(f'{TMP}/{tn}').topology
    else: top = top_fallback
    out = []
    for n in sorted(x for x in NL[key] if x.startswith(pre + 'trajs/') and x.endswith('.xtc')):
        ZB[key].extract(n, TMP); t = md.load(f'{TMP}/{n}', top=top); os.remove(f'{TMP}/{n}'); out.append(t)
    return top, out
for d in D16:
    if os.path.exists(f'{OUT}/{d}.json'): continue
    try:
        t0 = time.time(); rng = np.random.default_rng(37)
        topM, TM_ = runs('MSR_cath2', 'cath2_', d); topO, TO = runs('ONE_cath1', 'cath1_', d, top_fallback=topM)
        sM, iM = bb(topM); sO, iO = bb(topO)
        m = max(difflib.SequenceMatcher(None, sM, sO, autojunk=False).get_matching_blocks(), key=lambda b: b.size); L = m.size
        if L < 0.9 * min(len(sM), len(sO)): raise RuntimeError(f'짧은 일치 {L}')
        iM, iO = iM[m.a:m.a + L], iO[m.b:m.b + L]
        XM = [t.xyz[:, iM.ravel()].reshape(t.n_frames, L, 3, 3) * 10 for t in TM_]; XO = [t.xyz[:, iO.ravel()].reshape(t.n_frames, L, 3, 3) * 10 for t in TO]
        ca0 = XM[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); k_ = dn < 8.0; ii, jj, dn = ii[k_], jj[k_], dn[k_]
        Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
        fM = [X[Qf(X) >= QT] for X in XM]; fO = [X[Qf(X) >= QT] for X in XO]
        ref = fs.mean_point(frames(np.concatenate([X for X in fM if len(X)])[:600]))
        A_ = lambda X: fs.align_to_reference(frames(X), ref).numpy().reshape(len(X), -1)
        RM = [A_(X) for X in fM if len(X) >= 20]; RO = [A_(X) for X in fO if len(X) >= 20]; ne = L - 1
        if len(RM) < 2 or len(RO) < 2: raise RuntimeError(f'실행 부족 MSR {len(RM)} ONE {len(RO)}')
        res = dict(dom=d, L=int(L), ne=int(ne), n_runs_msr=len(RM), n_runs_one=len(RO), frames_msr=int(sum(len(x) for x in RM)), frames_one=int(sum(len(x) for x in RO)),
                   protocol=stats(RM, RO, ne, rng), sec=round(time.time() - t0, 1))
        json.dump(res, open(f'{OUT}/{d}.json', 'w')); print(time.strftime('%H:%M:%S'), d, res['protocol'], res['n_runs_msr'], res['n_runs_one'], flush=True)
    except Exception as e:
        import traceback; print('ERR', d, repr(e)[:200], traceback.format_exc()[-300:], flush=True)
    for f in glob.glob(TMP + '/*'): os.remove(f) if os.path.isfile(f) else shutil.rmtree(f, ignore_errors=True)
J = [json.load(open(q)) for q in glob.glob(OUT + '/*.json')]; rng = np.random.default_rng(3)
def mc(v, n=2000):
    v = np.asarray(v, float); b = [np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(n)]; return [round(float(np.median(v)), 4)] + [round(float(x), 4) for x in np.percentile(b, [2.5, 97.5])]
S = dict(n=len(J), protocol={k: mc([j['protocol'][k] for j in J]) for k in ('d2n', 'dmu', 'logspread')} if J else None, compare_R3_ff_plus_protocol=dict(d2n=0.076, dmu=0.0306, logspread=0.139, note='R3: SM 대 BioEmu MSR, 다른 54 도메인'))
S['written'] = time.strftime('%F %T'); json.dump(S, open(ROOT + '/results/r37_protocol.json', 'w'), indent=1, ensure_ascii=False); print(json.dumps(S, indent=1, ensure_ascii=False))
open(ROOT + '/results/logs/r37_protocol.done', 'w').write(time.strftime('%F %T'))
