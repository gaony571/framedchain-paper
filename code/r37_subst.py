"""(prereg round37 R5′) 치환 가지를 줄기에 붙이는 칸. MegaSim 야생형 ∩ mdCATH 중 서열 덮음 ≥ 0.875 인 5 계(1EM7, 1GJS, 1LP1, 1TG0, 2MI6).
같은 도메인의 mdCATH 온도 곡선(320 대 348/379/413/450 K, 전체·뒤 절반; Δ²/√변)을 내고, 그 단백질의 단일 치환 Δ²/√변(C4, MegaSim 좌표; derived/c4_sub)의
중앙값을 그 곡선에 보간해 온도로 환산한다. 역장(ff99SB-disp 295 K 대 CHARMM22* 320 K)이 섞인다. results/r37_subst.json."""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools
ROOT = sys.argv[1]; sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
import h5py
src = open(ROOT + '/code/r37_ff.py').read()
exec(src[src.index('AA = dict('):src.index('for x in [x for x in json.load')])  # AA, TT, md_atoms, at_atoms, block, stats
PAIRS = {'1EM7': '1em7A00', '1GJS': '1gjsA00', '1LP1': '1lp1A00', '1TG0': '1tg0A00', '2MI6': '2mi6A00'}
dT = {'348': 28, '379': 59, '413': 93, '450': 130}; OUTJ = {}
for wt, d in PAIRS.items():
    try:
        rng = np.random.default_rng(37)
        with h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r') as h:
            g = h[d]; s_, idx = md_atoms(g['pdbProteinAtoms'][()].decode()); L = len(s_)
            MD = {T: [g[T][r]['coords'][:][:, idx.ravel()].reshape(-1, L, 3, 3) for r in sorted(g[T].keys(), key=int)] for T in TT if T in g}
        ca0 = MD['320'][0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); k_ = dn < 8.0; ii, jj, dn = ii[k_], jj[k_], dn[k_]
        Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
        fold = {T: [X[Qf(X) >= QT] for X in MD[T]] for T in MD}
        ref = fs.mean_point(frames(np.concatenate([X for X in fold['320'] if len(X)])[:600]))
        A_ = lambda X: fs.align_to_reference(frames(X), ref).numpy().reshape(len(X), -1)
        R = {T: [A_(X) for X in fold[T] if len(X) >= 40] for T in fold}; ne = L - 1
        curve, curve_h = {}, {}
        for T in dT:
            if len(R.get(T, [])) >= 2: curve[T] = stats(R['320'], R[T], ne, rng)['d2n']
            H0 = [X[len(X) // 2:] for X in R['320'] if len(X) // 2 >= 40]; HT = [X[len(X) // 2:] for X in R.get(T, []) if len(X) // 2 >= 40]
            if len(H0) >= 2 and len(HT) >= 2: curve_h[T] = stats(H0, HT, ne, rng)['d2n']
        sub_vals = []
        for p in glob.glob(f'{ROOT}/derived/c4_sub/*.json'):
            j = json.load(open(p))
            if str(j.get('wt', '')).upper().startswith(wt) and 'mutants' in j and 'floor' in j:
                ne_s = j['nres'] - 1; sub_vals += [(m['ed'] - j['floor']) / np.sqrt(ne_s) for m in j['mutants'] if m.get('ed') is not None and m.get('single', True)]
        sv = float(np.median(sub_vals)) if sub_vals else None
        def eq(c):
            if sv is None or not c: return None
            xs = [0] + [dT[T] for T in c]; ys = np.maximum.accumulate([0.0] + [c[T] for T in c]); return float(np.interp(sv, ys, xs)) if sv <= ys[-1] else None
        OUTJ[wt] = dict(dom=d, ne=int(ne), curve=curve, curve_half=curve_h, n_single=len(sub_vals), subst_d2n_median=sv, equiv_dT_full=eq(curve), equiv_dT_half=eq(curve_h))
        print(wt, OUTJ[wt], flush=True)
    except Exception as e:
        import traceback; print('ERR', wt, repr(e)[:200], traceback.format_exc()[-300:], flush=True)
OUTJ['note'] = '치환 Δ² 는 MegaSim(ff99SB-disp 295 K) 좌표, 온도 곡선은 mdCATH(CHARMM22*) — 역장이 섞인 같은 단백질 비교. 사례 5 계.'
OUTJ['written'] = time.strftime('%F %T'); json.dump(OUTJ, open(ROOT + '/results/r37_subst.json', 'w'), indent=1, ensure_ascii=False)
open(ROOT + '/results/logs/r37_subst.done', 'w').write(time.strftime('%F %T'))
