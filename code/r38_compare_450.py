"""(prereg 7c 2차 4 용 보충) 67 비교(r38_compare.py)가 aSAMt 450 K 생성(5 단계) 전에 돌아 빠진 aSAMt450 칸만 더한다. 다른 칸은 바꾸지 않는다.
MD 실행은 r38_compare 와 같은 seed·같은 호출 순서로 다시 뽑아 같은 프레임이 되게 한다(점검: 다시 낸 MD450-MD320|full|0.75 가 저장값과 같아야 함).
모델 반쪽 나눔 seed 만 따로 둔다: crc32('r38_compare_450|도메인'). 칸: aSAMt450-MD320|창|q, aSAMt450-MD450|창|q (창 full·half, q 0.75·0.7·0.8·all).
기존 JSON 은 derived/r38/compare/_pre450/ 에 사본을 남기고 칸을 더해 다시 쓴다.
사용: python r38_compare_450.py <ROOT> <도메인 목록 파일>"""
import sys, os, json, shutil
ROOT_ = sys.argv[1]; LIST_ = sys.argv[2]
_c = open(ROOT_ + '/code/r38_compare.py').read(); _c = _c[:_c.index('\nfor d in doms:')]
sys.argv = ['r38_compare', ROOT_, LIST_]; exec(_c)
BK = OUT + '/_pre450'; os.makedirs(BK, exist_ok=True)
for d in doms:
    fo = f'{OUT}/{d}.json'; fx = f'{ROOT}/derived/r38/asamt/{d}.450'
    if not os.path.exists(fo) or not os.path.exists(fx + '.traj.dcd'): print('건너뜀', d, '비교 또는 450 K 표본 없음', flush=True); continue
    old = json.load(open(fo))
    if any(k.startswith('aSAMt450-') for k in old['cells']): print('이미 있음', d, flush=True); continue
    t0 = time.time(); rng = np.random.default_rng(zlib.crc32(f'r38_compare|{d}'.encode()))
    try:
        z = np.load(f'{ROOT}/derived/r38/md/{d}.npz'); ref = torch.as_tensor(z['ref']); ii, jj, dn = z['ii'], z['jj'], z['dn']
        h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode()); L = len(smd)
        MD = {T: [g[T][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g[T].keys(), key=int)] for T in TT if T in g}; h.close()
        assert np.array_equal(imd, z['imd']), 'imd 다름'
        ne = L - 1; Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
        A_ = lambda X: fs.align_to_reference(frames(X), ref).numpy().reshape(len(X), -1)
        mdrun = {}
        for T in MD:
            for W_ in ('full', 'half'):
                reps = [X if W_ == 'full' else X[len(X) // 2:] for X in MD[T]]
                Qv = [Qf(X) for X in reps]; Ar = [A_(X) for X in reps]
                for q in QS + ('all',):
                    runs = []
                    for X, qq in zip(Ar, Qv):
                        idx = np.arange(len(X)) if q == 'all' else np.where(qq >= q)[0]
                        if len(idx) >= N: runs.append(X[np.sort(rng.choice(idx, N, replace=False))])
                    mdrun[(T, W_, q)] = runs
        chk = None
        if len(mdrun.get(('320', 'full', 0.75), [])) >= 2 and len(mdrun.get(('450', 'full', 0.75), [])) >= 2 and 'MD450-MD320|full|0.75' in old['cells']:
            chk = abs(compare(mdrun[('320', 'full', 0.75)], mdrun[('450', 'full', 0.75)], ne)['D2'] - old['cells']['MD450-MD320|full|0.75']['D2'])
        rng2 = np.random.default_rng(zlib.crc32(f'r38_compare_450|{d}'.encode()))
        t = md.load(fx + '.traj.dcd', top=fx + '.top.pdb'); sm, im = at_atoms(t.topology); a0, b0, Lb = block(smd, sm)
        old['check']['aSAMt450'] = dict(L_model=len(sm), L_md=L, block=int(Lb), same_seq=bool(sm == smd), n_frames=int(t.n_frames))
        if Lb == L and sm == smd:
            Xm = t.xyz[:, im.ravel()].reshape(t.n_frames, L, 3, 3) * 10; qm = Qf(Xm); Amall = A_(Xm); old['check']['aSAMt450']['folded_frac_q075'] = float((qm >= 0.75).mean())
            for W_ in ('full', 'half'):
                for q in QS + ('all',):
                    Am = Amall if q == 'all' else Amall[np.where(qm >= q)[0]]
                    for key in (('320', W_, q), ('450', W_, q)):
                        P = mdrun.get(key, [])
                        if len(P) < 2 or len(Am) < 2 * MINH: continue
                        outs = []
                        for s in range(NSPLIT):
                            perm = rng2.permutation(len(Am)); hf = len(Am) // 2; outs.append(compare(P, [Am[perm[:hf]], Am[perm[hf:2 * hf]]], ne))
                        old['cells'][f'aSAMt450-MD{key[0]}|{W_}|{q}'] = avg(outs)
        else:
            old['check']['aSAMt450']['note'] = '서열 또는 길이 다름 — 칸을 내지 않음'
        old['add450'] = dict(written=time.strftime('%F %T'), md_runs_check_absD2=chk, seed='r38_compare_450', sec=round(time.time() - t0, 1))
        shutil.copy(fo, f'{BK}/{d}.json'); json.dump(old, open(fo + '.tmp', 'w')); os.replace(fo + '.tmp', fo)
        c = old['cells']; print(time.strftime('%H:%M:%S'), d, 'MD 실행 점검', chk, 'aSAMt450-MD450 d2n', round(c.get('aSAMt450-MD450|full|0.75', {}).get('d2n', np.nan), 4), 'aSAMt320-MD320 d2n', round(c.get('aSAMt320-MD320|full|0.75', {}).get('d2n', np.nan), 4), old['add450']['sec'], 's', flush=True)
    except Exception as e:
        import traceback; traceback.print_exc(); print('ERR', d, repr(e)[:200], flush=True)
print('COMPARE450_DONE', flush=True)
