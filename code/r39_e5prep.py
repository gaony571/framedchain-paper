"""(round39 E5 방향 검정의 재료, 맞춘 좌표) R3′·R3″ 칸의 실행 평균(9 변 성분)과 변별 퍼짐을 저장한다. 7c 칸은 derived/r38/compare/<dom>.npz 에 이미 있다.
R3′ (역장 + 규약 + ΔT 20 K; 44 사슬 = r38_dff 와 같은 목록): r38_dff.py 와 같은 읽기·접힌 프레임·기준 평균 점·rng 순서로 MD 320 K 실행(P)과 ATLAS 실행(X)을 다시 뽑는다
  (같은 프레임이라 compare(P, X)['dmu'] 가 results/r38_dff.json 과 같아야 함 — 저장하는 dmu_check). 같은 rng 를 이어 mdCATH 348–450 K 실행(전체 창, Q ≥ 0.75, 반복마다 N = 100)을 더 뽑는다.
R3″ (규약만; r37_protocol 의 16 도메인): MSR_cath2 실행(X 의 기준 Xref)과 ONE_cath1 실행(X) — r37_protocol 과 같은 방향(ONE 대 MSR). 잔기는 MSR∩ONE 가장 긴 일치 블록과 mdCATH 서열의 가장 긴 일치 블록.
  접힌 프레임 기준은 이 판에서 mdCATH 320 K 첫 프레임의 천연 접촉으로 통일한다(r37_protocol 은 MSR 첫 프레임). 실행마다 N = 100 무작위(접힌 프레임이 100 보다 적고 20 이상이면 모두, 수 기록).
  mdCATH 320·348–450 K 실행은 R3′ 와 같은 규칙. 모두 mdCATH 320 K 접힌 프레임 앞 600 의 평균 점에 맞춘다.
사용: python r39_e5prep.py <ROOT> r3p|r3pp|r3p,r3pp 출력: derived/r39/e5/<사슬 또는 도메인>.npz (키 X|mean, X|spr, Xref|mean, Xref|spr, MD<T>|mean, MD<T>|spr; 메타는 .json)"""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, glob, time, itertools, zipfile, shutil, difflib, zlib, traceback
ROOT = sys.argv[1]; PARTS = sys.argv[2].split(','); sys.argv = [sys.argv[0], ROOT, '__none__']
exec(open(ROOT + '/code/c1_core.py').read().split('doms = ONLY')[0])
src = open(ROOT + '/code/r37_ff.py').read(); exec(src[src.index('AA = dict('):src.index('def stats')])
cs = open(ROOT + '/code/r38_compare.py').read(); exec(cs[cs.index('def pdist_mean'):cs.index('\nfor d in doms:')])
import h5py, mdtraj as md
OUT = ROOT + '/derived/r39/e5'; os.makedirs(OUT, exist_ok=True); TMP = '/content/e5_tmp'; os.makedirs(TMP, exist_ok=True); N = 100; MINF = 20
assert QT == 0.75


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def save(tag, prof, meta):
    np.savez_compressed(f'{OUT}/{tag}.npz', **{k: np.asarray(v) for k, v in prof.items()}); json.dump(meta, open(f'{OUT}/{tag}.json', 'w'), ensure_ascii=False)


def md_T(g, imd, L, Qf, A_, rng, prof, meta):
    """mdCATH 348–450 K: 전체 창, 접힌 프레임 ≥ N 인 반복만, 반복마다 N 무작위."""
    for T in TT[1:]:
        if T not in g: continue
        runs = []
        for r in sorted(g[T].keys(), key=int):
            X = g[T][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3); X = X[Qf(X) >= QT]
            if len(X) >= N: runs.append(A_(X[np.sort(rng.choice(len(X), N, replace=False))]))
        meta['nruns'][f'MD{T}'] = len(runs)
        if len(runs) >= 2: prof[f'MD{T}|mean'], prof[f'MD{T}|spr'] = edge_profile(runs, L - 1)


def part_r3p():
    dff = {r['chain']: r for r in json.load(open(ROOT + '/results/r38_dff.json'))['rows']}
    for x in [x for x in json.load(open(ROOT + '/results/overlap_mdcath_atlas.json'))['list_cover90'] if x.get('atlas_local')]:
        c, d = x['atlas'], x['mdcath']
        if os.path.exists(f'{OUT}/{c}.npz') or c not in dff: continue
        try:
            t0 = time.time(); rng = np.random.default_rng(zlib.crc32(f'r38_dff|{c}'.encode()))
            h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode())
            z = zipfile.ZipFile(f'{ROOT}/raw/atlas/{c}_protein.zip'); names = z.namelist(); pdbn = [n for n in names if n.endswith('.pdb')][0]; z.extract(pdbn, TMP); top = md.load(f'{TMP}/{pdbn}').topology; sat, iat = at_atoms(top)
            a0, b0, L = block(smd, sat)
            if L < 0.9 * min(len(smd), len(sat)): raise RuntimeError(f'짧은 일치 {L}/{len(smd)},{len(sat)}')
            imd, iat = imd[a0:a0 + L], iat[b0:b0 + L]
            MD = [g['320'][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g['320'].keys(), key=int)]
            AT = []
            for xt in sorted(n for n in names if n.endswith('.xtc')):
                z.extract(xt, TMP); t = md.load(f'{TMP}/{xt}', top=top, atom_indices=iat.ravel(), stride=20); os.remove(f'{TMP}/{xt}'); AT.append(t.xyz.reshape(t.n_frames, L, 3, 3) * 10)
            ca0 = MD[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
            Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
            fMD = [X[Qf(X) >= QT] for X in MD]; fAT = [X[Qf(X) >= QT] for X in AT]
            allf = np.concatenate([X for X in fMD if len(X)])[:600]; ref = fs.mean_point(frames(allf))
            A_ = lambda X: fs.align_to_reference(frames(X), ref).numpy().reshape(len(X), -1)
            pick = lambda rs: [A_(X[np.sort(rng.choice(len(X), N, replace=False))]) for X in rs if len(X) >= N]
            P = pick(fMD); Q = pick(fAT); ne = L - 1
            meta = dict(cell='R3p', tag=c, dom=d, L=int(L), ne=int(ne), nruns=dict(MD320=len(P), X=len(Q)), frames_rule='r38_dff 와 같은 프레임(rng 순서 같음)')
            if len(P) < 2 or len(Q) < 2: raise RuntimeError(f'실행 부족 {len(P)} {len(Q)}')
            ff = compare(P, Q, ne); meta['dmu_check'] = dict(now=ff['dmu'], r38_dff=dff[c]['dmu'], rel=abs(ff['dmu'] - dff[c]['dmu']) / max(abs(dff[c]['dmu']), 1e-12))
            prof = {}; prof['MD320|mean'], prof['MD320|spr'] = edge_profile(P, ne); prof['X|mean'], prof['X|spr'] = edge_profile(Q, ne)
            prof['Xref|mean'], prof['Xref|spr'] = prof['MD320|mean'], prof['MD320|spr']; meta['xref_is_md320'] = True
            md_T(g, imd, L, Qf, A_, rng, prof, meta); h.close(); meta['sec'] = round(time.time() - t0, 1); save(c, prof, meta)
            log('e5 R3′', c, d, L, meta['nruns'], 'dmu 재현 상대 차', f"{meta['dmu_check']['rel']:.2e}", meta['sec'], 's')
        except Exception as e:
            traceback.print_exc(); log('ERR e5 R3′', c, repr(e)[:200])
        for f in glob.glob(TMP + '/*'): os.remove(f) if os.path.isfile(f) else shutil.rmtree(f, ignore_errors=True)
    log('E5PREP_R3P_DONE')


def part_r3pp():
    ps = open(ROOT + '/code/r37_protocol.py').read(); exec(ps[ps.index('def bb(top):'):ps.index('def stats(')], globals())
    ZB = {k: zipfile.ZipFile(f'{ROOT}/raw/bioemu_cath/{k}.zip') for k in ('MSR_cath2', 'ONE_cath1')}; NL = {k: z.namelist() for k, z in ZB.items()}
    def runs(key, pref, d, top_fallback=None):
        pre = f'{key}/{pref}{d}/'; tn = pre + 'topology.pdb'
        if tn in NL[key]: ZB[key].extract(tn, TMP); top = md.load(f'{TMP}/{tn}').topology
        else: top = top_fallback
        out = []
        for n in sorted(x for x in NL[key] if x.startswith(pre + 'trajs/') and x.endswith('.xtc')):
            ZB[key].extract(n, TMP); t = md.load(f'{TMP}/{n}', top=top); os.remove(f'{TMP}/{n}'); out.append(t)
        return top, out
    for d in sorted(os.path.basename(f)[:-5] for f in glob.glob(ROOT + '/derived/r37_protocol/*.json')):
        if os.path.exists(f'{OUT}/{d}.npz'): continue
        try:
            t0 = time.time(); rng = np.random.default_rng(zlib.crc32(f'r39_e5|{d}'.encode()))
            topM, TM_ = runs('MSR_cath2', 'cath2_', d); topO, TO = runs('ONE_cath1', 'cath1_', d, top_fallback=topM)
            sM, iM = bb(topM); sO, iO = bb(topO)
            m = max(difflib.SequenceMatcher(None, sM, sO, autojunk=False).get_matching_blocks(), key=lambda b: b.size); L1 = m.size
            if L1 < 0.9 * min(len(sM), len(sO)): raise RuntimeError(f'짧은 일치 MSR·ONE {L1}')
            sC, iM, iO = sM[m.a:m.a + L1], iM[m.a:m.a + L1], iO[m.b:m.b + L1]
            fh = f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5'
            if not os.path.exists(fh): raise RuntimeError('mdCATH h5 없음')
            h = h5py.File(fh, 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode())
            a0, b0, L = block(sC, smd)
            if L < 0.9 * min(len(sC), len(smd)): raise RuntimeError(f'짧은 일치 mdCATH {L}/{len(sC)},{len(smd)}')
            iM, iO, imd = iM[a0:a0 + L], iO[a0:a0 + L], imd[b0:b0 + L]
            MD = [g['320'][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g['320'].keys(), key=int)]
            XM = [t.xyz[:, iM.ravel()].reshape(t.n_frames, L, 3, 3) * 10 for t in TM_]; XO = [t.xyz[:, iO.ravel()].reshape(t.n_frames, L, 3, 3) * 10 for t in TO]
            ca0 = MD[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
            Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
            fMD = [X[Qf(X) >= QT] for X in MD]; allf = np.concatenate([X for X in fMD if len(X)])[:600]; ref = fs.mean_point(frames(allf))
            A_ = lambda X: fs.align_to_reference(frames(X), ref).numpy().reshape(len(X), -1)
            pickN = lambda rs: [A_(X[np.sort(rng.choice(len(X), N, replace=False))]) for X in rs if len(X) >= N]
            def pick_bio(rs):
                out, nf = [], []
                for X in rs:
                    X = X[Qf(X) >= QT]
                    if len(X) >= N: out.append(A_(X[np.sort(rng.choice(len(X), N, replace=False))])); nf.append(N)
                    elif len(X) >= MINF: out.append(A_(X)); nf.append(len(X))
                return out, nf
            P = pickN(fMD); RM, nM = pick_bio(XM); RO, nO = pick_bio(XO); ne = L - 1
            meta = dict(cell='R3pp', tag=d, dom=d, L=int(L), L_msr_one=int(L1), ne=int(ne), nruns=dict(MD320=len(P), Xref=len(RM), X=len(RO)), frames_msr=nM, frames_one=nO,
                        direction='X = ONE_cath1, Xref = MSR_cath2 (r37_protocol 과 같은 방향; 규약 차의 부호는 물리적 방향이 없음)', xref_is_md320=False)
            if len(P) < 2 or len(RM) < 2 or len(RO) < 2: raise RuntimeError(f'실행 부족 MD320 {len(P)} MSR {len(RM)} ONE {len(RO)}')
            prof = {}; prof['MD320|mean'], prof['MD320|spr'] = edge_profile(P, ne); prof['X|mean'], prof['X|spr'] = edge_profile(RO, ne); prof['Xref|mean'], prof['Xref|spr'] = edge_profile(RM, ne)
            md_T(g, imd, L, Qf, A_, rng, prof, meta); h.close(); meta['sec'] = round(time.time() - t0, 1); save(d, prof, meta)
            log('e5 R3″', d, L, meta['nruns'], meta['sec'], 's')
        except Exception as e:
            traceback.print_exc(); log('ERR e5 R3″', d, repr(e)[:200])
        for f in glob.glob(TMP + '/*'): os.remove(f) if os.path.isfile(f) else shutil.rmtree(f, ignore_errors=True)
    log('E5PREP_R3PP_DONE')


if __name__ == '__main__':
    for part in PARTS:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('E5PREP_DONE')
