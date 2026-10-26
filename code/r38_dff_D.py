"""(prereg round36 7g (b), 2026-10-20 UTC) R3′ 역장 차 d_ff 의 D 판: Δμ²₂_D = (C₂ − (B₂_P + B₂_Q)/2)/변 (지수 2 에너지 통계, r38_p2_D 와 같은 양). AlphaFlow 2차 3 의 D 판 재료.
자료·창·접힌 프레임·실행 표본은 r38_dff.py 와 같다(같은 seed crc32('r38_dff|사슬'), 같은 호출 순서 → 같은 프레임). P = mdCATH 320 K 실행, Q = ATLAS 실행(stride 20), 실행마다 N = 100.
D 교차·바닥 = 실행 쌍마다 프레임 쌍 2,000 (r38_s0 의 줄인 계산; 쌍 고르기 seed = r38_s0.seeded('r38_dff_D', 사슬)).
맞춘 좌표는 전수(실행 쌍마다 모든 프레임 쌍): 쌍선형 항등식으로 r38_dff 의 dmu 와 같아야 한다(바름 점검; 표 identity_max_rel).
출력: derived/r38/dff_D/<사슬>.json, results/r38_dff_D.json (rows: chain, dom, L, dmu = Δμ²₂_D/변, d2n, D2, C, B_P, B_Q, wbar, ratio, dmu_A, d2n_A, dmu_r38dff; 중앙 dmu 와 사슬 부트스트랩 2,000 구간)
사용: python r38_dff_D.py <ROOT>"""
import sys, os, json, glob, time, itertools, zipfile, shutil, zlib
ROOT = sys.argv[1]
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)
import mdtraj as md
OUTF = ROOT + '/derived/r38/dff_D'; os.makedirs(OUTF, exist_ok=True); TMP = '/content/dffD_tmp'; os.makedirs(TMP, exist_ok=True)
OLD = {r['chain']: r for r in json.load(open(RES + '/r38_dff.json'))['rows']}


def one(c, d):
    rng = np.random.default_rng(zlib.crc32(f'r38_dff|{c}'.encode())); rng2 = seeded('r38_dff_D', c)
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode())
    z = zipfile.ZipFile(f'{ROOT}/raw/atlas/{c}_protein.zip'); names = z.namelist(); pdbn = [n for n in names if n.endswith('.pdb')][0]; z.extract(pdbn, TMP); top = md.load(f'{TMP}/{pdbn}').topology; sat, iat = at_atoms(top)
    a0, b0, L = block(smd, sat)
    if L < 0.9 * min(len(smd), len(sat)): raise RuntimeError(f'짧은 일치 {L}/{len(smd)},{len(sat)}')
    imd, iat = imd[a0:a0 + L], iat[b0:b0 + L]
    MD = [g['320'][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g['320'].keys(), key=int)]; h.close()
    AT = []
    for xt in sorted(n for n in names if n.endswith('.xtc')):
        z.extract(xt, TMP); t = md.load(f'{TMP}/{xt}', top=top, atom_indices=iat.ravel(), stride=20); os.remove(f'{TMP}/{xt}'); AT.append(t.xyz.reshape(t.n_frames, L, 3, 3) * 10)
    ca0 = MD[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
    fMD = [X[Qf(X) >= QT] for X in MD]; fAT = [X[Qf(X) >= QT] for X in AT]
    allf = np.concatenate([X for X in fMD if len(X)])[:600]; ref = torch.as_tensor(fs.mean_point(frames(allf)))
    pick = lambda rs: [X[np.sort(rng.choice(len(X), N, replace=False))] for X in rs if len(X) >= N] # r38_dff 와 같은 호출 순서
    P = pick(fMD); Q = pick(fAT); ne = L - 1
    res = dict(chain=c, dom=d, L=int(L), ne=int(ne), n_runs_md=len(P), n_runs_at=len(Q))
    if len(P) < 2 or len(Q) < 2: return res
    E = {'P': coords(np.concatenate(P)), 'Q': coords(np.concatenate(Q))}
    A = {k: fs.align_to_reference(frames(np.concatenate(v)), ref).numpy().reshape(len(v) * N, -1) for k, v in (('P', P), ('Q', Q))}
    def go(a, b, rr):
        o = runpair_means(MeanStart(E[a], E[b]), rr, rng2); out = {}
        for g_, rp in enumerate(rr):
            s_ = o['GR'] == g_; full = cost(A[a][rp[0] * N:(rp[0] + 1) * N], A[b][rp[1] * N:(rp[1] + 1) * N])
            out[rp] = dict(d=float(o['dv'][s_].mean()), d2=float((o['dv'][s_] ** 2).mean()), a=float(np.sqrt(np.maximum(full, 0)).mean()), a2=float(full.mean()))
        return out
    fp = go('P', 'P', list(itertools.combinations(range(len(P)), 2))); fq = go('Q', 'Q', list(itertools.combinations(range(len(Q)), 2)))
    cr = go('P', 'Q', list(itertools.product(range(len(P)), range(len(Q)))))
    m = lambda xs, k: float(np.mean([v[k] for v in xs.values()]))
    Ds = estats(m(cr, 'd'), m(fp, 'd'), m(fq, 'd'), ne); As = estats(m(cr, 'a'), m(fp, 'a'), m(fq, 'a'), ne)
    res['ff'] = dict(dmu=(m(cr, 'd2') - (m(fp, 'd2') + m(fq, 'd2')) / 2) / ne, d2n=Ds['d2n'], D2=Ds['D2'], C=Ds['C'], B_P=Ds['B_P'], B_Q=Ds['B_Q'], wbar=Ds['wbar'], ratio=Ds['ratio'],
                     dmu_A=(m(cr, 'a2') - (m(fp, 'a2') + m(fq, 'a2')) / 2) / ne, d2n_A=As['d2n'])
    return res


if __name__ == '__main__':
    for x in [x for x in json.load(open(ROOT + '/results/overlap_mdcath_atlas.json'))['list_cover90'] if x.get('atlas_local')]:
        c, d = x['atlas'], x['mdcath']; fo = f'{OUTF}/{c}.json'
        if os.path.exists(fo): continue
        try:
            t0 = time.time(); r = one(c, d); r['sec'] = round(time.time() - t0, 1); jdump(r, fo); log('dff_D', c, r['L'], r.get('ff', {}).get('dmu'), r['sec'], 's')
        except Exception as e:
            traceback.print_exc(); log('ERR dff_D', c, repr(e)[:200])
        for f in glob.glob(TMP + '/*'): os.remove(f) if os.path.isfile(f) else shutil.rmtree(f, ignore_errors=True)
    J = [json.load(open(q)) for q in sorted(glob.glob(OUTF + '/*.json'))]
    rows = [dict(chain=j['chain'], dom=j['dom'], L=j['L'], **{k: j['ff'][k] for k in ('dmu', 'd2n', 'D2', 'C', 'B_P', 'B_Q', 'wbar', 'ratio', 'dmu_A', 'd2n_A')},
                 dmu_r38dff=OLD.get(j['chain'], {}).get('dmu')) for j in J if 'ff' in j]
    v = np.array([r['dmu'] for r in rows]); rb = np.random.default_rng(20261063); bs = [np.median(v[rb.integers(0, len(v), len(v))]) for _ in range(2000)] if len(v) else []
    idr = [abs(r['dmu_A'] - r['dmu_r38dff']) / max(abs(r['dmu_r38dff']), 1e-12) for r in rows if r['dmu_r38dff'] is not None]
    S = dict(written=time.strftime('%F %T'), n=len(J), n_ff=len(rows), rows=rows, median_dmu=float(np.median(v)) if len(v) else None,
             ci_median=[float(x) for x in np.percentile(bs, [2.5, 97.5])] if len(v) else None, median_dmu_A=float(np.median([r['dmu_A'] for r in rows])) if rows else None,
             identity_max_rel=float(max(idr)) if idr else None, n_identity=len(idr),
             rule='prereg 7g (b): Δμ²₂_D/변, 사슬 중앙값, 사슬 부트스트랩 2,000 (seed 20261063). 맞춘 좌표 전수 dmu_A 는 r38_dff 의 dmu 와 항등(같은 프레임)이어야 함.')
    jdump(S, RES + '/r38_dff_D.json', bk=True); log('DFF_D_DONE', len(rows), '중앙 dmu D', S['median_dmu'], S['ci_median'], '맞춘 좌표', S['median_dmu_A'], '항등 최대 상대 차', S['identity_max_rel'])
