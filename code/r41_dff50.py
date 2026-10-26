"""(round41 g4, 2026-10-21 UTC; Drive JSON round41_2026-10-21) R3′ 창 맞춘 d_ff: 양쪽 50–100 ns.
r38_dff_D.one() 과 같은 사슬·같은 짝·같은 접힌 프레임 규칙·같은 기준점(mdCATH 320 K 접힌 프레임 앞 600 의 평균 점)에서 창만 바꾼다:
  mdCATH 320 K = 프레임 번호 50–99 (저장 1 ns), ATLAS = mdtraj 시각 50,000–100,000 ps (stride 20). N = 50 (E2 의 창 규칙; 창 안 접힌 프레임 ≥ 50 인 실행만).
프레임 고르기 seed = crc32('r41_dff50|사슬'), 쌍 고르기 seed = seeded('r41_dff50', 사슬). 실행 쌍마다 프레임 쌍 2,000 (줄인 계산). 맞춘 좌표는 전수.
출력: derived/r41/dff50/<사슬>.json, results/r41_dff50.json (기존 전체 창 r38_dff_D 와 사슬별 비)
사용: python r41_dff50.py <ROOT>"""
import sys, os, json, glob, time, itertools, zipfile, shutil, zlib
ROOT = sys.argv[1]
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]; assert _src.count('N = 100; NPAIR') == 1
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src.replace('N = 100; NPAIR', 'N = 50; NPAIR'))
import mdtraj as md
assert N == 50 and upairs.__defaults__[1] == 50
OUTF = ROOT + '/derived/r41/dff50'; os.makedirs(OUTF, exist_ok=True); TMP = '/content/dff50_tmp'; os.makedirs(TMP, exist_ok=True)
OLD = {r['chain']: r for r in json.load(open(RES + '/r38_dff_D.json'))['rows']}
W_MD = (50, 100); W_AT = (50000.0, 100000.0)


def one(c, d):
    rng = np.random.default_rng(zlib.crc32(f'r41_dff50|{c}'.encode())); rng2 = seeded('r41_dff50', c)
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode())
    z = zipfile.ZipFile(f'{ROOT}/raw/atlas/{c}_protein.zip'); names = z.namelist(); pdbn = [n for n in names if n.endswith('.pdb')][0]; z.extract(pdbn, TMP); top = md.load(f'{TMP}/{pdbn}').topology; sat, iat = at_atoms(top)
    a0, b0, L = block(smd, sat)
    if L < 0.9 * min(len(smd), len(sat)): raise RuntimeError(f'짧은 일치 {L}/{len(smd)},{len(sat)}')
    imd, iat = imd[a0:a0 + L], iat[b0:b0 + L]
    MD = [g['320'][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g['320'].keys(), key=int)]; h.close()
    AT, ATt = [], []
    for xt in sorted(n for n in names if n.endswith('.xtc')):
        z.extract(xt, TMP); t = md.load(f'{TMP}/{xt}', top=top, atom_indices=iat.ravel(), stride=20); os.remove(f'{TMP}/{xt}'); AT.append(t.xyz.reshape(t.n_frames, L, 3, 3) * 10); ATt.append(np.asarray(t.time, float))
    ca0 = MD[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
    fMDall = [X[Qf(X) >= QT] for X in MD]; allf = np.concatenate([X for X in fMDall if len(X)])[:600]; ref = torch.as_tensor(fs.mean_point(frames(allf)))   # 기준점은 r38_dff_D 와 같음
    MDw = [X[W_MD[0]:W_MD[1]] for X in MD]; ATw = [X[(tt >= W_AT[0]) & (tt <= W_AT[1])] for X, tt in zip(AT, ATt)]
    fMD = [X[Qf(X) >= QT] for X in MDw]; fAT = [X[Qf(X) >= QT] for X in ATw]
    pick = lambda rs: [X[np.sort(rng.choice(len(X), N, replace=False))] for X in rs if len(X) >= N]
    P = pick(fMD); Q = pick(fAT); ne = L - 1
    res = dict(chain=c, dom=d, L=int(L), ne=int(ne), n_runs_md=len(P), n_runs_at=len(Q), n_win_md=[len(X) for X in MDw], n_fold_md=[len(X) for X in fMD],
               n_win_at=[len(X) for X in ATw], n_fold_at=[len(X) for X in fAT], at_time_range=[[float(tt.min()), float(tt.max())] for tt in ATt])
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
                     dmu_A=(m(cr, 'a2') - (m(fp, 'a2') + m(fq, 'a2')) / 2) / ne, d2n_A=As['d2n'], ratio_A=As['ratio'])
    return res


if __name__ == '__main__':
    for x in [x for x in json.load(open(ROOT + '/results/overlap_mdcath_atlas.json'))['list_cover90'] if x.get('atlas_local')]:
        c, d = x['atlas'], x['mdcath']; fo = f'{OUTF}/{c}.json'
        if os.path.exists(fo): continue
        try:
            t0 = time.time(); r = one(c, d); r['sec'] = round(time.time() - t0, 1); jdump(r, fo); log('dff50', c, r['L'], r['n_runs_md'], r['n_runs_at'], r.get('ff', {}).get('dmu'), r['sec'], 's')
        except Exception as e:
            traceback.print_exc(); log('ERR dff50', c, repr(e)[:200])
        for f in glob.glob(TMP + '/*'): os.remove(f) if os.path.isfile(f) else shutil.rmtree(f, ignore_errors=True)
    J = [json.load(open(q)) for q in sorted(glob.glob(OUTF + '/*.json'))]
    rows = []
    for j in J:
        if 'ff' not in j: continue
        o = OLD.get(j['chain'], {}); r = dict(chain=j['chain'], dom=j['dom'], L=j['L'], n_runs_md=j['n_runs_md'], n_runs_at=j['n_runs_at'], **{k: j['ff'][k] for k in ('dmu', 'd2n', 'ratio', 'wbar', 'dmu_A', 'd2n_A', 'ratio_A')},
                                              dmu_full=o.get('dmu'), d2n_full=o.get('d2n'), ratio_full=o.get('ratio'))
        r['dmu_w50_over_full'] = r['dmu'] / r['dmu_full'] if r['dmu_full'] else None; rows.append(r)
    rb = np.random.default_rng(20261086)
    def med_ci(v):
        v = np.array([x for x in v if x is not None and np.isfinite(x)], float)
        if not len(v): return None
        bs = [np.median(v[rb.integers(0, len(v), len(v))]) for _ in range(2000)]; return dict(n=len(v), median=float(np.median(v)), ci95=[float(x) for x in np.percentile(bs, [2.5, 97.5])])
    S = dict(written=time.strftime('%F %T'), n=len(J), n_ff=len(rows), dropped=[dict(chain=j['chain'], n_runs_md=j['n_runs_md'], n_runs_at=j['n_runs_at'], n_fold_md=j.get('n_fold_md')) for j in J if 'ff' not in j], rows=rows,
             w50=dict(dmu=med_ci([r['dmu'] for r in rows]), d2n=med_ci([r['d2n'] for r in rows]), ratio=med_ci([r['ratio'] for r in rows]), dmu_A=med_ci([r['dmu_A'] for r in rows])),
             full_same_chains=dict(dmu=med_ci([r['dmu_full'] for r in rows]), d2n=med_ci([r['d2n_full'] for r in rows]), ratio=med_ci([r['ratio_full'] for r in rows])),
             ratio_w50_over_full=med_ci([r['dmu_w50_over_full'] for r in rows]),
             rule='round41 g4: 양쪽 50–100 ns (mdCATH 320 K 프레임 50–99, ATLAS 50–100 ns), N = 50, Δμ²₂_D/변 등; 사슬 중앙값, 사슬 부트스트랩 2,000 (seed 20261086). K 환산 안 함. 서술.')
    jdump(S, RES + '/r41_dff50.json', bk=True); log('DFF50_DONE', len(rows), S['w50']['dmu'], '전체 창 같은 사슬', S['full_same_chains']['dmu'], '비', S['ratio_w50_over_full'])
