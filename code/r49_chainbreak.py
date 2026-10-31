"""(등록 = Drive results/prereg_2026-10-31_addendum12.json "r49_chainbreak") 사슬 끊김(긴 Cα–Cα 결합) 점검을 2026-10-28 점검(줄기 336, H-풀 1,496, C 949)에 들지 않은 집합으로 넓힘.
문턱: Cα–Cα > 4.2 Å (s = 길이 / 3.8 Å > 1.1), 10-28 점검과 같음.
집합과 좌표(각 분석이 읽은 원자 목록 그대로):
  model  = 모델 오차 67 도메인(derived/r38/compare_D/<d>.json 의 이름).  mdCATH h5: md_atoms(pdbProteinAtoms) 의 CA, 320 K 반복마다 첫 프레임.
  b6     = 서열 간 비교 281 쌍(results/r37_trunk_pairs.tsv)의 도메인; 분석에 든 234 쌍(derived/r38/b6/*.json 가운데 skip 없고 D 있는 것)의 도메인은 따로 셈.  같은 좌표.
  r10a   = 자연 삽입·결실 주 19 쌍(results/r37_trunk_pairs.json 의 r10a_prime.main)의 도메인.  같은 좌표.
  bnb    = 인증 500 쌍(derived/r38/static/bnb.jsonl 의 ab)의 도메인.  같은 좌표.
  p1     = 첫 호열 연구 도메인(derived/r38/p1/pairs.json 의 domains).  같은 좌표.
  p1n    = 호열 재현 도메인(derived/r38/p1n/pairs.json 의 domains): r38_p1n 규칙대로 raw/mdcath_p1n/<d>.npz 가 있으면 그것(reps320 마다 h320_<r> 첫 프레임의 CA), 없으면 h5.
  atlas  = 역장 비교 37 사슬(results/r41_dff50.json 의 rows): r41_dff50 과 같은 공통 잔기(block = difflib 가장 긴 일치 구간, at_atoms 는 mdtraj 위상).
           ATLAS 쪽 = zip 안 PDB 의 좌표(공통 잔기), mdCATH 쪽 = 같은 공통 잔기의 h5 320 K 반복마다 첫 프레임.
  static = R2 의 정적 사슬(results/v_residue_unit_graded.tsv 의 A·B): hu1_rulers.bb 와 같은 좌표(derived/xprot_structs.npz 의 '<이름>|CA', 없으면 derived/atlas_pdb/<이름>.pdb 를 framespace.read_pdb_backbone).
생성 표본(aSAMt·AlphaFlow)은 넣지 않음: 모형이 도메인 전체 서열에서 만든 사슬이라 빠진 잔기가 없고, 그 기하는 재려는 모형 오차의 일부.
도메인·사슬마다 최대 Cα–Cα (반복 가운데 최대)와 4.2 Å 를 넘는 결합 수. 10-28 점검에 든 도메인(줄기 = derived/mdcath_5T, H-풀 = derived/r38/p3/dom 의 sbar 있는 npz, C = derived/r38/c/dom 의 sbar 있는 npz)과의 겹침을 셈.
판정(검산): 좌표를 찾은 모든 도메인·사슬에서 긴 결합 0 → holds, 아니면 fails (도메인·집합·긴 결합 수·최대 길이를 남김). 좌표를 찾지 못한 것은 missing 으로 따로 셈.
사용: python r49_chainbreak.py <ROOT>   출력: results/r49_chainbreak.json"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
import sys, json, glob, time, zipfile, difflib, tempfile
from multiprocessing import Pool
import numpy as np, pandas as pd, h5py
ROOT = sys.argv[1]; RES, DER = ROOT + '/results', ROOT + '/derived'; FO = RES + '/r49_chainbreak.json'
THR, NW = 4.2, 8
AA = dict(ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E', GLY='G', HIS='H', HSD='H', HSE='H', HSP='H', HIE='H', HID='H', HIP='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F', PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V')


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)
def md_atoms(pdbtxt):  # r37_ff.md_atoms 그대로
    res = {}; order = []
    for i, l in enumerate([l for l in pdbtxt.splitlines() if l.startswith(('ATOM', 'HETATM'))]):
        nm, rn, rs = l[12:16].strip(), l[17:21].strip(), (l[21], l[22:27].strip())
        if rs not in res: res[rs] = {'rn': rn}; order.append(rs)
        res[rs][nm] = i
    keep = [r for r in order if all(a in res[r] for a in ('CA', 'C', 'O'))]
    return ''.join(AA.get(res[r]['rn'], 'X') for r in keep), np.array([[res[r]['CA'], res[r]['C'], res[r]['O']] for r in keep])
def at_atoms(top):  # r37_ff.at_atoms 그대로
    keep = [r for r in top.residues if r.is_protein and all(any(a.name == n for a in r.atoms) for n in ('CA', 'C', 'O'))]
    idx = np.array([[next(a.index for a in r.atoms if a.name == n) for n in ('CA', 'C', 'O')] for r in keep])
    return ''.join(AA.get(r.name, 'X') for r in keep), idx
def block(s1, s2):  # r37_ff.block 그대로
    m = max(difflib.SequenceMatcher(None, s1, s2, autojunk=False).get_matching_blocks(), key=lambda b: b.size); return m.a, m.b, m.size
def bonds(cas):
    """CA 좌표 목록(반복마다 (L, 3)) → (최대 결합 길이, 4.2 Å 넘는 결합 수의 반복 가운데 최대)."""
    ls = [np.linalg.norm(np.diff(np.asarray(c, float), axis=0), axis=1) for c in cas]
    return float(max(l.max() for l in ls)), int(max(int((l > THR).sum()) for l in ls))


def ca_h5(d, a0=0, L=None):
    with h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r') as h:
        g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode())
        if L is not None: imd = imd[a0:a0 + L]
        return [g['320'][r]['coords'][0][imd[:, 0]] for r in sorted(g['320'].keys(), key=int)]
def one_md(d, npz_first=False):
    try:
        f = f'{ROOT}/raw/mdcath_p1n/{d}.npz'
        if npz_first and os.path.exists(f):
            z = np.load(f); cas = [z[f'h320_{int(r)}'][0, :, 0] for r in z['reps320']]; src = 'p1n_npz'
        elif os.path.exists(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5'):
            cas = ca_h5(d); src = 'h5'
        else:
            return d, dict(src='missing')
        mx, nl = bonds(cas); return d, dict(src=src, max=mx, n_long=nl, n_rep=len(cas))
    except Exception as e:
        return d, dict(src='error', err=repr(e)[:200])
def one_atlas(cd):
    c, d = cd
    try:
        import mdtraj as md
        with h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r') as h: smd, _ = md_atoms(h[d]['pdbProteinAtoms'][()].decode())
        z = zipfile.ZipFile(f'{ROOT}/raw/atlas/{c}_protein.zip'); pdbn = [n for n in z.namelist() if n.endswith('.pdb')][0]
        with tempfile.TemporaryDirectory() as tmp:
            z.extract(pdbn, tmp); t = md.load(f'{tmp}/{pdbn}')
        sat, iat = at_atoms(t.topology); a0, b0, L = block(smd, sat)
        at = bonds([t.xyz[0][iat[b0:b0 + L, 0]] * 10]); mdc = bonds(ca_h5(d, a0, L))
        return c, dict(dom=d, L=int(L), atlas_max=at[0], atlas_n_long=at[1], md_max=mdc[0], md_n_long=mdc[1], max=max(at[0], mdc[0]), n_long=max(at[1], mdc[1]))
    except Exception as e:
        return c, dict(dom=d, src='error', err=repr(e)[:200])


def main():
    t0 = time.time()
    model = sorted(os.path.basename(f)[:-5] for f in glob.glob(DER + '/r38/compare_D/*.json'))
    T = pd.read_csv(RES + '/r37_trunk_pairs.tsv', sep='\t'); b6_all = sorted(set(T.a) | set(T.b))
    J = [(os.path.basename(f)[:-5].split('__'), json.load(open(f))) for f in sorted(glob.glob(DER + '/r38/b6/*__*.json'))]; J = [ab for ab, j in J if 'skip' not in j and j.get('D')]
    b6_used = sorted({x for ab in J for x in ab})
    r10 = json.load(open(RES + '/r37_trunk_pairs.json'))['r10a_prime']['main']; r10a = sorted({x['a'] for x in r10} | {x['b'] for x in r10})
    bnbp = [json.loads(l)['ab'] for l in open(DER + '/r38/static/bnb.jsonl') if l.strip()]; bnb = sorted({x for p in bnbp for x in p})
    p1 = sorted(json.load(open(DER + '/r38/p1/pairs.json'))['domains']); p1n = sorted(json.load(open(DER + '/r38/p1n/pairs.json'))['domains'])
    dff = json.load(open(RES + '/r41_dff50.json')); atl = [(r['chain'], r['dom']) for r in dff['rows'] if r.get('d2n') is not None]
    log('r41_dff50 n', dff.get('n'), 'rows', len(dff['rows']), '쓴 사슬', len(atl))
    G = pd.read_csv(RES + '/v_residue_unit_graded.tsv', sep='\t'); stat = sorted(set(G.A) | set(G.B))
    sets = dict(model=model, b6=b6_all, b6_used=b6_used, r10a=r10a, bnb=bnb, p1=p1, p1n=p1n)
    log('집합', {k: len(v) for k, v in sets.items()}, 'atlas', len(atl), 'static', len(stat), '| b6 쌍', len(T), '분석', len(J), '| bnb 쌍', len(bnbp), '| r10a 쌍', len(r10))
    # 10-28 점검과의 겹침
    stem = {os.path.basename(f)[:-4] for f in glob.glob(DER + '/mdcath_5T/*.npz')}
    sb = lambda p: {os.path.basename(f)[:-4] for f in glob.glob(p + '/*.npz') if 'sbar' in np.load(f).files}
    hp, cc = sb(DER + '/r38/p3/dom'), sb(DER + '/r38/c/dom'); old = stem | hp | cc
    log('10-28 집합', len(stem), len(hp), len(cc))
    # mdCATH 도메인(p1n 은 npz 우선)
    md_only = sorted(set(model) | set(b6_all) | set(r10a) | set(bnb) | set(p1)); nz = sorted(set(p1n))
    with Pool(NW) as pool:
        A = dict(pool.starmap(one_md, [(d, False) for d in md_only], chunksize=4)); log('mdCATH h5', len(A), round(time.time() - t0))
        B = dict(pool.starmap(one_md, [(d, True) for d in nz], chunksize=4)); log('p1n', len(B), round(time.time() - t0))
        C = dict(pool.map(one_atlas, atl, chunksize=1)); log('atlas', len(C), round(time.time() - t0))
    # 정적 사슬(hu1_rulers.bb 와 같은 읽기)
    sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
    import framespace as fs
    z = np.load(DER + '/xprot_structs.npz', allow_pickle=True); S = {}
    for nm in stat:
        try:
            k = nm + '|CA'
            if k in z.files: ca = z[k][0] if z[k].ndim == 3 else z[k]; src = 'xprot'
            else: ca = fs.read_pdb_backbone(DER + f'/atlas_pdb/{nm}.pdb')['CA'][0]; src = 'atlas_pdb'
            mx, nl = bonds([ca]); S[nm] = dict(src=src, max=mx, n_long=nl)
        except Exception as e:
            S[nm] = dict(src='error', err=repr(e)[:200])
    log('static', len(S), round(time.time() - t0))

    def summ(names, D, doms=None):
        v = [D[x] for x in names]; ok = [x for x in v if 'max' in x]; mxs = [x['max'] for x in ok]
        longs = [(n, D[n]['n_long'], round(D[n]['max'], 3)) for n in names if 'max' in D[n] and D[n]['n_long'] > 0]
        return dict(n=len(names), n_checked=len(ok), n_missing=sum(x.get('src') == 'missing' for x in v), n_error=sum(x.get('src') == 'error' for x in v),
                    n_with_long=len(longs), long=longs, max_median=float(np.median(mxs)) if mxs else None, max_max=float(max(mxs)) if mxs else None,
                    in_1028=None if doms is None else sum(x in old for x in doms))
    out = {k: summ(v, B if k == 'p1n' else A, v) for k, v in sets.items()}
    out['atlas'] = summ([c for c, _ in atl], C, [d for _, d in atl]); out['atlas']['atlas_side_max'] = max((x['atlas_max'] for x in C.values() if 'atlas_max' in x), default=None)
    out['atlas']['md_side_max'] = max((x['md_max'] for x in C.values() if 'md_max' in x), default=None)
    out['static'] = summ(stat, S)
    out['p1n']['src'] = {s: sum(B[d].get('src') == s for d in nz) for s in ('p1n_npz', 'h5', 'missing', 'error')}
    allD = {**{('md', k): v for k, v in A.items()}, **{('p1n', k): v for k, v in B.items()}, **{('atlas', k): v for k, v in C.items()}, **{('static', k): v for k, v in S.items()}}
    nlong = sum(1 for v in allD.values() if v.get('n_long', 0) > 0); nchk = sum(1 for v in allD.values() if 'max' in v)
    errs = {f'{a}|{b}': v.get('err', v.get('src')) for (a, b), v in allD.items() if v.get('src') in ('missing', 'error')}
    res = dict(written=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()), threshold_A=THR, sets=out, n_checked=nchk, n_with_long=nlong, not_checked=errs,
               old_1028=dict(stem=len(stem), hpool_sbar=len(hp), c_sbar=len(cc)),
               verdict='holds' if nlong == 0 else 'fails',
               rule='addendum12 r49_chainbreak: 좌표를 찾은 모든 도메인·사슬에서 Cα–Cα > 4.2 Å 결합 0 → holds, 아니면 fails. 찾지 못한 것은 not_checked 로 따로.',
               seconds=round(time.time() - t0))
    json.dump(res, open(FO, 'w'), indent=1, ensure_ascii=False); log('판정', res['verdict'], '| 점검', nchk, '| 긴 결합', nlong, '| 못 봄', len(errs)); log('R49_DONE')


if __name__ == '__main__':
    main()
