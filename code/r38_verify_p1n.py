"""(prereg 7h) r38_p1n.py 구현 검증. 새 자료의 판정 값은 보지 않는다(1차 P1 자료와 시험용 사본만 씀). 출력: results/r38_verify_p1n.json
사용: CUDA_VISIBLE_DEVICES= python r38_verify_p1n.py <ROOT>"""
import sys, os, json, time, traceback, shutil, hashlib
ROOT = sys.argv[1]
NS = {'__name__': 'vp1n'}; _s = open(ROOT + '/code/r38_p1n.py').read(); _s = _s[:_s.index("ALLN = [")]
sys.argv = ['r38_p1n', ROOT, 'none', '2']; exec(_s, NS)
np, h5py = NS['np'], NS['h5py']; res = {}
def check(name, ok, info=''):
    res[name] = dict(ok=bool(ok), info=info); print(('통과 ' if ok else '실패 ') + name, info, flush=True)

try:
    pj = json.load(open(NS['OUT1'] + '/pairs.json')); n = pj['n']; md5 = hashlib.md5(open(NS['OUT1'] + '/pairs.json', 'rb').read()).hexdigest()
    check('V0 쌍 목록 = 등록(483·59·1,014·35·531·311, md5 b2d5f22b)', (n['HT'], n['HT_sf'], n['MM'], n['sf_both'], n['domains'], n['download']) == (483, 59, 1014, 35, 531, 311) and md5.startswith('b2d5f22b'), f'{n} {md5[:8]}')
    # V1 npz 입력 = h5 입력 (1차 P1 의 h5 도메인 하나를 받은 꼴로 잘라 npz 로)
    old = json.load(open(ROOT + '/derived/r38/p1/pairs.json'))
    d = next(x for x in sorted(old['domains']) if os.path.exists(f'{ROOT}/derived/r38/p1/mean/{x}.npz') and 'm' in np.load(f'{ROOT}/derived/r38/p1/mean/{x}.npz').files) # 칸이 있는 도메인
    tmp = '/content/vp1n'; shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = NS['md_atoms'](g['pdbProteinAtoms'][()].decode()); L = len(smd); out = dict(smd=np.array(list(smd)))
    for T in ('320', '348', '379', '413'):
        if T not in g: continue
        reps = sorted(g[T].keys(), key=int); nf = []
        for k, r in enumerate(reps):
            ds = g[T][r]['coords']; nn = ds.shape[0]; nf.append(nn)
            if T == '320' and k == 0: out['first'] = ds[0][imd.ravel()].reshape(L, 3, 3).astype(np.float32)
            out[f'h{T}_{r}'] = ds[nn // 2:][:, imd.ravel()].reshape(-1, L, 3, 3).astype(np.float32)
            if nn >= 100: out[f'a{T}_{r}'] = ds[50:100][:, imd.ravel()].reshape(-1, L, 3, 3).astype(np.float32)
        out[f'nframes{T}'] = np.array(nf); out[f'reps{T}'] = np.array([int(r) for r in reps])
    h.close(); np.savez(f'{tmp}/{d}.npz', **out)
    keep = NS['RAWN']; NS['RAWN'] = tmp
    try: c1, l1, s1 = NS['build_cells_npz'](d); NS['OUT1'] = tmp; os.makedirs(tmp + '/fold', exist_ok=True); NS['_fold_worker'](d); f1 = json.load(open(f'{tmp}/fold/{d}.json'))['f']
    finally: NS['RAWN'] = keep; NS['OUT1'] = ROOT + '/derived/r38/p1n'
    c0, l0, s0 = NS['build_cells_h5'](d)
    same = sorted(c1) == sorted(c0) and all(np.array_equal(c1[k][0], c0[k][0]) and np.array_equal(c1[k][1], c0[k][1]) for k in c0) and l1 == l0 and s1 == s0
    check('V1 npz 입력 칸 = h5 입력 칸 (모든 칸·같은 프레임·같은 라벨·궤적 길이)', same and len(c0) >= 4, f'{d}, 칸 {len(c0)}')
    NS['RAWN'] = '/nonexistent'; os.makedirs('/content/vp1n_h/fold', exist_ok=True); NS['OUT1'] = '/content/vp1n_h'
    try: NS['_fold_worker'](d); f0 = json.load(open(f'/content/vp1n_h/fold/{d}.json'))['f']
    finally: NS['RAWN'] = keep; NS['OUT1'] = ROOT + '/derived/r38/p1n'
    check('V2 접힌 비율: npz = h5', f1 == f0, str(f1.get('379')))
    # V3 ustats 를 1차 P1 자료에 돌려 결과 뒤 진단(r38_p1_diag.json)의 점추정을 재현
    NS['OUT1'] = ROOT + '/derived/r38/p1'; SPo = NS['load_spread'](old['domains']); NS['OUT1'] = ROOT + '/derived/r38/p1n'
    U = NS['ustats'](old['HT'], old['MM'], SPo, np.random.default_rng(1), B=50); dg = json.load(open(ROOT + '/results/r38_p1_diag.json'))['unc']
    pairs_ = [('half|0.70', 'q070'), ('half|0.75', 'main'), ('half|0.80', 'q080'), ('abs|0.75', 'abs_50_100ns')]
    ok3 = all(abs(U[a]['U1']['median'] - dg[b]['U1']['median']) < 1e-12 and abs(U[a]['U2']['median'] - dg[b]['U2']['median']) < 1e-12 for a, b in pairs_)
    check('V3 ustats 점추정 = 결과 뒤 진단 값 (1차 자료, 네 칸)', ok3, f"U2(0.7) {U['half|0.70']['U2']['median']:.4f} / 진단 {dg['q070']['U2']['median']:.4f}; 분해 차 {U['decomp']['diff']:.4f}")
    # V4 병렬 정렬 작업자 = align_batch 직접 (1차 P1 평균 점, 쌍 넷)
    ps = [(r['a'], r['b']) for r in old['HT'][:4]]; Ms = {}
    for p in ps:
        for x in p:
            z = np.load(f'{ROOT}/derived/r38/p1/mean/{x}.npz'); Ms[x] = z['m'] if 'm' in z.files else None
    ps = [p for p in ps if Ms[p[0]] is not None and Ms[p[1]] is not None]; NS['MSA'].update(Ms)
    a = NS['_align_worker'](ps); b = NS['align_batch'](ps, Ms, 'cpu'); al = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(ROOT + '/derived/r38/p1/align.jsonl')}
    ok4 = all(abs(x[1][0] - y[0]) < 1e-9 and x[1][1] == y[1] for x, y in zip(a, b)) and all(abs(x[1][0] - al[tuple(x[0])]['S']) < 1e-6 * abs(al[tuple(x[0])]['S']) for x in a if tuple(x[0]) in al)
    check('V4 병렬 정렬 = align_batch 직접 = 1차 P1 align.jsonl', ok4, f'{len(ps)} 쌍')
except Exception as e:
    traceback.print_exc(); check('스크립트 오류', False, repr(e)[:300])
json.dump(dict(written=time.strftime('%F %T'), results=res, n_ok=sum(v['ok'] for v in res.values()), n=len(res)), open(ROOT + '/results/r38_verify_p1n.json', 'w'), ensure_ascii=False, indent=1)
print('VERIFY_P1N_DONE', sum(v['ok'] for v in res.values()), '/', len(res))
