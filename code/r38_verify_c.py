"""확인 실험 C (r38_c.py) 구현 검증. 판정 값(C-1·C-2 의 편 Spearman)은 보지 않는다. 출력: results/r38_verify_c.json
사용: CUDA_VISIBLE_DEVICES= python r38_verify_c.py <ROOT>"""
import sys, os, json, time, traceback, shutil
ROOT = sys.argv[1]
CN = {'__name__': 'cv'}; _s = open(ROOT + '/code/r38_c.py').read(); _s = _s[:_s.index("ALLC = [")]
sys.argv = ['r38_c', ROOT, 'none', '1']; exec(_s, CN)
np, torch, fs, h5py = CN['np'], CN['torch'], CN['fs'], CN['h5py']
res = {}
def check(name, ok, info=''):
    res[name] = dict(ok=bool(ok), info=info); print(('통과 ' if ok else '실패 ') + name, info, flush=True)

try:
    pj = json.load(open(CN['OUT'] + '/pairs.json'))
    check('V0 쌍 목록 = 등록 수(4,395 쌍·268 상과·1,014 도메인)', pj['n_pairs'] == 4395 and pj['n_sf'] == 268 and pj['n_domains'] == 1014, f"{pj['n_pairs']} {pj['n_sf']} {pj['n_domains']}")
    # V1 npz 입력 = h5 입력 (H-풀의 작은 도메인 하나를 받은 꼴과 같게 잘라 npz 로)
    hp = json.load(open(CN['P3OUT'] + '/pairs.json'))['domains']; d = None
    for x in hp:
        f = f"{CN['P3OUT']}/dom/{x}.npz"
        if os.path.exists(f) and 'm' in np.load(f).files: d = x; break
    tmp = '/content/cverify'; shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = CN['md_atoms'](g['pdbProteinAtoms'][()].decode()); L = len(smd); out = dict(smd=np.array(list(smd)))
    reps = sorted(g['320'].keys(), key=int)
    for k, r in enumerate(reps):
        ds = g['320'][r]['coords']; n = ds.shape[0]
        if k == 0: out['first'] = ds[0][imd.ravel()].reshape(L, 3, 3).astype(np.float32)
        out[f'run{r}'] = ds[n // 2:][:, imd.ravel()].reshape(-1, L, 3, 3).astype(np.float32)
    h.close(); out['reps'] = np.array([int(r) for r in reps]); np.savez(f'{tmp}/{d}.npz', **out)
    keep = CN['RAWC']; CN['RAWC'] = tmp
    try: r1, i1 = CN['dom_build_npz'](d)
    finally: CN['RAWC'] = keep
    r0, i0 = CN['dom_build_h5'](d)
    same = all(np.array_equal(r1[k], r0[k]) for k in ('A', 'rep', 'degen')) and np.abs(r1['m'] - r0['m']).max() < 1e-12 and i1 == i0
    check('V1 npz 입력 앙상블 = h5 입력 앙상블 (같은 프레임·같은 seed)', same, f'{d}, 프레임 {len(r0["A"])}')
    # V2 ANM 성질
    med = np.load(f"{CN['P3OUT']}/dom2/{d}.npz")['medoid'].astype(np.float64); ca = np.asarray(fs.ca_trace(torch.as_tensor(med))); p0, w8 = CN['anm_edge_spread'](ca)
    q = np.random.default_rng(0).normal(size=4); q /= np.linalg.norm(q); a, b, c, e = q
    Rm = np.array([[1 - 2 * (c * c + e * e), 2 * (b * c - e * a), 2 * (b * e + c * a)], [2 * (b * c + e * a), 1 - 2 * (b * b + e * e), 2 * (c * e - b * a)], [2 * (b * e - c * a), 2 * (c * e + b * a), 1 - 2 * (b * b + c * c)]])
    p1, _ = CN['anm_edge_spread'](ca @ Rm.T + np.array([3.0, -7.0, 11.0]))
    check('V2 ANM: 0 고유값 6 개·회전·평행이동 불변·양수·길이', np.abs(w8[:6]).max() < 1e-6 * w8[6] and w8[6] > 1e-8 and np.abs(p0 - p1).max() < 1e-6 * p0.max() and (p0 > 0).all() and len(p0) == len(med),
          f'w6 {w8[6]:.2e}, 0 고유값 최대 {np.abs(w8[:6]).max():.1e}, 회전 차 {np.abs(p0 - p1).max():.1e}')
    md = CN['gpa_profile'](r0['A'].astype(np.float64), np.arange(r0['A'].shape[1])); rho = CN['spearman'](p0, md)
    check('V3 ANM 예측 윤곽과 MD 윤곽의 Spearman (바름 점검, 양수)', rho > 0, f'{d}: {rho:.3f}')
    # V4 analyse_sf 점추정 = r38_p3.analyse 점추정, 상과만 구간이 등록 1b 구간보다 좁음
    rows = [json.loads(l) for l in open(CN['P3OUT'] + '/pairs_out.jsonl')]; CN['NB'] = 200
    a_sf = CN['analyse_sf'](rows, seed=1, B=200); a_1b = CN['analyse_sf'](rows, seed=1, B=200, sf_only=False); a_p3 = CN['analyse'](rows, boot=False)
    check('V4 analyse_sf 점추정 = r38_p3.analyse, 상과만 구간 < 1b 구간', abs(a_sf['est'] - a_p3['est']) < 1e-12 and (a_sf['hi'] - a_sf['lo']) < (a_1b['hi'] - a_1b['lo']),
          f"폭 {a_sf['hi'] - a_sf['lo']:.3f} 대 {a_1b['hi'] - a_1b['lo']:.3f}")
    # V5 anm_pairs 의 공유 변 = r38_p3 pairs_out 의 공유 변 (H-풀 쌍 몇 개)
    po = {tuple(r['ab']): r for r in rows}; sub = [p for p in json.load(open(CN['P3OUT'] + '/pairs.json'))['pairs'] if tuple(p[:2]) in po][:6]
    t5 = '/content/cverify5'; shutil.rmtree(t5, ignore_errors=True); os.makedirs(t5 + '/dom')
    json.dump(dict(pairs=sub), open(t5 + '/pairs.json', 'w'))
    with open(t5 + '/pibar.jsonl', 'w') as fh:
        for l in open(CN['P3OUT'] + '/pibar.jsonl'):
            if tuple(json.loads(l)['ab']) in {tuple(p[:2]) for p in sub}: fh.write(l)
    for p in sub:
        for x in p[:2]: shutil.copy(f"{CN['P3OUT']}/dom/{x}.npz", f'{t5}/dom/{x}.npz')
    CN['anm_pairs'](t5, t5 + '/anm.jsonl'); an = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(t5 + '/anm.jsonl')}
    ok5 = [an.get(tuple(p[:2]), {}).get(tag, {}).get('n_shared') == po[tuple(p[:2])][tag]['n_shared'] # 30 변 이상인 판만 anm_pairs 가 남긴다
           for p in sub for tag in ('main', 'inter') if po[tuple(p[:2])][tag].get('r_dyn') is not None]
    check('V5 anm_pairs 공유 변 = r38_p3 one_pair 공유 변 (주·교집합)', len(ok5) > 0 and all(ok5), str(ok5))
    # V7 C 의 pibar(묶음 크기만 다름): 묶음이 같으면 r38_p3 pibar 와 같은 출력, 묶음 1 이면 S 가 EPS 안
    t7 = '/content/cverify7'; shutil.rmtree(t7, ignore_errors=True); keep = CN['OUT']; txt = {}
    try:
        for tag, fn in (('p3', lambda: CN['part_pibar_p3']()), ('cbig', lambda: CN['part_pibar'](cap=1e12)), ('c1', lambda: CN['part_pibar'](cap=1))):
            os.makedirs(f'{t7}/{tag}/dom'); json.dump(dict(pairs=sub[:4]), open(f'{t7}/{tag}/pairs.json', 'w'))
            for p in sub[:4]:
                for x in p[:2]: shutil.copy(f"{CN['P3OUT']}/dom/{x}.npz", f'{t7}/{tag}/dom/{x}.npz')
            CN['OUT'] = f'{t7}/{tag}'; fn(); txt[tag] = open(f'{t7}/{tag}/pibar.jsonl').read()
    finally: CN['OUT'] = keep
    S_ = lambda s: {tuple(json.loads(l)['ab']): json.loads(l)['S'] for l in s.splitlines()}
    s0, s1 = S_(txt['p3']), S_(txt['c1']); rel = max(abs(s1[k] - s0[k]) / abs(s0[k]) for k in s0)
    check('V7 C pibar: 같은 묶음 = r38_p3 출력, 묶음 1 의 S 상대 차 < EPS', txt['p3'] == txt['cbig'] and len(s0) == 4 and rel < CN['EPS'], f"같은 묶음 동일 {txt['p3'] == txt['cbig']}, 묶음 1 상대 차 최대 {rel:.1e}")
    B = CN['blosum']()
    check('V6 BLOSUM62 대칭·같은 잔기 양수·알려진 값', B is not False and B['A', 'R'] == B['R', 'A'] == -1 and B['W', 'W'] == 11, 'biopython 없음' if B is False else f"A/R {B['A', 'R']}, W/W {B['W', 'W']}")
except Exception as e:
    traceback.print_exc(); check('스크립트 오류', False, repr(e)[:300])
json.dump(dict(written=time.strftime('%F %T'), results=res, n_ok=sum(v['ok'] for v in res.values()), n=len(res)), open(ROOT + '/results/r38_verify_c.json', 'w'), ensure_ascii=False, indent=1)
print('VERIFY_C_DONE', sum(v['ok'] for v in res.values()), '/', len(res))
