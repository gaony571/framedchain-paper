"""(prereg 7c 2차 2·3 의 D 판; 7d 로 P2 주 판이 D 이므로) AlphaFlow 증류판 대 MD 320 K 를 D 로. r38_p2_D.py 의 모델 칸과 같은 계산이다:
모델 표본 = derived/r38/alphaflow/<dom>.pdb 의 접힌 표본(Q ≥ 0.75, MD 첫 프레임 천연 접촉), 서열·길이가 MD 와 같을 때만; 무작위 반쪽 둘(각 ≥ 40), 20 나눔 평균;
MD 320 K 실행 = r38_s0 의 이름 seed 표본(N = 100), MD 바닥 = 실행 쌍마다 프레임 쌍 2,000 (줄인 계산). 창 = 전체, Q = 0.75. 쌍 고르기·나눔 seed = seeded('r38_p2_D_af', 도메인).
출력: derived/r38/compare_D2/<dom>.json = compare_D/<dom>.json 에 'AlphaFlow-MD320|full|0.75' 칸을 더한 사본(compare_D 는 그대로).
집계: r38_p2_aggregate_D2.py (동결 집계를 입력 compare_D2 와 d_ff D 판 results/r38_dff_D.json 으로만 바꿔 돌림).
사용: python r38_p2_D_af.py <ROOT> <도메인 목록 파일>"""
import sys, os, json, glob, time, itertools
ROOT = sys.argv[1]; LIST = sys.argv[2]
_src = open(ROOT + '/code/r38_p2_D.py').read(); _src = _src[:_src.index('\nfor d in DOMS:')]
sys.argv = ['r38_p2_D', ROOT, LIST]; exec(_src)
OUTD2 = ROOT + '/derived/r38/compare_D2'; os.makedirs(OUTD2, exist_ok=True)


def af_frames(d, Qf, L, smd):
    f = f'{ROOT}/derived/r38/alphaflow/{d}.pdb'
    if not os.path.exists(f): return None, 'AlphaFlow 없음'
    t = md.load(f); sm, im = at_atoms(t.topology); a0, b0, Lb = block(smd, sm)
    if Lb != L or sm != smd: return None, '서열 또는 길이 다름'
    X = t.xyz[:, im.ravel()].reshape(t.n_frames, L, 3, 3) * 10
    return X[np.where(Qf(X) >= QT)[0]], None


def one_af(d):
    rng = seeded('r38_p2_D_af', d); _runs.clear(); _md.clear()
    MD, Qf, L = md_load(d, TT); ne = L - 1; smd = seq_of(d)
    ref = torch.as_tensor(np.load(f'{ROOT}/derived/r38/md/{d}.npz')['ref'])
    R320 = dom_runs(d, '320', W, TT)[0]; out = dict(n_md320=len(R320))
    if len(R320) < 2: out['skip'] = '320 K 실행 < 2'; return out
    E320 = coords(np.concatenate([X for _, X in R320])); A320 = fs.align_to_reference(frames(np.concatenate([X for _, X in R320])), ref).numpy().reshape(len(R320) * N, -1)
    n320 = len(R320); rr = list(itertools.combinations(range(n320), 2))
    o = runpair_means(MeanStart(E320, E320), rr, rng); ea2 = ((A320[o['IA']] - A320[o['IB']]) ** 2).sum(1)
    fP = [dict(d=float(o['dv'][o['GR'] == g].mean()), d2=float((o['dv'][o['GR'] == g] ** 2).mean()), a=float(np.sqrt(ea2[o['GR'] == g]).mean()), a2=float(ea2[o['GR'] == g].mean())) for g in range(len(rr))]
    Xm, why = af_frames(d, Qf, L, smd)
    if Xm is None or len(Xm) < 2 * MINH:
        out['model_skip'] = why or f'접힌 표본 {len(Xm)} < {2 * MINH}'; return out
    EM = coords(Xm); AM = fs.align_to_reference(frames(Xm), ref).numpy().reshape(len(Xm), -1); nm = len(Xm)
    ia, ib = np.meshgrid(np.arange(nm), np.arange(len(E320)), indexing='ij')
    dX = red_pairs(MeanStart(EM, E320), ia.ravel(), ib.ravel())[0].reshape(nm, len(E320)); aX2 = cost(AM, A320)
    iu, ju = np.triu_indices(nm, 1); dm = red_pairs(MeanStart(EM, EM), iu, ju)[0]; DM = np.zeros((nm, nm)); DM[iu, ju] = dm; DM = DM + DM.T; aM2 = cost(AM, AM)
    blk = lambda r: np.arange(r * N, (r + 1) * N)
    one_ = lambda dd, aa: dict(d=float(dd.mean()), d2=float((dd ** 2).mean()), a=float(np.sqrt(np.maximum(aa, 0)).mean()), a2=float(aa.mean()))
    outs = []
    for s in range(NSPLIT):
        perm = rng.permutation(nm); h = nm // 2; H = [perm[:h], perm[h:2 * h]]
        cr = [one_(dX[np.ix_(hh, blk(r))], aX2[np.ix_(hh, blk(r))]) for hh in H for r in range(n320)]
        fq = [one_(DM[np.ix_(H[0], H[1])], aM2[np.ix_(H[0], H[1])])]
        outs.append(stats_of(cr, fP, fq, ne))
    c = {k: float(np.mean([x[k] for x in outs])) for k in outs[0]}; c.update(n_model=int(nm), nsplit=NSPLIT)
    out['cell'] = c; return out


if __name__ == '__main__':
    for d in DOMS:
        fo = f'{OUTD2}/{d}.json'
        if os.path.exists(fo): continue
        try:
            t0 = time.time(); base = json.load(open(f'{OUTD}/{d}.json')); r = one_af(d)
            base['af_D'] = {k: v for k, v in r.items() if k != 'cell'}
            if 'cell' in r: base['cells'][f'AlphaFlow-MD320|{W}|{QK}'] = r['cell']
            base['af_sec'] = time.time() - t0; jdump(base, fo); log('p2_D_af', d, r.get('cell', {}).get('dmu'), r.get('model_skip', ''), r.get('skip', ''), round(base['af_sec']), 's')
        except Exception as e:
            traceback.print_exc(); log('ERR p2_D_af', d, repr(e)[:200])
    log('P2_D_AF_DONE')
