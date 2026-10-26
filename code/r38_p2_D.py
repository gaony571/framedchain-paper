"""(prereg round36 7c P2 의 D 판; 고침 round36_amend_p2_D_2026-10-19, 17:07:53 UTC) aSAMt 320 K 대 MD 320 K 와 MD 온도 반응을 D 로 다시 낸다.
양: Δμ²₂_D = (C₂ − (B₂_P + B₂_Q)/2)/변 (지수 2 에너지 통계; 맞춘 좌표에서는 등록된 Δμ² 와 같음), Δ²_D/√변 (2차 1),
    같은 표본의 맞춘 좌표 Δμ²(지수 2 꼴)·Δ²/√변 (비교용).
창 = 전체(aSAMt 학습 창), Q ≥ 0.75. MD 실행 = 반복(N = 100; r38_s0 의 이름 seed 표본), MD 교차·바닥 = 실행 쌍마다 프레임 쌍 2,000 (줄인 계산).
모델 = derived/r38/asamt/<dom>.320.traj.dcd 의 접힌 표본(MD 첫 프레임 천연 접촉 Q ≥ 0.75), 서열·길이가 MD 와 같을 때만(r38_compare 와 같은 점검).
모델 실행 = 무작위 반쪽 둘(각 ≥ 40), 20 나눔 평균(r38_compare 와 같은 구성). 모델 쪽 쌍(모델 × MD 320, 모델 × 모델)은 모든 프레임 쌍을 한 번 내고 나눔마다 부분 평균.
귀무 = MD 320 K 반복 2 대 3 (10 가지).
출력: derived/r38/compare_D/<dom>.json — compare 와 같은 칸 이름(aSAMt320-MD320|full|0.75, MD{T}-MD320|full|0.75, null320d2n|full|0.75, null320|full|0.75),
      dmu = Δμ²₂_D/변, d2n = Δ²_D/√변 (그래서 동결된 r38_p2_aggregate.py 를 입력 경로만 바꿔 쓴다: r38_p2_aggregate_D.py). 맞춘 좌표 값은 dmu_A·d2n_A.
사용: python r38_p2_D.py <ROOT> <도메인 목록 파일>"""
import sys, os, json, glob, time, itertools
ROOT = sys.argv[1]; LIST = sys.argv[2]
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)
import mdtraj as md
OUTD = ROOT + '/derived/r38/compare_D'; os.makedirs(OUTD, exist_ok=True)
DOMS = open(LIST).read().split(); W = 'full'; QK = '0.75'; MINH = 40; NSPLIT = 20


def seq_of(d):
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); s, _ = md_atoms(h[d]['pdbProteinAtoms'][()].decode()); h.close(); return s


def model_frames(d, Qf, L, smd):
    f = f'{ROOT}/derived/r38/asamt/{d}.320'
    if not os.path.exists(f + '.traj.dcd'): return None, 'aSAMt 없음'
    t = md.load(f + '.traj.dcd', top=f + '.top.pdb'); sm, im = at_atoms(t.topology); a0, b0, Lb = block(smd, sm)
    if Lb != L or sm != smd: return None, '서열 또는 길이 다름'
    X = t.xyz[:, im.ravel()].reshape(t.n_frames, L, 3, 3) * 10
    return X[np.where(Qf(X) >= QT)[0]], None


def stats_of(cr, fp, fq, ne):
    """cr, fp, fq: 실행 쌍마다 dict(d, d2, a, a2) 목록. P = 앞쪽(MD 320), Q = 뒤쪽."""
    m = lambda xs, k: float(np.mean([v[k] for v in xs]))
    D = estats(m(cr, 'd'), m(fp, 'd'), m(fq, 'd'), ne); Aa = estats(m(cr, 'a'), m(fp, 'a'), m(fq, 'a'), ne)
    return dict(dmu=(m(cr, 'd2') - (m(fp, 'd2') + m(fq, 'd2')) / 2) / ne, d2n=D['d2n'], D2=D['D2'], C=D['C'], B_P=D['B_P'], B_Q=D['B_Q'], wbar=D['wbar'], ratio=D['ratio'],
                dmu_A=(m(cr, 'a2') - (m(fp, 'a2') + m(fq, 'a2')) / 2) / ne, d2n_A=Aa['d2n'])


def one(d):
    rng = seeded('r38_p2_D', d); _runs.clear(); _md.clear()
    MD, Qf, L = md_load(d, TT); ne = L - 1; smd = seq_of(d)
    ref = torch.as_tensor(np.load(f'{ROOT}/derived/r38/md/{d}.npz')['ref'])
    R = {T: dom_runs(d, T, W, TT)[0] for T in TT}
    res = dict(dom=d, L=L, ne=ne, nruns={f'MD{T}|{W}|{QK}': len(R[T]) for T in TT}, cells={}, window=W, q=QK,
               quantity='dmu = Δμ²₂_D/변 (지수 2 에너지 통계, D), d2n = Δ²_D/√변; dmu_A·d2n_A = 같은 표본의 맞춘 좌표')
    if len(R['320']) < 2: res['skip'] = '320 K 실행 < 2'; return res
    Ts = [T for T in TT if len(R[T]) >= 2]
    E = {T: coords(np.concatenate([X for _, X in R[T]])) for T in Ts}
    A = {T: fs.align_to_reference(frames(np.concatenate([X for _, X in R[T]])), ref).numpy().reshape(len(R[T]) * N, -1) for T in Ts}
    def go(Ta, Tb, rr):
        o = runpair_means(MeanStart(E[Ta], E[Tb]), rr, rng); ea2 = ((A[Ta][o['IA']] - A[Tb][o['IB']]) ** 2).sum(1); out = {}
        for g, rp in enumerate(rr):
            s_ = o['GR'] == g; out[rp] = dict(d=float(o['dv'][s_].mean()), d2=float((o['dv'][s_] ** 2).mean()), a=float(np.sqrt(ea2[s_]).mean()), a2=float(ea2[s_].mean()))
        return out
    n320 = len(R['320']); f320 = go('320', '320', list(itertools.combinations(range(n320), 2)))
    for T in Ts[1:]:
        nT = len(R[T]); fT = go(T, T, list(itertools.combinations(range(nT), 2))); cr = go('320', T, list(itertools.product(range(n320), range(nT))))
        res['cells'][f'MD{T}-MD320|{W}|{QK}'] = stats_of(list(cr.values()), list(f320.values()), list(fT.values()), ne)
    if n320 == 5:
        ff = lambda a, b: f320[(min(a, b), max(a, b))]; nl = []
        for g2 in itertools.combinations(range(5), 2):
            g3 = [k for k in range(5) if k not in g2]
            nl.append(stats_of([ff(a, b) for a in g2 for b in g3], [ff(*g2)], [ff(a, b) for a, b in itertools.combinations(g3, 2)], ne))
        res['cells'][f'null320d2n|{W}|{QK}'] = [x['d2n'] for x in nl]; res['cells'][f'null320|{W}|{QK}'] = [x['dmu'] * ne for x in nl]
        res['cells'][f'null320dmu|{W}|{QK}'] = [x['dmu'] for x in nl]; res['cells'][f'null320dmu_A|{W}|{QK}'] = [x['dmu_A'] for x in nl]
    Xm, why = model_frames(d, Qf, L, smd)
    if Xm is None or len(Xm) < 2 * MINH:
        res['model_skip'] = why or f'접힌 표본 {len(Xm)} < {2 * MINH}'; return res
    EM = coords(Xm); AM = fs.align_to_reference(frames(Xm), ref).numpy().reshape(len(Xm), -1); nm = len(Xm); E320, A320 = E['320'], A['320']
    ia, ib = np.meshgrid(np.arange(nm), np.arange(len(E320)), indexing='ij')
    dX = red_pairs(MeanStart(EM, E320), ia.ravel(), ib.ravel())[0].reshape(nm, len(E320)); aX2 = cost(AM, A320)
    iu, ju = np.triu_indices(nm, 1); dm = red_pairs(MeanStart(EM, EM), iu, ju)[0]; DM = np.zeros((nm, nm)); DM[iu, ju] = dm; DM = DM + DM.T; aM2 = cost(AM, AM)
    fP = list(f320.values()); outs = []
    blk = lambda r: np.arange(r * N, (r + 1) * N)
    one_ = lambda dd, aa: dict(d=float(dd.mean()), d2=float((dd ** 2).mean()), a=float(np.sqrt(np.maximum(aa, 0)).mean()), a2=float(aa.mean()))
    for s in range(NSPLIT):
        perm = rng.permutation(nm); h = nm // 2; H = [perm[:h], perm[h:2 * h]]
        cr = [one_(dX[np.ix_(hh, blk(r))], aX2[np.ix_(hh, blk(r))]) for hh in H for r in range(n320)]
        fq = [one_(DM[np.ix_(H[0], H[1])], aM2[np.ix_(H[0], H[1])])]
        outs.append(stats_of(cr, fP, fq, ne))
    c = {k: float(np.mean([o[k] for o in outs])) for k in outs[0]}; c.update(n_model=int(nm), nsplit=NSPLIT)
    res['cells'][f'aSAMt320-MD320|{W}|{QK}'] = c
    return res


for d in DOMS:
    fo = f'{OUTD}/{d}.json'
    if os.path.exists(fo): continue
    try:
        t0 = time.time(); r = one(d); r['sec'] = time.time() - t0; jdump(r, fo); log('p2_D', d, len(r['cells']), '칸', round(r['sec']), 's', r.get('model_skip', ''), r.get('skip', ''))
    except Exception as e:
        traceback.print_exc(); log('ERR p2_D', d, repr(e)[:200])
log('P2_D_DONE')
