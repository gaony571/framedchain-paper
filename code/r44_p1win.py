"""(round44, 등록 = Drive results/prereg_2026-10-28_addendum2.json "r44_p1win", 2026-10-28 UTC) 확인 실험 P1′ 의 U2 를 절대 시간 창 여럿에서.
P1′ (r38_p1n) 의 주 통계 U2(Q 0.7) 는 뒤 절반 창(반복마다 길이가 달라 절대 시각이 다름)이고, 50–100 ns 창(Q 0.75)에서는 −0.033 으로 작았다(외부 검토).
같은 쌍·같은 공유 결합(derived/r38/p1n/align.jsonl 의 best)에서 320 K 의 절대 창 셋을 낸다: 50–100, 300–400, 400–500 ns (프레임 = 1 ns).
프레임: raw/mdcath_p1n/<d>.npz 가 있으면 그 뒤 절반(h320_r, 절대 시작 = 길이//2)과 50–100 ns(a320_r), 없으면 raw/mdcath h5 (전체).
창 안이 모두 궤적 안이어야 하고(창 끝 ≤ 길이), 창이 npz 의 뒤 절반 안에 있어야 한다. 창마다 Q ≥ 0.7 인 접힌 프레임 ≥ 50 인 반복에서 50 개(고정 seed), 반복 ≥ 2.
퍼짐 = r38_p1 과 같은 gpa_profile 합(공유 결합 위). U2 = 상과마다 [호열–중온 쌍의 log(퍼짐_호열/퍼짐_중온) 중앙 − 중온–중온 대조 쌍의 같은 값 중앙] 의 상과 중앙.
구간 = 상과 재추출 2,000 (seed 20261096), 한쪽 97.5% 상한. 세 창 모두 값이 있는 같은 상과 집합 판도 낸다.
부분: spread (CPU 작업자) → stats   사용: python r44_p1win.py <ROOT> spread,stats [작업자=4]   출력: derived/r44/p1win/<d>.json, results/r44_p1win.json"""
import os, sys
ROOT = sys.argv[1]; PARTSW = sys.argv[2].split(','); NWW = int(sys.argv[3]) if len(sys.argv) > 3 else 4
_s = open(ROOT + '/code/r38_p1n.py').read(); _s = _s[:_s.index("ALLN = [")]
sys.argv = ['r38_p1n', ROOT, 'none', str(NWW)]; exec(_s)   # r38_p1 + r38_p1n 의 함수: frames, fs, gpa_profile, md_atoms, _npz_q, OUT1, RAWN, sf_med, MINSH1, log, jdump
OUTW = ROOT + '/derived/r44/p1win'; os.makedirs(OUTW, exist_ok=True)
WINS = (('w050_100', 50, 100), ('w300_400', 300, 400), ('w400_500', 400, 500)); QW = 0.7; NWF = 50; MINRUNW = 2


def runs320(d):
    """[(절대 시작 프레임, 프레임 배열 (k, L, 3, 3))] 반복마다 (뒤 절반과 50–100 ns 조각), Qf, 반복 길이."""
    f = f'{RAWN}/{d}.npz'
    if os.path.exists(f):
        z = np.load(f); smd, Qf = _npz_q(z); reps = [int(r) for r in z['reps320']]; lens = [int(n) for n in z['nframes320']]
        segs = [[(n // 2, z[f'h320_{r}'])] + ([(50, z[f'a320_{r}'])] if f'a320_{r}' in z.files else []) for r, n in zip(reps, lens)]
        return segs, Qf, lens
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode()); L = len(smd)
    reps = [g['320'][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g['320'].keys(), key=int)]; h.close()
    ca0 = reps[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
    return [[(0, X)] for X in reps], Qf, [len(X) for X in reps]


def window_frames(segs, Qf, n, lo, hi):
    """한 반복의 [lo, hi) 창 프레임(그 창을 온전히 담는 조각에서)."""
    if hi > n: return None
    for s0, X in segs:
        if s0 <= lo and s0 + len(X) >= hi: return X[lo - s0:hi - s0]
    return None


def _win_worker(d):
    fo = f'{OUTW}/{d}.json'
    if os.path.exists(fo): return d, 'skip'
    try:
        mine = [(k, v) for k, v in AL.items() if d in k]
        if not mine: json.dump(dict(dom=d, skip='쌍 없음'), open(fo, 'w')); return d, 'none'
        segs, Qf, lens = runs320(d); cells = {}
        for w, lo, hi in WINS:
            rng = np.random.default_rng(zlib.crc32(f'r44_p1win|{d}|{w}|{QW:.2f}'.encode())); runs, lab = [], []
            for r, (sg, n) in enumerate(zip(segs, lens)):
                X = window_frames(sg, Qf, n, lo, hi)
                if X is None: continue
                idx = np.where(Qf(X) >= QW)[0]
                if len(idx) >= NWF: runs.append(X[np.sort(rng.choice(idx, NWF, replace=False))]); lab += [r] * NWF
            if len(runs) >= MINRUNW:
                E = frames(np.concatenate(runs)); m = np.asarray(fs.mean_point(E)); cells[w] = np.asarray(fs.align_to_reference(E, torch.as_tensor(m)), dtype=np.float64)
        res = dict(dom=d, lens=lens, cells=sorted(cells), pairs={})
        for (a, b), al in mine:
            side = 0 if d == a else 1; e0 = np.array(sorted({x[side] for x in al['best']}), int)
            if len(e0) < MINSH1: continue
            res['pairs'][f'{a}__{b}'] = dict(n_shared=int(len(e0)), gpa={w: float(gpa_profile(A, e0).sum()) for w, A in cells.items()})
        json.dump(res, open(fo, 'w')); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


def part_spread():
    global AL
    from multiprocessing import Pool
    AL = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(OUT1 + '/align.jsonl')}
    doms = sorted({x for p in AL for x in p}); cnt = {}
    with Pool(NWW) as p:
        for k, (d, s_) in enumerate(p.imap_unordered(_win_worker, doms)):
            cnt[s_[:10]] = cnt.get(s_[:10], 0) + 1
            if s_ not in ('skip', 'ok') or k % 50 == 0: log('p1win spread', k, '/', len(doms), d, s_)
    log('R44_P1WIN_SPREAD', cnt); log('R44_P1WIN_SPREAD_DONE')


def part_stats():
    pj = json.load(open(OUT1 + '/pairs.json')); SP = {}
    for d in pj['domains']:
        f = f'{OUTW}/{d}.json'
        if os.path.exists(f):
            z = json.load(open(f))
            if 'pairs' in z: SP[d] = z
    def lr(r, w):
        key = f"{r['a']}__{r['b']}"; H_, R_ = SP.get(r['hot']), SP.get(r['ref'])
        if H_ is None or R_ is None or key not in H_['pairs'] or key not in R_['pairs']: return None
        a, b = H_['pairs'][key]['gpa'].get(w), R_['pairs'][key]['gpa'].get(w)
        return float(np.log(a / b)) if a and b else None
    res = dict(written=time.strftime('%F %T'))
    for nm, sel in (('all', lambda r: True), ('indep', lambda r: r['set'] == 'indep')):
        HT = [r for r in pj['HT'] if sel(r)]; MM = [r for r in pj['MM'] if sel(r)]; o = {}
        D = {}
        for w, _, _ in WINS:
            Ts, Cs = sf_med(HT, lambda r: lr(r, w)), sf_med(MM, lambda r: lr(r, w)); D[w] = {h: Ts[h] - Cs[h] for h in sorted(set(Ts) & set(Cs))}
        common = sorted(set.intersection(*[set(v) for v in D.values()])) if all(D.values()) else []
        for tag, keys in (('each', None), ('common', common)):
            for w, _, _ in WINS:
                dv = np.array([D[w][h] for h in (D[w] if keys is None else keys) if h in D[w]]); rng = np.random.default_rng(20261096)
                bd = [np.median(dv[rng.integers(0, len(dv), len(dv))]) for _ in range(2000)] if len(dv) else []
                o[f'{tag}|{w}'] = dict(n_sf=int(len(dv)), U2=float(np.median(dv)) if len(dv) else None, hi975=float(np.percentile(bd, 97.5)) if bd else None, lo025=float(np.percentile(bd, 2.5)) if bd else None)
        res[nm] = o
    a = res['all']; late = [a['each|w300_400']['U2'], a['each|w400_500']['U2']]
    res['predictions'] = dict(pa1_late_both_negative=bool(all(v is not None and v < 0 for v in late)),
                              pa2_300_400_upper_negative=bool(a['each|w300_400']['hi975'] is not None and a['each|w300_400']['hi975'] < 0))
    res['sentence'] = ('늦은 절대 창(300–400·400–500 ns)에서도 호열 쪽이 덜 흔들린다 — 창에 매인 결과가 아니다' if res['predictions']['pa1_late_both_negative'] else
                       '절대 창에 따라 U2 가 달라진다 — 뒤 절반 창의 결과는 창 조건에 매인다')
    res['rule'] = '등록 r44_p1win: 서술 + 예측 둘(PA1 늦은 두 창의 U2 점추정 < 0, PA2 300–400 ns 창의 97.5% 상한 < 0). 50–100 ns 는 같은 Q 0.7 로 다시 낸 비교값.'
    jdump(res, RES + '/r44_p1win.json', bk=True); log('R44_P1WIN', res['predictions'], res['sentence']); log('R44_P1WIN_STATS_DONE')


if __name__ == '__main__':
    for part in PARTSW:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('R44_P1WIN_DONE', ','.join(PARTSW))
