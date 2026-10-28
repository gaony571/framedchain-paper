"""(round44, 등록 = Drive results/prereg_2026-10-28_addendum6.json "r44_redsame", 2026-10-28 11:50:14 UTC) 같은 서열 온도 칸에서 줄인 계산의 과대평가가 Δ²_D 의 압축을 만드는가.
E7 100 도메인(derived/r38/s0/e7_100/·e7_30/ 의 도메인 파일) 가운데 "MD379-MD320|full" 칸이 있는 것을 정렬해 default_rng(20261102).permutation 순서로 16 개.
도메인마다 r38_e7_30.one 과 같은 실행·프레임(dom_runs, N = 100, Q ≥ 0.75, 전체 창)과 같은 MeanStart(E320, E379)·(E320, E320)·(E379, E379).
쌍 종류 셋(cross 320–379, floor320, floor379) × 20: 실행 무작위(바닥은 서로 다른 두 실행), 프레임 무작위, seed = seeded('r44_redsame', d, kind, t).
err = d_red − d_ref, d_ref = √(‖x‖² + ‖y‖² − 2 max(S_lb, S_red)), S_lb = fs.bnb(rel_tol = REL) — r38_s0.part_bias 와 같은 정의.
b = 2·ē_cross − ē_floor320 − ē_floor379 (Δ²_D 에 주는 치우침), g = Δ²(A2000) − Δ²(D) (저장된 칸), 몫 R = −Σb / Σg, 도메인 재추출 2,000 (seed 20261103).
사용: python r44_redsame.py <ROOT> run,sum   출력: derived/r44/redsame/<d>.jsonl, results/r44_redsame.json"""
import os, sys, json, glob, time, traceback
ROOT = sys.argv[1]; PARTS_RS = sys.argv[2].split(',')   # r38_s0 가 PARTS 를 덮어쓰므로 다른 이름
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)   # dom_runs, coords, MeanStart, red_pairs, sq_norm, seeded, fs, N, REL, DEV, TT, _runs, _md, log
OUTR = ROOT + '/derived/r44/redsame'; os.makedirs(OUTR, exist_ok=True); NPK = 20; KINDS = ('cross', 'floor320', 'floor379'); CELL = 'MD379-MD320|full'


def cells():
    out = {}
    for f in sorted(glob.glob(ROOT + '/derived/r38/s0/e7_100/*.json') + glob.glob(ROOT + '/derived/r38/s0/e7_30/*.json')):
        if os.path.basename(f) == 'pick.json': continue
        z = json.load(open(f))
        if isinstance(z, dict) and CELL in z.get('cells', {}): out[z['dom']] = z['cells'][CELL]
    return out


def pick():
    C = cells(); doms = sorted(C); order = np.random.default_rng(20261102).permutation(len(doms))
    return [doms[int(o)] for o in order[:16]], C


def part_run():
    sel, _ = pick(); log('redsame 도메인', len(sel), sel)
    for d in sel:
        fo = f'{OUTR}/{d}.jsonl'; done = set()
        if os.path.exists(fo):
            for l in open(fo):
                try: j = json.loads(l); done.add((j['kind'], j['t']))
                except Exception: pass
        if len(done) >= NPK * len(KINDS): continue
        t0 = time.time(); _runs.clear(); _md.clear()
        R = {T: dom_runs(d, T, 'full', TT)[0] for T in ('320', '379')}
        E = {T: coords(np.concatenate([X for _, X in R[T]])) for T in R}; nr = {T: len(R[T]) for T in R}
        MS = dict(cross=MeanStart(E['320'], E['379']), floor320=MeanStart(E['320'], E['320']), floor379=MeanStart(E['379'], E['379']))
        for t in range(NPK):
            for kind in KINDS:
                if (kind, t) in done: continue
                rng = seeded('r44_redsame', d, kind, t)
                if kind == 'cross':
                    i = int(rng.integers(nr['320'])) * N + int(rng.integers(N)); j = int(rng.integers(nr['379'])) * N + int(rng.integers(N))
                else:
                    T = kind[-3:]; ra, rb = rng.choice(nr[T], 2, replace=False); i = int(ra) * N + int(rng.integers(N)); j = int(rb) * N + int(rng.integers(N))
                ms = MS[kind]; d_red, _, S_red, ku, _ = red_pairs(ms, [i], [j])
                x, y = ms.EA[i], ms.EB[j]; nx, ny = float(sq_norm(x)), float(sq_norm(y))
                cert = fs.bnb(x, y, rel_tol=REL, device=DEV); S_ref = max(cert['S_lb'], float(S_red[0])); d_ref = float(np.sqrt(max(nx + ny - 2 * S_ref, 0)))
                rec = dict(dom=d, kind=kind, t=t, i=i, j=j, d_red=float(d_red[0]), d_ref=d_ref, err=float(d_red[0]) - d_ref, S_red=float(S_red[0]), S_lb=cert['S_lb'],
                           S_ub=cert['S_ub'], certified=bool(cert['certified']), k_used=int(ku[0]), bnb_sec=cert['sec'])
                with open(fo, 'a') as fh: fh.write(json.dumps(rec) + '\n')
        log('redsame', d, round(time.time() - t0), 's')
    log('R44_REDSAME_RUN_DONE')


def part_sum():
    sel, C = pick(); rows = []
    for d in sel:
        f = f'{OUTR}/{d}.jsonl'
        if not os.path.exists(f): continue
        J = [json.loads(l) for l in open(f) if l.strip()]; e = {k: [j['err'] for j in J if j['kind'] == k] for k in KINDS}
        if any(len(v) < NPK for v in e.values()): continue
        m = {k: float(np.mean(v)) for k, v in e.items()}; b = 2 * m['cross'] - m['floor320'] - m['floor379']
        g = float(C[d]['A2000']['D2'] - C[d]['D']['D2'])
        rows.append(dict(dom=d, e=m, b=b, g=g, D2_D=float(C[d]['D']['D2']), D2_A2000=float(C[d]['A2000']['D2']), n_uncert=int(sum(not j['certified'] for j in J)),
                         d_ref_mean=float(np.mean([j['d_ref'] for j in J]))))
    b = np.array([r['b'] for r in rows]); g = np.array([r['g'] for r in rows]); rng = np.random.default_rng(20261103)
    Rhat = float(-b.sum() / g.sum()); bs = []
    for _ in range(2000):
        ix = rng.integers(0, len(rows), len(rows)); bs.append(-b[ix].sum() / g[ix].sum())
    em = {k: float(np.mean([r['e'][k] for r in rows])) for k in KINDS}
    out = dict(written=time.strftime('%F %T UTC', time.gmtime()), n_dom=len(rows), pairs_per_kind=NPK, cell=CELL,
               R=dict(est=Rhat, lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5))),
               per_domain_share_median=float(np.median(-b / g)), mean_err=em, mean_b=float(b.mean()), mean_g=float(g.mean()),
               compression_ratio_median=float(np.median([r['D2_D'] / r['D2_A2000'] for r in rows])), n_uncertified=int(sum(r['n_uncert'] for r in rows)), rows=rows)
    out['predictions'] = dict(pr1_upper_lt_01=bool(out['R']['hi'] < 0.1), pr2_cross_ge_floor=bool(em['cross'] >= (em['floor320'] + em['floor379']) / 2))
    out['rule'] = '등록 r44_redsame (prereg_2026-10-28_addendum6.json): PR1 몫 R = −Σb/Σg 의 95% 상한 < 0.1; PR2 ē_교차 ≥ (ē_바닥320 + ē_바닥379)/2.'
    json.dump(out, open(ROOT + '/results/r44_redsame.json', 'w'), ensure_ascii=False, indent=1)
    log('R44_REDSAME', {k: out[k] for k in ('n_dom', 'R', 'per_domain_share_median', 'mean_err', 'mean_b', 'mean_g', 'compression_ratio_median', 'n_uncertified', 'predictions')})
    log('R44_REDSAME_SUM_DONE')


if __name__ == '__main__':
    for part in PARTS_RS:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
