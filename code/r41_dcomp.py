"""(round41 g2, 2026-10-21 UTC; Drive JSON round41_2026-10-21) D 압축 기전: E7 100 의 같은 실행·프레임 쌍에서 네 수준의 거리.
  (a) 프레임마다 320 K 기준 평균 점에 회전 하나로 맞춘 좌표의 유클리드(저장 A2000 과 같은 양), (b) d_fix (쌍마다 회전, 항등 대응, SVD 한 번),
  (c1) 항등 대응 + 변마다 min(짝, 빈칸 둘) = ‖x‖² + ‖y‖² − 2 Σ_k max(0, ⟨x_k, U y_k⟩) (U 는 (b) 에서 출발해 남긴 변으로 Procrustes 되풀이, 최대 10), (c2) D (저장 칸 값).
부분: cpu (a·b·c1 칸 값; 작업자 여럿, GPU 안 씀) → gpu (실행 쌍마다 앞 200 프레임 쌍의 줄인 계산 다시: 쌍마다 c2 와 c1, 빈칸 변·i ≠ j 짝 비율을 끝 10 변/중심으로;
  첫 3 도메인은 2,000 전부로 저장 칸 C·B 와 대조) → sum
프레임 쌍 = r38_e7_30.one() 과 같은 rng(seeded('r38_e7_30', d))·같은 upairs 호출 순서. (a) 칸 값이 저장 A2000 과 상대 1e-9 안이어야 한다(재현 점검).
사용: python r41_dcomp.py <ROOT> <부분들 쉼표> [조각 수=1] [조각 번호=0]   (cpu 는 조각마다 프로세스 하나: 도메인 i 중 i % 조각 수 == 번호)
출력: derived/r41/dcomp/<dom>.json, derived/r41/dcomp_gpu/<dom>.json, results/r41_dcomp.json"""
import os, sys
ROOT = sys.argv[1]; PARTS41 = sys.argv[2].split(','); NW41 = int(sys.argv[3]) if len(sys.argv) > 3 else 1; SH41 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
if 'gpu' not in PARTS41: os.environ['CUDA_VISIBLE_DEVICES'] = ''
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]; assert _src.count("os.environ[k_] = '4'") == 1
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src.replace("os.environ[k_] = '4'", "os.environ[k_] = '2'"))
_e = open(ROOT + '/code/r38_e7_30.py').read(); exec(_e[_e.index('def ref_of'):_e.index('def one')])   # ref_of (그대로)
from framespace.geometry import procrustes
OUT41 = ROOT + '/derived/r41/dcomp'; OUTG = ROOT + '/derived/r41/dcomp_gpu'; os.makedirs(OUT41, exist_ok=True); os.makedirs(OUTG, exist_ok=True)
FO = RES + '/r41_dcomp.json'; NSUB = 200; CUTE = 10; LV3 = ('a', 'b', 'c1')


def doms():
    pk = json.load(open(OUT + '/e7_100/pick.json')); e30 = json.load(open(OUT + '/e7_30.json'))['e7_30']
    L = [(d, f'{OUT}/e7_30/{d}.json') for d in e30] + [(d, f'{OUT}/e7_100/{d}.json') for d in pk['e7_5'] + pk['new65']]
    return [(d, f) for d, f in L if os.path.exists(f) and any(k.startswith('MD') for k in json.load(open(f)).get('cells', {}))]


def lv_pairs(x, y, Aa=None, Ab=None):
    """x, y (P, n, 3, 3) torch float64 (같은 n). 쌍마다 (a), (b), (c1) 거리, c1 의 포화 변 비율, 점검 값."""
    da = np.sqrt(((Aa - Ab) ** 2).sum(1)) if Aa is not None else None
    nx = (x ** 2).sum((1, 2, 3)); ny = (y ** 2).sum((1, 2, 3))
    H = torch.einsum('pkcb,pkab->pca', y, x); sv = torch.linalg.svdvals(H); Sb = sv[:, 0] + sv[:, 1] + torch.sign(torch.det(H)) * sv[:, 2]
    db = (nx + ny - 2 * Sb).clamp_min(0).sqrt()
    w = torch.ones(x.shape[:2], dtype=x.dtype, device=x.device); best = torch.full((len(x),), -1e30, dtype=x.dtype, device=x.device); sat = torch.zeros(len(x), dtype=x.dtype, device=x.device); S0 = None
    for it in range(11):
        U = procrustes(torch.einsum('pk,pkxz,pkyz->pxy', w, y, x)); K = torch.einsum('pkxz,pxy,pkyz->pk', x, U, y)
        if it == 0: S0 = K.sum(1)
        S = K.clamp_min(0).sum(1); w2 = (K > 0).to(x.dtype); imp = S > best
        best = torch.where(imp, S, best); sat = torch.where(imp, 1 - w2.mean(1), sat)
        if torch.equal(w2, w): break
        w = w2
    dc = (nx + ny - 2 * best).clamp_min(0).sqrt()
    chk = float(((S0 - Sb).abs() / Sb.abs().clamp_min(1e-12)).max())     # 첫 걸음(모든 변) = d_fix 의 SVD 최적이어야 함
    return da, db.cpu().numpy(), dc.cpu().numpy(), sat.cpu().numpy(), chk


def replay(d, fstored, fn):
    """one() 과 같은 호출 순서로 실행·프레임 쌍을 다시 뽑고, 실행 쌍마다 fn(Ta, Tb, rp, ia, ib, E, A) 를 부른다. 칸 정의도 one() 과 같다."""
    rng = seeded('r38_e7_30', d); _runs.clear(); _md.clear(); ref = ref_of(d); per = {}
    for W_ in ('full', 'half'):
        R = {T: dom_runs(d, T, W_, TT)[0] for T in TT}; L = dom_runs(d, '320', W_, TT)[1]; ne = L - 1
        Ts = [T for T in TT if len(R[T]) >= 2]
        if '320' not in Ts: continue
        E = {T: coords(np.concatenate([X for _, X in R[T]])) for T in Ts}
        A = {T: fs.align_to_reference(frames(np.concatenate([X for _, X in R[T]])), ref).numpy().reshape(len(R[T]) * N, -1) for T in Ts}
        def go(Ta, Tb, rr):
            return {rp: fn(Ta, Tb, rp, *upairs(rp[0], rp[1], rng), E, A) for rp in rr}
        pm = {}; n320 = len(R['320']); pm[('320', '320')] = go('320', '320', list(itertools.combinations(range(n320), 2)))
        for T in Ts[1:]:
            nT = len(R[T]); pm[(T, T)] = go(T, T, list(itertools.combinations(range(nT), 2))); pm[('320', T)] = go('320', T, list(itertools.product(range(n320), range(nT))))
            per[f'MD{T}-MD320|{W_}'] = dict(cr=list(pm[('320', T)].values()), fp=list(pm[('320', '320')].values()), fq=list(pm[(T, T)].values()), ne=ne)
        if n320 == 5:
            f = lambda a, b: pm[('320', '320')][(min(a, b), max(a, b))]
            for i, g2 in enumerate(itertools.combinations(range(5), 2)):
                g3 = [k for k in range(5) if k not in g2]
                per[f'null320|{W_}|{i}'] = dict(cr=[f(a, b) for a in g2 for b in g3], fp=[f(*g2)], fq=[f(a, b) for a, b in itertools.combinations(g3, 2)], ne=ne)
    return per


def cpu_one(job):
    d, fstored = job; fo = f'{OUT41}/{d}.json'
    if os.path.exists(fo): return d, '있음'
    torch.set_num_threads(2); t0 = time.time(); st = json.load(open(fstored)); CH = dict(Sb_iter=0.0, viol_ab=0, viol_bc1=0, n_pairs=0)
    def fn(Ta, Tb, rp, ia, ib, E, A):
        da, db, dc, sat, chk = lv_pairs(E[Ta][ia], E[Tb][ib], A[Ta][ia], A[Tb][ib])
        CH['Sb_iter'] = max(CH['Sb_iter'], chk); CH['viol_ab'] += int((db > da + 1e-9 * (1 + da)).sum()); CH['viol_bc1'] += int((dc > db + 1e-9 * (1 + db)).sum()); CH['n_pairs'] += len(da)
        return dict(a=float(da.mean()), b=float(db.mean()), c1=float(dc.mean()), a2=float((da ** 2).mean()), b2=float((db ** 2).mean()), c12=float((dc ** 2).mean()), sat=float(sat.mean()))
    try:
        per = replay(d, fstored, fn)
    except Exception as e:
        traceback.print_exc(); return d, 'ERR ' + repr(e)[:200]
    cells = {}; worst = 0.0
    for k, p in per.items():
        if k not in st['cells']: continue
        m = lambda xs, v: float(np.mean([x[v] for x in xs])); c = {}
        for v in LV3:
            c[v] = estats(m(p['cr'], v), m(p['fp'], v), m(p['fq'], v), p['ne']); c[v]['dmu2'] = float((m(p['cr'], v + '2' if v != 'c1' else 'c12') - (m(p['fp'], v + '2' if v != 'c1' else 'c12') + m(p['fq'], v + '2' if v != 'c1' else 'c12')) / 2) / p['ne'])
        s = st['cells'][k]; c['c2'] = {q: s['D'][q] for q in ('C', 'B_P', 'B_Q', 'D2', 'wbar', 'ratio')}; c['c2']['dmu2'] = s['dmu2_D']
        c['sat_cross'] = m(p['cr'], 'sat'); c['sat_within'] = float(np.mean([x['sat'] for x in p['fp'] + p['fq']]))
        for q in ('C', 'B_P', 'B_Q'): worst = max(worst, abs(c['a'][q] - s['A2000'][q]) / max(abs(s['A2000'][q]), 1e-12))
        cells[k] = c
    r = dict(dom=d, cells=cells, check=dict(a_vs_A2000_maxrel=worst, **CH), sec=round(time.time() - t0, 1))
    jdump(r, fo, bk=False); return d, f"칸 {len(cells)} a 재현 {worst:.1e} 순서 위반 {CH['viol_ab']}/{CH['viol_bc1']} {r['sec']} s"


def part_cpu():
    J = doms(); my = [x for i, x in enumerate(J) if i % NW41 == SH41]; log('dcomp cpu', len(J), '도메인 중 조각', SH41, '/', NW41, '=', len(my))
    for x in my:
        d, msg = cpu_one(x); log('dcomp', d, msg)
    log(f'G2_CPU_DONE_{SH41}')


# ---------- gpu ----------
def red_pairs_rec(ms, ia, ib, k0=K0, iters=10, membudget=3e9):
    """r38_s0.red_pairs 와 같은 계산·같은 연산 순서(want_gap 없음). 쌍마다 마지막 풀이의 짝에서 끝(양 끝 CUTE 변)/중심 별 빈칸 변 수·짝 수·i ≠ j 짝 수를 남긴다."""
    EA, EB = ms.EA, ms.EB; n, m = EA.shape[1], EB.shape[1]
    ia = torch.as_tensor(np.asarray(ia), device=EA.device); ib = torch.as_tensor(np.asarray(ib), device=EA.device)
    S = EA.new_zeros(len(ia)); ST_ = torch.zeros(len(ia), 8, dtype=torch.float64, device=EA.device); kused = np.zeros(len(ia), int)
    tA = (torch.arange(n, device=EA.device) < CUTE) | (torch.arange(n, device=EA.device) >= n - CUTE); tB = (torch.arange(m, device=EA.device) < CUTE) | (torch.arange(m, device=EA.device) >= m - CUTE)
    todo = torch.arange(len(ia), device=EA.device); kk = k0
    while len(todo):
        lo, w = band_windows(ms.center, m, kk); step = max(1, int(membudget / (n * w * 9 * 8 * 4)))
        edge_all = torch.zeros(len(todo), dtype=torch.bool, device=EA.device)
        for c in range(0, len(todo), step):
            sel = todo[c:c + step]; a, b = EA[ia[sel]], EB[ib[sel]]
            sc, _, M, edge = align_band(a, b, lo, w, ms.starts(ia[sel], ib[sel], n_local=1), iters=iters)
            S[sel] = sc; edge_all[c:c + step] = edge
            uA, uB = ~M.any(2), ~M.any(1); p_, i_, j_ = M.nonzero(as_tuple=True); P_ = len(sel); sh = (i_ != j_).double(); te = tA[i_].double()
            ST_[sel] = torch.stack([(uA & tA).sum(1).double(), (uA & ~tA).sum(1).double(), (uB & tB).sum(1).double(), (uB & ~tB).sum(1).double(),
                                    torch.bincount(p_, te, P_).double(), torch.bincount(p_, 1 - te, P_).double(), torch.bincount(p_, sh * te, P_).double(), torch.bincount(p_, sh * (1 - te), P_).double()], 1)
        kused[todo.cpu().numpy()] = kk
        if w >= m: break
        todo = todo[edge_all]; kk *= 2
    d2 = (sq_norm(EA)[ia] + sq_norm(EB)[ib] - 2 * S).clamp_min(0)
    return d2.sqrt().cpu().numpy(), ST_.cpu().numpy(), dict(ne_A=n, ne_B=m, nt_A=int(tA.sum()), nt_B=int(tB.sum()))


def gpu_one(d, fstored, full):
    fo = f'{OUTG}/{d}.json'
    if os.path.exists(fo): return '있음'
    t0 = time.time(); ms_cache = {}; st = json.load(open(fstored))
    def fn(Ta, Tb, rp, ia, ib, E, A):
        k = (Ta, Tb, id(E))     # 창마다 E 가 다르다; 캐시가 E 를 붙잡아 id 재사용을 막는다
        if k not in ms_cache: ms_cache.clear(); ms_cache[k] = (MeanStart(E[Ta], E[Tb]), E)
        if not full: ia, ib = ia[:NSUB], ib[:NSUB]
        dv, stt, nm = red_pairs_rec(ms_cache[k][0], ia, ib); _, db, dc, sat, _ = lv_pairs(E[Ta][ia], E[Tb][ib])
        return dict(d=float(dv.mean()), d2=float((dv ** 2).mean()), c2_gt_c1=int((dv > dc + 1e-9 * (1 + dc)).sum()), c2_lt_c1=int((dv < dc - 1e-9 * (1 + dc)).sum()), n=len(dv),
                    st=stt.sum(0).tolist(), c1=float(dc.mean()), b=float(db.mean()), **nm)
    per = replay(d, fstored, fn); cells = {}; chk = 0.0
    for k, p in per.items():
        if k not in st['cells']: continue
        S_ = lambda xs: np.sum([x['st'] for x in xs], 0); n_ = lambda xs: sum(x['n'] for x in xs); nm = p['cr'][0]
        o = {}
        for part, xs in (('cross', p['cr']), ('within', p['fp'] + p['fq'])):
            s, nn = S_(xs), n_(xs)
            o[part] = dict(n_pairs=nn, gap_frac_term=float((s[0] + s[2]) / (nn * (nm['nt_A'] + nm['nt_B']))), gap_frac_core=float((s[1] + s[3]) / (nn * (nm['ne_A'] + nm['ne_B'] - nm['nt_A'] - nm['nt_B']))),
                           shift_frac_term=float(s[6] / max(s[4], 1)), shift_frac_core=float(s[7] / max(s[5], 1)), frac_c2_gt_c1=float(sum(x['c2_gt_c1'] for x in xs) / nn), frac_c2_lt_c1=float(sum(x['c2_lt_c1'] for x in xs) / nn),
                           mean_c2=float(np.mean([x['d'] for x in xs])), mean_c1=float(np.mean([x['c1'] for x in xs])))
        if full:
            m = lambda xs: float(np.mean([x['d'] for x in xs])); s = st['cells'][k]['D']
            for q, v in (('C', m(p['cr'])), ('B_P', m(p['fp'])), ('B_Q', m(p['fq']))): chk = max(chk, abs(v - s[q]) / max(abs(s[q]), 1e-12))
        cells[k] = o
    r = dict(dom=d, full=full, cells=cells, check_full_vs_stored_maxrel=chk if full else None, sec=round(time.time() - t0, 1)); jdump(r, fo, bk=False)
    return f"칸 {len(cells)} {'전부 대조 ' + format(chk, '.1e') if full else '부분 200'} {r['sec']} s"


def part_gpu():
    J = doms(); log('dcomp gpu', len(J), DEV)
    for i, (d, f) in enumerate(J):
        try: log('dcomp gpu', d, gpu_one(d, f, i < 3))
        except Exception as e: traceback.print_exc(); log('ERR dcomp gpu', d, repr(e)[:200])
        if DEV == 'cuda': torch.cuda.empty_cache()
    log('G2_GPU_DONE')


# ---------- sum ----------
def cath_sf():
    h = {}
    for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
        if l.startswith('#') or not l.strip(): continue
        f = l.split(); h[f[0]] = '.'.join(f[1:5])
    return h


def part_sum():
    sf = cath_sf(); J = [json.load(open(f)) for f in sorted(glob.glob(OUT41 + '/*.json'))]; G = {j['dom']: j for j in (json.load(open(f)) for f in sorted(glob.glob(OUTG + '/*.json')))}
    rows = []
    for j in J:
        for k, c in j['cells'].items():
            if not k.startswith('MD'): continue
            T, W = k[2:5], k.split('|')[1]; r = dict(dom=j['dom'], sf=sf.get(j['dom'], j['dom']), T=T, W=W)
            for q in ('D2', 'wbar', 'ratio', 'dmu2'): r.update({f'{q}_{v}': c[v][q] for v in ('a', 'b', 'c1', 'c2')})
            r['sat_cross'] = c['sat_cross']; g = G.get(j['dom'], {}).get('cells', {}).get(k)
            if g: r.update({f'g_{x}': g['cross'][x] for x in ('gap_frac_term', 'gap_frac_core', 'shift_frac_term', 'shift_frac_core', 'frac_c2_gt_c1')})
            rows.append(r)
    rng = np.random.default_rng(20261085); sfl = sorted({r['sf'] for r in rows}); si = {s: i for i, s in enumerate(sfl)}
    BI = [rng.integers(0, len(sfl), len(sfl)) for _ in range(2000)]
    def steps(sub, q):
        ok = [r for r in sub if all(np.isfinite(r[f'{q}_{v}']) and r[f'{q}_{v}'] > 0 for v in ('a', 'b', 'c1', 'c2'))]
        if len(ok) < 5: return dict(n=len(ok))
        Ls = np.array([[np.log(r[f'{q}_b'] / r[f'{q}_a']), np.log(r[f'{q}_c1'] / r[f'{q}_b']), np.log(r[f'{q}_c2'] / r[f'{q}_c1'])] for r in ok]); g = np.array([si[r['sf']] for r in ok])
        def summ(Lx):
            med = np.median(Lx, 0); return np.concatenate([med, med / med.sum() if abs(med.sum()) > 1e-12 else [np.nan] * 3, [np.exp(med.sum())]])
        pt = summ(Ls); mem = [np.where(g == s)[0] for s in range(len(sfl))]; bs = []
        for b in BI:
            ix = np.concatenate([mem[s] for s in b]);
            if len(ix): bs.append(summ(Ls[ix]))
        bs = np.array(bs); lo, hi = np.nanpercentile(bs, 2.5, 0), np.nanpercentile(bs, 97.5, 0)
        nm = ['ln_a_b', 'ln_b_c1', 'ln_c1_c2', 'share_a_b', 'share_b_c1', 'share_c1_c2', 'ratio_c2_over_a']
        return dict(n=len(ok), **{x: dict(est=float(pt[i]), ci95=[float(lo[i]), float(hi[i])]) for i, x in enumerate(nm)},
                    median_ratio_vs_a={v: float(np.median([r[f'{q}_{v}'] / r[f'{q}_a'] for r in ok])) for v in ('b', 'c1', 'c2')})
    out = dict(written=time.strftime('%F %T'), n_dom=len(J), n_cells=len(rows), checks=dict(a_vs_A2000_maxrel=max(j['check']['a_vs_A2000_maxrel'] for j in J), Sb_iter_maxrel=max(j['check']['Sb_iter'] for j in J),
               viol_ab=sum(j['check']['viol_ab'] for j in J), viol_bc1=sum(j['check']['viol_bc1'] for j in J), n_pairs=sum(j['check']['n_pairs'] for j in J),
               gpu_full_vs_stored_maxrel=max([g['check_full_vs_stored_maxrel'] for g in G.values() if g.get('full')], default=None), n_gpu=len(G)))
    for q in ('D2', 'wbar', 'ratio', 'dmu2'):
        out[q] = dict(all=steps(rows, q), by_T={T: steps([r for r in rows if r['T'] == T], q) for T in ('348', '379', '413', '450')}, by_W={W: steps([r for r in rows if r['W'] == W], q) for W in ('full', 'half')})
    rd = {}
    for q in ('D2', 'ratio'):
        a = out[q]['all']
        if 'share_a_b' in a:
            sh = {'a→b': a['share_a_b']['est'], 'b→c1': a['share_b_c1']['est'], 'c1→c2': a['share_c1_c2']['est']}; top = max(sh, key=sh.get)
            rd[q] = dict(shares=sh, reading={'a→b': '쌍별 회전이 집단 운동을 흡수(D 는 내부 좌표형)', 'b→c1': '120° 넘는 변의 포화(경첩)', 'c1→c2': '등록 미끄러짐(끝 재정렬·주기 구조)'}[top] if sh[top] > 0.5 else '여럿에 나뉨')
    out['reading'] = rd; out['sat_cross_median'] = float(np.median([r['sat_cross'] for r in rows]))
    if G:
        gx = [r for r in rows if 'g_gap_frac_term' in r]
        out['gpu'] = {x: dict(median=float(np.median([r[f'g_{x}'] for r in gx])), q25=float(np.percentile([r[f'g_{x}'] for r in gx], 25)), q75=float(np.percentile([r[f'g_{x}'] for r in gx], 75)))
                      for x in ('gap_frac_term', 'gap_frac_core', 'shift_frac_term', 'shift_frac_core', 'frac_c2_gt_c1')} | dict(n_cells=len(gx))
    out['rule'] = 'round41 g2 (Drive JSON round41_2026-10-21): 칸마다 단계 로그 비, 단계 중앙값과 그 합에 대한 몫, 상과 부트스트랩 2,000 (seed 20261085). 몫 > 0.5 인 단계로 읽음.'
    jdump(out, FO); log('G2_SUM', out['checks'], out['reading'])


if __name__ == '__main__':
    for p in PARTS41:
        log('부분 시작', p); dict(cpu=part_cpu, gpu=part_gpu, sum=part_sum)[p]()
    log('G2_ALL_DONE')
