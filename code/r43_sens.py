"""(round43, 2026-10-22 UTC; Drive JSON round43_2026-10-22 08:54:19) 사후 민감도 셋 — 원래 판정은 바꾸지 않는다. CPU 만.
  s4  P2 2차 1 의 구간: 도메인별 above_null 비율의 상과 군집 부트스트랩 (2,000, seed 20261090) — 주 판 D, 맞춘 좌표, 있는 민감도 D 판.
  s1  그림 5 칸의 길이 맞춤: 칸마다 전체 판과 n̄ ∈ [80, 150] 판의 상과 중앙값 [상과 부트스트랩 2,000, seed 20261091] — 비, Δ²_D/√n̄, (같은 서열 칸) Δμ²₂_D/변.
  s3  C-2 의 ANM 민감도: 변 퍼짐 예측 = ANM 10·13·15 Å (r38_c.anm_edge_spread 원문; 15 Å 는 저장값 재현 점검) · GNM 7.3 Å (C_ii + C_jj − 2C_ij);
      r_ANM 은 r38_c.anm_pairs 와 같은 규칙, C-2 = r38_c.judge(control + r_ANM_판) (seed C 20261074, H-풀 20261077).
사용: python r43_sens.py <ROOT> <부분들 쉼표> [작업자=10]   출력: results/r43_sens.json, derived/r43/anm/<d>.npz"""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
import sys, json, glob, time, shutil, warnings
import numpy as np
from scipy.stats import kendalltau
warnings.filterwarnings('ignore', category=RuntimeWarning)
ROOT = sys.argv[1]; PARTS = sys.argv[2].split(','); NW = int(sys.argv[3]) if len(sys.argv) > 3 else 10
RES, DER = ROOT + '/results', ROOT + '/derived'; P3, CC = DER + '/r38/p3', DER + '/r38/c'; OUTA = DER + '/r43/anm'; os.makedirs(OUTA, exist_ok=True); FO = RES + '/r43_sens.json'
NB = 2000; MINSH = 30; ANM_CUT = 15.0
sys.path[:0] = [ROOT + '/framespace_pkg']
from framespace.frames import ca_trace
_p3 = open(ROOT + '/code/r38_p3.py').read()
exec(_p3[_p3.index('def ranks'):_p3.index('def analyse(')]); exec(_p3[_p3.index('def spearman'):_p3.index('def one_pair')]); exec(_p3[_p3.index('def verdict'):_p3.index('def part_stats')])
_c = open(ROOT + '/code/r38_c.py').read()
exec(_c[_c.index('def anm_edge_spread'):_c.index('def _anm_worker')]); exec(_c[_c.index('def analyse_sf'):_c.index('def merged')]); exec(_c[_c.index('def judge'):_c.index('def part_stats')])


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def jload():
    return json.load(open(FO)) if os.path.exists(FO) else dict(rule='round43 (Drive JSON round43_2026-10-22): 사후 민감도, 원래 판정 그대로')


def jdump(o):
    json.dump(o, open(FO + '.tmp', 'w'), ensure_ascii=False, indent=1, default=float); os.replace(FO + '.tmp', FO)
    b = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(FO, b + '/' + time.strftime('%H%M%S_') + os.path.basename(FO))


def cath_sf():
    h = {}
    for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
        if l.startswith('#') or not l.strip(): continue
        f = l.split(); h[f[0]] = '.'.join(f[1:5])
    return h


def pct(bs, a, b):
    return [float(np.percentile(bs, a)), float(np.percentile(bs, b))]


# ---------- s4 ----------
def part_s4():
    rng = np.random.default_rng(20261090); res = {}
    fs_ = sorted(set(glob.glob(RES + '/r38_p2_p2_*_D.json')) | {RES + '/r38_p2_p2_67.json'})
    for f in fs_:
        tag = os.path.basename(f)[len('r38_p2_'):-5]; rows = [r for r in json.load(open(f)).get('rows', []) if r.get('above_null') is not None]
        if not rows: continue
        sfl = sorted({r['sf'] for r in rows}); mem = {s: np.array([r['above_null'] for r in rows if r['sf'] == s], float) for s in sfl}; pt = float(np.mean([r['above_null'] for r in rows]))
        bs = [np.concatenate([mem[sfl[i]] for i in rng.integers(0, len(sfl), len(sfl))]).mean() for _ in range(NB)]
        c95 = pct(bs, 2.5, 97.5); res[tag] = dict(n=len(rows), n_sf=len(sfl), frac=pt, ci95=c95, ci90=pct(bs, 5, 95), pass_point=pt >= 0.70, ci95_contains_070=c95[0] <= 0.70 <= c95[1])
        log('s4', tag, res[tag])
    out = jload(); out['s4'] = dict(written=time.strftime('%F %T'), sets=res, rule='도메인 비율, 상과 복원 2,000 (seed 20261090). 등록 판정(점추정 ≥ 0.70)은 그대로.'); jdump(out); log('S4_DONE')


# ---------- s1 ----------
def cells():
    C = {}
    J = [json.load(open(f)) for f in sorted(glob.glob(DER + '/r39/D/*.json'))]
    for T, dT in (('348', 28), ('379', 59), ('413', 93), ('450', 130)):
        C[f'온도 +{dT} K (뒤 절반)'] = [dict(dom=j['dom'], n=float(j['ne']), ratio=c['D']['ratio'], d2n=c['D']['d2n'], dmu=c['dmu2_D']) for j in J for c in [j['cells'].get(f'MD{T}-MD320|half')] if c]
    for tag, sub in (('R3′ (50–100 ns, 주)', 'r41/dff50'), ('R3′ (전체 창, 참고)', 'r38/dff_D')):
        C[tag] = [dict(dom=j['dom'], n=float(j['ne']), ratio=j['ff']['ratio'], d2n=j['ff']['d2n'], dmu=j['ff']['dmu']) for j in (json.load(open(f)) for f in sorted(glob.glob(f'{DER}/{sub}/*.json'))) if 'ff' in j]
    doms67 = open(DER + '/r38/p2_67.txt').read().split()
    for tag, key in (('7c aSAMt 320 K', 'aSAMt320-MD320|full|0.75'), ('7c AlphaFlow', 'AlphaFlow-MD320|full|0.75')):
        rows = []
        for d in doms67:
            f = f'{DER}/r38/compare_D2/{d}.json'
            if not os.path.exists(f): continue
            c = json.load(open(f)).get('cells', {}).get(key)
            if c and c.get('d2n'): rows.append(dict(dom=d, n=float((c['D2'] / c['d2n']) ** 2), ratio=c['ratio'], d2n=c['d2n'], dmu=c['dmu']))
        C[tag] = rows
    B6 = [json.load(open(f)) for f in glob.glob(DER + '/r38/b6/*.json')]; B6 = [j for j in B6 if 'skip' not in j and j.get('D') and j.get('level') == 'H']
    for b in sorted({j['bin'] for j in B6}):
        C[f'상동체 TM {b}'] = [dict(dom=j['a'], sf=j['Ha'], n=float(j['nbar']), ratio=j['D']['ratio'], d2n=j['D']['d2n']) for j in B6 if j['bin'] == b]
    R10 = [json.load(open(f)) for f in glob.glob(DER + '/r38/r10a/*.json')]
    C['자연 삽입·결실 (R10a′ 주)'] = [dict(dom=j['a'], sf=j['Ha'], n=float(j['nbar']), ratio=j['D']['ratio'], d2n=j['D']['d2n']) for j in R10 if 'skip' not in j and j.get('main') and j.get('D')]
    return {k: v for k, v in C.items() if v}


def sfmed_ci(rows, q, sfmap, rng):
    v = np.array([r[q] for r in rows], float); sf = [r.get('sf') or sfmap.get(r['dom'], r['dom']) for r in rows]; u = sorted(set(sf)); g = np.array([u.index(s) for s in sf])
    meds = np.array([np.nanmedian(v[g == s]) for s in range(len(u))]); pt = float(np.nanmedian(meds))
    bs = [np.nanmedian(meds[rng.integers(0, len(u), len(u))]) for _ in range(NB)]
    return dict(est=pt, ci95=pct(bs, 2.5, 97.5), n=len(rows), n_sf=len(u))


def region(r): return '가까움' if r < 0.5 else ('교차' if r <= 2 else '멂')


def part_s1():
    sfmap = cath_sf(); rng = np.random.default_rng(20261091); CL = cells(); res = {}
    for k, rows in CL.items():
        m = [r for r in rows if 80 <= r['n'] <= 150]; o = dict(n_units=len(rows), nbar_median=float(np.median([r['n'] for r in rows])), n_units_matched=len(m), nbar_median_matched=float(np.median([r['n'] for r in m])) if m else None)
        for q in ('ratio', 'd2n', 'dmu'):
            if q not in rows[0]: continue
            o[q] = sfmed_ci(rows, q, sfmap, rng); o[q + '_matched'] = sfmed_ci(m, q, sfmap, rng) if len(m) >= 5 else None
            if o[q + '_matched']: o[q + '_matched_over_full'] = o[q + '_matched']['est'] / o[q]['est'] if o[q]['est'] else None
        o['region_full'] = region(o['ratio']['est']); o['region_matched'] = region(o['ratio_matched']['est']) if o.get('ratio_matched') else '길이 맞춤 불가'
        res[k] = o; log('s1', k, o['n_units'], o['n_units_matched'], round(o['ratio']['est'], 3), round(o['ratio_matched']['est'], 3) if o.get('ratio_matched') else None)
    main = [k for k in res if '참고' not in k and res[k].get('ratio_matched')]
    tau = float(kendalltau([res[k]['ratio']['est'] for k in main], [res[k]['ratio_matched']['est'] for k in main]).correlation) if len(main) > 2 else None
    changed = [k for k in main if res[k]['region_full'] != res[k]['region_matched']]
    out = jload(); out['s1'] = dict(written=time.strftime('%F %T'), cells=res, kendall_tau_ratio_full_vs_matched=tau, region_changed=changed, n_cells_compared=len(main),
                                    reading='두 판을 그림 5 에 나란히' if (changed or (tau is not None and tau < 0.9)) else '전체 판 + n̄ 열 본문, 맞춘 판 보충',
                                    rule='n̄ ∈ [80, 150], 상과 중앙값, 상과 부트스트랩 2,000 (seed 20261091). 서술.')
    jdump(out); log('S1_DONE', out['s1']['reading'], tau, changed)


# ---------- s3 ----------
def gnm_edge(ca, cut=7.3):
    N = len(ca); dd = np.linalg.norm(ca[:, None] - ca[None], axis=-1); K = -(dd < cut).astype(float); np.fill_diagonal(K, 0.0); np.fill_diagonal(K, -K.sum(1))
    w, V = np.linalg.eigh(K); Cv = (V[:, 1:] / w[1:]) @ V[:, 1:].T; i = np.arange(N - 1)
    return Cv[i, i] + Cv[i + 1, i + 1] - 2 * Cv[i, i + 1], w[:3]


def anm_one(job):
    d, src = job; fo = f'{OUTA}/{d}.npz'
    if os.path.exists(fo): return d, 'skip', None
    try:
        if src == 'hpool':
            f = f'{P3}/dom2/{d}.npz'
            if not os.path.exists(f) or 'medoid' not in np.load(f).files: return d, 'no medoid', None
            med = np.load(f)['medoid'].astype(np.float64)
        else:
            f = f'{CC}/medoid/{d}.npy'
            if not os.path.exists(f): return d, 'no medoid', None
            med = np.load(f).astype(np.float64)
        ca = np.asarray(ca_trace(med)); o = {}
        for cut in (10.0, 13.0, 15.0): o[f'anm{int(cut)}'], o[f'w{int(cut)}'] = anm_edge_spread(ca, cut)
        o['gnm'], o['wg'] = gnm_edge(ca); np.savez(fo, **o)
        a = f'{CC}/anm/{d}.npz'; rep = float(np.abs(np.load(a)['pred'] - o['anm15']).max()) if os.path.exists(a) else None
        return d, 'ok', rep
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150], None


VARS = ('anm10', 'anm13', 'anm15', 'gnm')


def pair_rows(outdir):
    P = json.load(open(outdir + '/pairs.json'))['pairs']; PB = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(outdir + '/pibar.jsonl')}; sm = {}
    PO = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(outdir + '/pairs_out.jsonl')}
    def small(d):
        if d not in sm:
            z = np.load(f'{outdir}/dom/{d}.npz'); a = f'{OUTA}/{d}.npz'; za = np.load(a) if os.path.exists(a) else None
            sm[d] = (z['degen'] if 'degen' in z.files else None, {v: za[v] for v in VARS} if za is not None else None)
        return sm[d]
    rows = []
    for p in P:
        k = (p[0], p[1]); rec, po = PB.get(k), PO.get(k)
        if rec is None or po is None: continue
        (ga, pa), (gb, pb) = small(p[0]), small(p[1])
        if ga is None or gb is None: continue
        row = dict(ab=[p[0], p[1]], H=p[2])
        for tag, prs in (('main', rec['best']), ('inter', rec['inter'])):
            m = dict(po.get(tag) or {}); prs = np.array(prs, int).reshape(-1, 2); keep = ~(ga[prs[:, 0]] | gb[prs[:, 1]]) if len(prs) else np.zeros(0, bool); sh = prs[keep]
            if len(sh) >= MINSH and pa is not None and pb is not None:
                for v in VARS: m['r_' + v] = spearman(pa[v][sh[:, 0]], pb[v][sh[:, 1]])
            row[tag] = m
        rows.append(row)
    return rows


def part_s3():
    from multiprocessing import Pool
    hd = json.load(open(P3 + '/pairs.json'))['domains']; cd = json.load(open(CC + '/pairs.json'))['domains']; jobs = [(d, 'hpool') for d in hd] + [(d, 'c') for d in cd]; cnt = {}; rep = []; t0 = time.time()
    with Pool(NW) as pool:
        for k, (d, st, r) in enumerate(pool.imap_unordered(anm_one, jobs, chunksize=8)):
            cnt[st[:20]] = cnt.get(st[:20], 0) + 1
            if r is not None: rep.append(r)
            if k % 500 == 0 or st.startswith('ERR'): log('s3 anm', k, '/', len(jobs), d, st, round(time.time() - t0), 's')
    log('s3 anm 상태', cnt, '15 Å 재현 최대 차', max(rep) if rep else None)
    res = dict(anm_status=cnt, anm15_reproduction_maxabs=max(rep) if rep else None, n_reproduction=len(rep))
    for nm, outdir, anmf, seed in (('C', CC, CC + '/anm_pairs.jsonl', 20261074), ('H', P3, CC + '/anm_pairs_hpool.jsonl', 20261077)):
        rows = pair_rows(outdir); st_ = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(anmf)}; diffs = []
        for r in rows:
            s = st_.get(tuple(r['ab']), {}).get('main', {}).get('r_anm'); v = r['main'].get('r_anm15')
            if s is not None and v is not None: diffs.append(abs(s - v))
        o = dict(n_rows=len(rows), r_anm15_vs_stored_maxabs=float(max(diffs)) if diffs else None, n_compared=len(diffs), judge={})
        for v in VARS:
            o['judge'][v] = judge(rows, ('control', 'r_' + v), seed); m = o['judge'][v]['main']; log('s3', nm, v, round(m.get('est', np.nan), 4), [round(m.get('lo', np.nan), 4), round(m.get('hi', np.nan), 4)], o['judge'][v]['passed'])
        res[nm] = o
    allpass = all(res['C']['judge'][v]['passed'] for v in VARS)
    res['reading'] = 'C-2 는 탄성망 매개변수에 매이지 않는다' if allpass else '매개변수 의존 — 서지 않은 판: ' + ', '.join(v for v in VARS if not res['C']['judge'][v]['passed'])
    out = jload(); out['s3'] = dict(written=time.strftime('%F %T'), **res, rule='r38_c.judge (상과 Exp(1) 2,000; 주 판 하한 > 0 이고 교집합 관문 하한 > 0), 통제 = control + r_ANM_판. GNM 절단 7.3 Å.'); jdump(out); log('S3_DONE', res['reading'])


if __name__ == '__main__':
    for p in PARTS:
        log('부분 시작', p); dict(s4=part_s4, s1=part_s1, s3=part_s3)[p]()
    log('R43_ALL_DONE')
