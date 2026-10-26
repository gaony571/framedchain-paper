"""(round41 g1, 2026-10-21 UTC; Drive JSON round41_2026-10-21) 정적 축 분해와 빈칸 몫의 정보량 (주장 1).
부분: static (S1 서술, S2 편 Spearman, S3 ΔAUC; 메도이드 정렬 표) → tm (H-풀·C 쌍의 메도이드 Cα tmtools TM, 긴 쪽 정규화) → rdyn (R1: H-풀 탐색, C 확인)
짝 몫 = share × n_match, 빈칸 몫 = D² − 짝 몫 (= 짝 없는 변의 노름² 합; 표에 짝이 저장돼 있지 않아 대수 항등으로 쓰고, 점검은 부호와 변당 평균 노름 어림).
순위: S2 의 CATH 수준은 동점이 많아(H 12,003) 평균 순위(rankdata)를 쓴다 — argsort 순위는 동점 안 순위를 행 순서가 정한다. R1 은 등록대로 r38_c.analyse_sf(그대로) 가 주이고 평균 순위 판을 점검으로 함께 낸다.
사용: python r41_gapinfo.py <ROOT> <부분들 쉼표> [작업자=4]
출력: derived/r41/tm_med.jsonl, results/r41_gapinfo.json"""
import os
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import sys, json, time, shutil, warnings
import numpy as np
from scipy.stats import spearmanr, rankdata
warnings.filterwarnings('ignore', category=RuntimeWarning)
ROOT = sys.argv[1]; PARTS = sys.argv[2].split(','); NW = int(sys.argv[3]) if len(sys.argv) > 3 else 4
RES = ROOT + '/results'; DER = ROOT + '/derived'; OUT = DER + '/r41'; os.makedirs(OUT, exist_ok=True)
P3, CC, ST = DER + '/r38/p3', DER + '/r38/c', DER + '/r38/static'
NB = 2000; SEED = dict(s2=20261080, cv=20261081, auc=20261082, rh=20261083, rc=20261084)
LV = dict(H=0, T=1, A=2, C=3, X=4); L0 = 3.83; FO = RES + '/r41_gapinfo.json'
_p3 = open(ROOT + '/code/r38_p3.py').read(); exec(_p3[_p3.index('def ranks'):_p3.index('def analyse(')])   # ranks, wls_resid, wcorr, pcorr (그대로)
_c = open(ROOT + '/code/r38_c.py').read(); _af = _c[_c.index('def analyse_sf'):_c.index('def merged')]; exec(_af)   # analyse_sf (1c 등록 구간, 그대로)
aranks = lambda v: rankdata(v) - 1.0
exec(_af.replace('def analyse_sf(', 'def analyse_sf_avg(').replace('ranks(', 'aranks('))


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def jdump(o, p, bk=True):
    json.dump(o, open(p + '.tmp', 'w'), ensure_ascii=False, indent=1, default=float); os.replace(p + '.tmp', p)
    if bk:
        b = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + os.path.basename(p))


def jload():
    return json.load(open(FO)) if os.path.exists(FO) else dict(rule='round41 g1 (Drive JSON round41_2026-10-21)')


def cath_sf():
    h = {}
    for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
        if l.startswith('#') or not l.strip(): continue
        f = l.split(); h[f[0]] = '.'.join(f[1:5])
    return h


def ci(bs, q=(2.5, 97.5)):
    bs = np.asarray(bs, float); bs = bs[np.isfinite(bs)]; return [float(x) for x in np.percentile(bs, q)] if len(bs) else [None, None]


# ---------- static ----------
def static_rows():
    sf = cath_sf(); rows = []; bad = 0
    for l in open(ST + '/table.jsonl'):
        t = json.loads(l); a, b = t['ab']; na, nb = t['ne_a'], t['ne_b']; nbar = (na + nb) / 2; pair = t['share'] * t['n_match']; gap = t['D2'] - pair
        r = dict(ab=[a, b], level=t['level'], lv=LV[t['level']], TM=t['TM'], dfrac=abs(na - nb) / max(na, nb), mfrac=t['n_match'] / min(na, nb), na=na, nb=nb, n_match=t['n_match'],
                 pair_n=pair / np.sqrt(nbar), gap_n=gap / np.sqrt(nbar), D2n=t['D2n'], gap=gap, sa=sf.get(a), sb=sf.get(b))
        if r['TM'] is None or not r['sa'] or not r['sb']: bad += 1; continue
        rows.append(r)
    return rows, bad


def pig_w(ia, ib, v):
    """같은 상과 쌍 = v_s, 다른 상과 쌍 = v_a·v_b (상과 Exp(1) 무게)."""
    return np.where(ia == ib, v[ia], v[ia] * v[ib])


def wauc(s, y, w):
    o = np.argsort(s, kind='mergesort'); s, y, w = s[o], y[o], w[o]; u, inv = np.unique(s, return_inverse=True)
    wp = np.bincount(inv, w * y, len(u)); wn = np.bincount(inv, w * (1 - y), len(u)); cn = np.cumsum(wn) - wn
    return float((wp * (cn + 0.5 * wn)).sum() / (wp.sum() * wn.sum()))


def part_static():
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    rows, bad = static_rows(); out = jload(); n = len(rows)
    sfl = sorted({r['sa'] for r in rows} | {r['sb'] for r in rows}); si = {s: i for i, s in enumerate(sfl)}
    ia = np.array([si[r['sa']] for r in rows]); ib = np.array([si[r['sb']] for r in rows])
    A = lambda k: np.array([r[k] for r in rows], float)
    lv, TM, dfrac, mfrac, gap_n, pair_n, gap = A('lv'), A('TM'), A('dfrac'), A('mfrac'), A('gap_n'), A('pair_n'), A('gap')
    # 점검: 빈칸 몫 ≥ 0, 변당 평균 노름 어림(짝 없는 변 수 × 그 사슬의 변당 평균 노름²)과의 상대 차
    chk = dict(n_rows=n, n_dropped_no_TM_or_sf=bad, gap_min=float(gap.min()), n_gap_neg=int((gap < -1e-9).sum()),
               H_level_sf_mismatch=int(sum(1 for r in rows if r['level'] == 'H' and r['sa'] != r['sb'])), O_level_same_sf=int(sum(1 for r in rows if r['level'] != 'H' and r['sa'] == r['sb'])))
    rng0 = np.random.default_rng(SEED['s2']); pick = rng0.choice(n, min(100, n), replace=False); rel = []
    for q in pick:
        r = rows[q]; xa, xb = medoid(r['ab'][0]), medoid(r['ab'][1])
        if xa is None or xb is None: continue
        ea, eb = (xa ** 2).sum((-1, -2)).mean(), (xb ** 2).sum((-1, -2)).mean(); appx = ea * (r['na'] - r['n_match']) + eb * (r['nb'] - r['n_match'])
        if appx > 0: rel.append(abs(r['gap'] - appx) / appx)
    chk.update(gap_vs_count_approx_rel=dict(n=len(rel), median=float(np.median(rel)) if rel else None, p95=float(np.percentile(rel, 95)) if rel else None,
                                            note='등록한 "상대 1e-9" 확인은 짝이 표에 없어 할 수 없음 → 빈칸 몫 = D² − 짝 몫 은 대수 항등(같은 정렬의 S·U·짝)으로 쓰고, 여기서는 변당 평균 노름 어림과의 차만 본다'))
    # S1 서술
    bins = [0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.01]; S1 = {}
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (TM >= lo) & (TM < hi); key = f'{lo}-{hi}'; S1[key] = dict(levels={})
        for L_, c in LV.items():
            mm = m & (lv == c)
            if mm.sum(): S1[key]['levels'][L_] = {k: dict(median=float(np.median(v[mm])), q25=float(np.percentile(v[mm], 25)), q75=float(np.percentile(v[mm], 75))) for k, v in (('pair_n', pair_n), ('gap_n', gap_n), ('D2n', A('D2n')))} | dict(n=int(mm.sum()))
        if m.sum() >= 10 and len(set(lv[m])) >= 2:
            S1[key]['spearman_level'] = {k: float(spearmanr(lv[m], v[m]).correlation) for k, v in (('pair_n', pair_n), ('gap_n', gap_n))}
    # S2 편 Spearman (평균 순위)
    R = {k: aranks(v) for k, v in dict(lv=lv, TM=TM, dfrac=dfrac, mfrac=mfrac, gap_n=gap_n, pair_n=pair_n).items()}; one = np.ones(n)
    def s2(y, zs, seed):
        est = pcorr(R[y], R['lv'], [R[z] for z in zs], one); rng = np.random.default_rng(seed); bs = []
        for _ in range(NB): bs.append(pcorr(R[y], R['lv'], [R[z] for z in zs], pig_w(ia, ib, rng.exponential(1, len(sfl)))))
        return dict(est=float(est), ci95=ci(bs), ci90=ci(bs, (5, 95)), se=float(np.std(bs)), z=list(zs))
    S2 = dict(gap_main=s2('gap_n', ('TM', 'dfrac'), SEED['s2']), pair=s2('pair_n', ('TM', 'dfrac'), SEED['s2']), gap_plus_mfrac=s2('gap_n', ('TM', 'dfrac', 'mfrac'), SEED['s2']),
              n=n, n_sf=len(sfl), note='수준 = H<T<A<C<X (0–4), 평균 순위, 점추정 무게 1, 상과 무게(같은 상과 v_s, 다른 상과 v_a·v_b) 2,000. gap_plus_mfrac 은 "짝 없는 변의 변당 노름" 만 보는 판.')
    log('S2', S2['gap_main'], 'pair', S2['pair']['est'])
    # S3 상과 나눈 교차 검증 AUC
    y = (lv == 0).astype(float); feats = dict(M0=[TM], M1=[TM, dfrac], M2=[TM, dfrac, gap_n], M2p=[TM, dfrac, pair_n])
    X = {m: np.column_stack(f) for m, f in feats.items()}; ssum = {m: np.zeros(n) for m in X}; cnt = np.zeros(n); rng = np.random.default_rng(SEED['cv'])
    for rep in range(20):
        fold = rng.integers(0, 5, len(sfl)); fa, fb = fold[ia], fold[ib]
        for k in range(5):
            te = (fa == k) & (fb == k); tr = (fa != k) & (fb != k)
            if te.sum() == 0 or len(set(y[tr])) < 2: continue
            for m in X:
                sc = StandardScaler().fit(X[m][tr]); cl = LogisticRegression(max_iter=2000).fit(sc.transform(X[m][tr]), y[tr]); ssum[m][te] += cl.predict_proba(sc.transform(X[m][te]))[:, 1]
            cnt[te] += 1
    ok = cnt > 0; pr = {m: ssum[m][ok] / cnt[ok] for m in X}; yo = y[ok]; iao, ibo = ia[ok], ib[ok]
    auc = {m: wauc(pr[m], yo, np.ones(ok.sum())) for m in X}; rng = np.random.default_rng(SEED['auc']); bd, bdp = [], []
    for _ in range(NB):
        w = pig_w(iao, ibo, rng.exponential(1, len(sfl))); a1 = wauc(pr['M1'], yo, w); bd.append(wauc(pr['M2'], yo, w) - a1); bdp.append(wauc(pr['M2p'], yo, w) - a1)
    S3 = dict(auc=auc, dAUC_gap=dict(est=auc['M2'] - auc['M1'], ci95=ci(bd), ci90=ci(bd, (5, 95))), dAUC_pair=dict(est=auc['M2p'] - auc['M1'], ci95=ci(bdp)),
              n_tested=int(ok.sum()), n_never_tested=int((~ok).sum()), n_H_tested=int(yo.sum()), mean_times_tested=float(cnt[ok].mean()))
    log('S3', S3['auc'], S3['dAUC_gap'])
    v2 = S2['gap_main']['ci95'][0] is not None and S2['gap_main']['ci95'][0] > 0; v3 = S3['dAUC_gap']['ci95'][0] is not None and S3['dAUC_gap']['ci95'][0] > 0
    out['static'] = dict(written=time.strftime('%F %T'), checks=chk, S1=S1, S2=S2, S3=S3, pred_S2=bool(v2), pred_S3=bool(v3),
                         sentence='빈칸 몫은 TM·길이 차 너머의 정보를 담는다' if (v2 and v3) else '주장 1 은 "관행이 버리는 양의 크기" 로만 쓴다')
    jdump(out, FO); log('G1_STATIC_DONE', out['static']['sentence'])


# ---------- tm ----------
_MED = {}
def medoid(d):
    if d in _MED: return _MED[d]
    x = None; f = f'{P3}/dom2/{d}.npz'
    if os.path.exists(f):
        z = np.load(f)
        if 'medoid' in z.files: x = z['medoid'].astype(np.float64)
    if x is None:
        for f in (f'{ST}/medoid/{d}.npy', f'{CC}/medoid/{d}.npy'):
            if os.path.exists(f): x = np.load(f).astype(np.float64); break
    if len(_MED) > 400: _MED.clear()
    _MED[d] = x; return x


def ca(x):
    e = 2.0 * L0 * x[:, :, 0]; return np.ascontiguousarray(np.concatenate([np.zeros((1, 3)), np.cumsum(e, 0)], 0))


def tm_one(ab):
    from tmtools import tm_align
    a, b = ab; xa, xb = medoid(a), medoid(b)
    if xa is None or xb is None: return dict(ab=[a, b], err='메도이드 없음')
    pa, pb = ca(xa), ca(xb); r = tm_align(pa, pb, 'A' * len(pa), 'A' * len(pb))
    return dict(ab=[a, b], tm_long=float(r.tm_norm_chain1 if len(pa) >= len(pb) else r.tm_norm_chain2), tm1=float(r.tm_norm_chain1), tm2=float(r.tm_norm_chain2), La=len(pa), Lb=len(pb))


def part_tm():
    from multiprocessing import Pool
    fo = OUT + '/tm_med.jsonl'; done = {tuple(json.loads(l)['ab']) for l in open(fo)} if os.path.exists(fo) else set()
    prs = []
    for base in (P3, CC):
        for l in open(base + '/pibar.jsonl'):
            k = tuple(json.loads(l)['ab'])
            if k not in done: prs.append(k)
    prs = sorted(set(prs)); log('tm', len(prs), '남음', len(done), '끝남'); t0 = time.time()
    with Pool(NW) as pool, open(fo, 'a') as fh:
        for i, r in enumerate(pool.imap(tm_one, prs, chunksize=32)):
            fh.write(json.dumps(r) + '\n')
            if i % 2000 == 0: fh.flush(); log('tm', i, round(time.time() - t0), 's')
    log('G1_TM_DONE')


# ---------- rdyn ----------
def rdyn_rows(base, tmd):
    pb = {}
    for l in open(base + '/pibar.jsonl'):
        j = json.loads(l); pb[tuple(j['ab'])] = j
    nm = {}
    def norms(d):
        if d not in nm: nm[d] = (np.load(f'{base}/dom/{d}.npz')['m'].astype(np.float64) ** 2).sum((-1, -2))
        return nm[d]
    rows = []
    for l in open(base + '/pairs_out.jsonl'):
        r = json.loads(l); k = tuple(r['ab'])
        if k not in pb or not r.get('main') or r['main'].get('r_dyn') is None: continue
        na_, nb_ = norms(k[0]), norms(k[1]); best = pb[k]['best']; ua = np.ones(len(na_), bool); ub = np.ones(len(nb_), bool)
        for i, j in best: ua[i] = False; ub[j] = False
        nbar = (len(na_) + len(nb_)) / 2
        ex = dict(gap_n=float((na_[ua].sum() + nb_[ub].sum()) / np.sqrt(nbar)), dfrac=abs(len(na_) - len(nb_)) / max(len(na_), len(nb_)), tm_med=tmd.get(k), tm_table=r.get('tm_long'))
        for tag in ('main', 'inter'):
            if r.get(tag) is not None: r[tag] = {**r[tag], **ex}
        rows.append(r)
    return rows


def part_rdyn():
    tmd = {}
    for l in open(OUT + '/tm_med.jsonl'):
        j = json.loads(l)
        if 'tm_long' in j: tmd[tuple(j['ab'])] = j['tm_long']
    out = jload(); z = ('tm_med', 'dfrac'); res = {}
    for nm_, base, seed in (('H', P3, SEED['rh']), ('C', CC, SEED['rc'])):
        rows = rdyn_rows(base, tmd); o = dict(n_rows=len(rows), n_no_tm=sum(1 for r in rows if r['main'].get('tm_med') is None))
        o['main'] = analyse_sf(rows, 'main', x='gap_n', z=z, seed=seed); o['main_avgrank'] = analyse_sf_avg(rows, 'main', x='gap_n', z=z, seed=seed)
        o['sens'] = dict(plus_shared_frac=analyse_sf(rows, 'main', x='gap_n', z=z + ('shared_frac',), seed=seed + 10), plus_control=analyse_sf(rows, 'main', x='gap_n', z=z + ('control',), seed=seed + 11),
                         inter=analyse_sf(rows, 'inter', x='gap_n', z=z, seed=seed + 12))
        if nm_ == 'H':
            o['sens']['tm_table'] = analyse_sf(rows, 'main', x='gap_n', z=('tm_table', 'dfrac'), seed=seed + 13)
            pr = [(r['main']['tm_med'], r['main']['tm_table']) for r in rows if r['main'].get('tm_med') is not None and r['main'].get('tm_table') is not None]
            o['check_tm_med_vs_table'] = dict(n=len(pr), spearman=float(spearmanr(*zip(*pr)).correlation) if len(pr) > 2 else None, median_diff=float(np.median([a - b for a, b in pr])) if pr else None)
        res[nm_] = o; log('R1', nm_, o['main'], 'avg', o['main_avgrank'].get('est'))
    h, c = res['H']['main'], res['C']['main']; sg = np.sign(h['est']) if h.get('est') is not None else 0
    conf = bool(c.get('lo') is not None and sg != 0 and ((sg > 0 and c['lo'] > 0) or (sg < 0 and c['hi'] < 0)))
    res['confirm'] = dict(direction='양' if sg > 0 else '음', C_ci95=[c.get('lo'), c.get('hi')], verdict='독립 상과에서 확인' if conf else '확인 안 됨')
    out['rdyn'] = dict(written=time.strftime('%F %T'), **res, note='x = 평균 점 π̄ 빈칸 몫/√n̄, y = r_dyn, z = TM_med(메도이드 tmtools 긴 쪽) + dfrac; r38_c.analyse_sf 그대로(주), 평균 순위 판은 점검')
    jdump(out, FO); log('G1_RDYN_DONE', res['confirm'])


if __name__ == '__main__':
    for p in PARTS:
        log('부분 시작', p); dict(static=part_static, tm=part_tm, rdyn=part_rdyn)[p]()
    log('G1_ALL_DONE')
