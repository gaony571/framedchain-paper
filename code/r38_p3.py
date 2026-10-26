"""(prereg round36 1b, 1차 P3) 동역학 보존: H-풀(상동체) 쌍에서 편 Spearman(r_dyn, 서열 동일성 | 공유 핵심부 거리).
부분: pairs → dom → pibar → pair → stats   (사용: python r38_p3.py <ROOT> <부분들 쉼표 | all> [작업자 수=3])
  pairs : results/pilot_mdcath_homolog_pairs.tsv 14,189 쌍 중 두 도메인 모두 raw/mdcath 에 h5 가 있는 쌍(수·날짜 박은 도메인 목록·상과 수),
          상과당 50 쌍 상한 판(상과마다 정렬 뒤 default_rng(20261027).permutation 앞 50).
  dom   : 도메인마다 320 K 각 반복의 뒤 절반에서 접힌 프레임(첫 320 K 프레임 천연 접촉 Q ≥ 0.75; r37_ff·r38 과 같음).
          접힌 프레임 ≥ 40 인 반복마다 최대 100 개 무작위(도메인 이름 seed), 그런 반복이 둘 이상이어야 씀(아니면 뺌, 수 기록).
          평균 점 m = fs.mean_point, 사영 평균 점 = 변마다 m_i 의 SVD 로 가장 가까운 회전 R̄_i × 평균 길이 비 s̄_i / 2,
          퇴화 변 = m_i 의 가장 작은 특이값 < 0.5 × s̄_i / 2. 평균에 맞춘 프레임은 /content/p3cache 에(지워지면 같은 seed 로 다시 만듦).
  pibar : 쌍마다 두 평균 점 사이 정렬(항등 + 무작위 127 출발, 번갈아 풀기 10 번; fs.align 과 같은 계산을 묶음으로, GPU).
          π̄ = 최적 출발의 짝. 거의 최적 = S ≥ (1 − 0.005) S_최적 인 출발들, 교집합 = 그 모두에서 같은 짝.
  pair  : 공유 변 = π̄ 짝 − 어느 쪽이든 퇴화 변. 30 개 미만이면 뺌. 도메인마다 프레임을 공유 변 위 GPA(자기 평균 둘레, 3 번)로 다시 맞춰
          변별 퍼짐 윤곽 → r_dyn = 두 윤곽의 Spearman. 통제 = 사영 평균 점 사이 d_fix 짝 몫(공유 짝 위 회전만) ÷ 짝 수.
          동일성 = 공유 짝(변 i 의 앞 잔기) 위 서열 동일성. 교집합 판도 같은 방식(교집합 변 위). 민감도용 값도 함께(전체 변 중첩 윤곽,
          평균 점 그대로의 짝 몫, 평균 점 D/√n̄, 공유 변 분율, 퇴화 비율).
  stats : 세 변수를 순위로 바꾼 가중 편 Pearson (쌍 무게 = 1/상과 쌍 수). 교차 군집 부트스트랩 2,000 (도메인 w ~ Exp(1), 상과 v ~ Exp(1),
          쌍 무게 = v_s w_i w_j / 상과 쌍 수; seed 20261028). 판정 표(최소 효과 0.15; ±0.10 동등은 보고), 공동 관문(교집합 ≥ 30 변 쌍에서
          주 판·교집합 판), 설계 효과, Freedman–Lane null 1,000(점검용), 민감도(상과 고정 효과, 다섯 층, 공유 변 분율 공변량, 상한 50,
          통제 바꾸기, 전체 변 중첩 윤곽, 퇴화 변 포함).
          아직 안 넣은 민감도: 메도이드·결정 구조 통제, Q 0.7·0.8, 반복 나눔 신뢰도 보정, 로그 퍼짐 비 반응, 281 쌍 판.
출력: derived/r38/p3/{pairs.json, dom/<d>.npz, pibar.jsonl, pairs_out.jsonl}, results/r38_p3.json (사본 results/_backup/<날짜>/)"""
import os, sys
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
import json, glob, time, zlib, itertools, shutil, traceback, difflib
ROOT = sys.argv[1]; PARTS = sys.argv[2].split(',') if len(sys.argv) > 2 else ['all']; NW = int(sys.argv[3]) if len(sys.argv) > 3 else 3
sys.path[:0] = [ROOT + '/framespace_pkg']
import numpy as np, torch, h5py
import framespace as fs
from framespace.distance import dp_score
from framespace.geometry import procrustes, random_rotations
torch.set_num_threads(1)
src = open(ROOT + '/code/r37_ff.py').read(); exec(src[src.index('AA = dict('):src.index('def stats')])   # md_atoms
OUT = ROOT + '/derived/r38/p3'; RES = ROOT + '/results'; CACHE = '/content/p3cache'
for p in (OUT, OUT + '/dom', CACHE): os.makedirs(p, exist_ok=True)
QT = 0.75; NMAX = 100; NMIN = 40; MINSH = 30; EPS = 0.005; NB = 2000; NPERM = 1000


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)
def backup(p):
    b = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + os.path.basename(p))
def jdump(o, p, bk=False):
    json.dump(o, open(p + '.tmp', 'w'), ensure_ascii=False, indent=1, default=float); os.replace(p + '.tmp', p)
    if bk: backup(p)
def frames(X): return np.asarray(fs.coordinates_from_backbone(torch.as_tensor(X[:, :, 0]), torch.as_tensor(X[:, :, 1]), torch.as_tensor(X[:, :, 2])))


# ---------- pairs ----------
def part_pairs():
    import pandas as pd
    hp = pd.read_csv(RES + '/pilot_mdcath_homolog_pairs.tsv', sep='\t')
    h5 = {os.path.basename(x)[len('mdcath_dataset_'):-3] for x in glob.glob(ROOT + '/raw/mdcath/mdcath_dataset_*.h5')}
    ok = hp[hp.a.isin(h5) & hp.b.isin(h5)].reset_index(drop=True)
    rng = np.random.default_rng(20261027); cap = []
    for H, g in sorted(ok.groupby('H'), key=lambda x: x[0]):
        g = g.sort_values(['a', 'b']); idx = list(g.index)
        cap += idx if len(idx) <= 50 else [idx[i] for i in rng.permutation(len(idx))[:50]]
    doms = sorted(set(ok.a) | set(ok.b))
    jdump(dict(written=time.strftime('%F %T'), n_pairs_table=len(hp), n_pairs=len(ok), n_sf=int(ok.H.nunique()), n_domains=len(doms), domains=doms,
               pairs=[[r.a, r.b, r.H, float(r.ident), int(r.La), int(r.Lb), float(max(r.tm_a, r.tm_b))] for r in ok.itertuples()], cap50=sorted(int(i) for i in cap),
               rule='두 도메인 모두 raw/mdcath h5 가 있는 쌍. 상한 판 = 상과마다 (a, b) 정렬 뒤 default_rng(20261027).permutation 앞 50'), OUT + '/pairs.json', bk=True)
    log('pairs', len(hp), '→', len(ok), '쌍,', ok.H.nunique(), '상과,', len(doms), '도메인, 상한 판', len(cap))


# ---------- dom ----------
def dom_build(d):
    """도메인 앙상블(평균에 맞춘 프레임)과 요약. 같은 seed 라 다시 만들어도 같다."""
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); g = h[d]; smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode()); L = len(smd)
    reps = [g['320'][r]['coords'][:][:, imd.ravel()].reshape(-1, L, 3, 3) for r in sorted(g['320'].keys(), key=int)]; h.close()
    ca0 = reps[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
    rng = np.random.default_rng(zlib.crc32(f'r38_p3|{d}'.encode())); runs, lab = [], []
    for r, X in enumerate(reps):
        Xh = X[len(X) // 2:]; idx = np.where(Qf(Xh) >= QT)[0]
        if len(idx) >= NMIN:
            sel = np.sort(rng.choice(idx, min(NMAX, len(idx)), replace=False)); runs.append(Xh[sel]); lab += [r] * len(sel)
    if len(runs) < 2: return None, dict(dom=d, L=L, skip=f'접힌 프레임 ≥ {NMIN} 인 반복 {len(runs)} 개')
    E = frames(np.concatenate(runs)); m = np.asarray(fs.mean_point(E)); A = np.asarray(fs.align_to_reference(E, torch.as_tensor(m)), dtype=np.float32)
    s = 2 * np.linalg.norm(A[:, :, :, 0], axis=-1); sbar = s.mean(0)
    U, S, Vt = np.linalg.svd(m); dt = np.sign(np.linalg.det(U @ Vt)); Dg = np.ones((len(m), 3)); Dg[:, 2] = dt
    Rbar = np.einsum('iab,ib,ibc->iac', U, Dg, Vt); Pm = Rbar * (sbar / 2)[:, None, None]
    degen = S[:, 2] < 0.5 * sbar / 2
    info = dict(dom=d, L=L, ne=L - 1, n_frames=int(len(A)), n_runs=len(runs), seq=smd)
    return dict(A=A, rep=np.array(lab), m=m, Pm=Pm, sbar=sbar, degen=degen, sv=S), info


def dom_get(d):
    f = f'{CACHE}/{d}.npz'
    if os.path.exists(f): z = np.load(f); return {k: z[k] for k in z.files}
    r, info = dom_build(d)
    if r is None: return None
    np.savez(f, **r); return r


def _dom_worker(d):
    fo = f'{OUT}/dom/{d}.npz'
    if os.path.exists(fo): return d, 'skip'
    try:
        r, info = dom_build(d)
        if r is None: np.savez(fo, info=json.dumps(info, ensure_ascii=False)); return d, info['skip']
        np.savez(f'{CACHE}/{d}.npz', **r)
        np.savez(fo, m=r['m'], Pm=r['Pm'], sbar=r['sbar'], degen=r['degen'], sv=r['sv'], info=json.dumps(info, ensure_ascii=False)); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


def part_dom():
    doms = json.load(open(OUT + '/pairs.json'))['domains']; from multiprocessing import Pool
    with Pool(NW) as p:
        for k, (d, st) in enumerate(p.imap_unordered(_dom_worker, doms)):
            if st != 'skip' and (k % 50 == 0 or st != 'ok'): log('dom', k, d, st)
    n_ok = sum(1 for d in doms if 'm' in np.load(f'{OUT}/dom/{d}.npz').files)
    log('DOM_DONE', n_ok, '/', len(doms))


def dom_small(d):
    z = np.load(f'{OUT}/dom/{d}.npz'); info = json.loads(str(z['info']))
    return (dict(m=z['m'], Pm=z['Pm'], sbar=z['sbar'], degen=z['degen'], **info) if 'm' in z.files else None), info


# ---------- pibar ----------
def part_pibar():
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'; P = json.load(open(OUT + '/pairs.json'))['pairs']; fo = OUT + '/pibar.jsonl'
    done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    D = {}
    for d in sorted({p[0] for p in P} | {p[1] for p in P}):
        s_, _ = dom_small(d)
        if s_ is not None: D[d] = s_['m']
    todo = [p for p in P if (p[0], p[1]) not in done and p[0] in D and p[1] in D]
    todo.sort(key=lambda p: (len(D[p[0]]) // 16, len(D[p[1]]) // 16)); log('pibar', len(todo), '쌍 남음', dev)
    BS = 48; t0 = time.time()
    for b0 in range(0, len(todo), BS):
        B = todo[b0:b0 + BS]; n = max(len(D[p[0]]) for p in B); m = max(len(D[p[1]]) for p in B)
        x = fs.pad_stack([D[p[0]] for p in B], n).to(dev, torch.float64); y = fs.pad_stack([D[p[1]] for p in B], m).to(dev, torch.float64); Pn = len(B)
        I = torch.eye(3, dtype=x.dtype, device=dev).expand(Pn, 1, 3, 3)
        U = torch.cat([I, random_rotations(Pn * 127, seed=b0, dtype=x.dtype, device=dev).view(Pn, 127, 3, 3)], 1)
        best = torch.full((Pn, 128), -1e30, dtype=x.dtype, device=dev); bG = None
        for it in range(11):
            with torch.no_grad(): K = torch.einsum('pnxz,psmxz->psnm', x, torch.einsum('psxy,pmyz->psmxz', U, y))
            K = K.detach().requires_grad_(True)
            with torch.enable_grad(): Sc = dp_score(K); G, = torch.autograd.grad(Sc.sum(), K)
            with torch.no_grad():
                imp = Sc > best + 1e-7 * (1 + best.abs()); best = torch.where(imp, Sc, best)
                bG = G if bG is None else torch.where(imp[..., None, None], G, bG)
                if it == 10 or not imp.any(): break
                U = torch.where(imp[..., None, None], procrustes(torch.einsum('psnm,pmxz,pnyz->psxy', G, y, x)), U)
        Mt = (bG > 0.5).cpu().numpy(); bs = best.cpu().numpy()
        with open(fo, 'a') as fh:
            for q, p in enumerate(B):
                na, nb = len(D[p[0]]), len(D[p[1]]); Smax = float(bs[q].max()); k = int(bs[q].argmax())
                near = np.where(bs[q] >= (1 - EPS) * Smax)[0]; Mq = Mt[q][:, :na, :nb]
                inter = np.logical_and.reduce(Mq[near], axis=0); bi, bj = np.nonzero(Mq[k]); ii_, jj_ = np.nonzero(inter)
                nd = len({Mq[s_].tobytes() for s_ in near})
                fh.write(json.dumps(dict(ab=[p[0], p[1]], S=Smax, n_near=int(len(near)), n_distinct=nd, best=[[int(a), int(c)] for a, c in zip(bi, bj)], inter=[[int(a), int(c)] for a, c in zip(ii_, jj_)])) + '\n')
        if (b0 // BS) % 20 == 0: log('pibar', b0 + len(B), '/', len(todo), round(time.time() - t0), 's')
    log('PIBAR_DONE')


# ---------- pair ----------
def gpa_profile(F, idx, iters=3):
    """F (n, ne, 3, 3) 평균에 맞춘 프레임, idx 공유 변 → 공유 변 위 GPA(자기 평균 둘레) 뒤의 변별 퍼짐(len(idx),)."""
    X = F[:, idx].astype(np.float64); mu = X.mean(0)
    for _ in range(iters):
        A = np.einsum('kiab,icb->kac', X, mu); W, S, Zt = np.linalg.svd(A); dt = np.sign(np.linalg.det(np.einsum('kba,kcb->kac', Zt, W)))
        Dg = np.ones((len(X), 3)); Dg[:, 2] = dt; Uk = np.einsum('kba,kb,kcb->kac', Zt, Dg, W)   # U = Z D Wᵀ
        X = np.einsum('kac,kicb->kiab', Uk, X); mu = X.mean(0)
    return ((X - mu) ** 2).sum((-1, -2)).mean(0)


def rot_share(Pa, Pb):
    """Σ‖Pa_i − U Pb_i‖² 를 회전만으로 최소 (짝 순서대로) → 짝 수로 나눈 값."""
    A = np.einsum('iab,icb->ac', Pb, Pa); W, S, Zt = np.linalg.svd(A); D_ = np.diag([1, 1, np.sign(np.linalg.det(Zt.T @ W.T))]); U = Zt.T @ D_ @ W.T
    return float(((Pa - np.einsum('ac,icb->iab', U, Pb)) ** 2).sum() / len(Pa))


def spearman(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b)); return float(np.corrcoef(ra, rb)[0, 1])


def one_pair(a, b, rec, Da, Db, Fa, Fb):
    out = dict(ab=[a, b], S_mean=rec['S'], n_near=rec['n_near'], n_distinct=rec['n_distinct'])
    for tag, prs in (('main', rec['best']), ('inter', rec['inter'])):
        prs = np.array(prs, int).reshape(-1, 2); keep = ~(Da['degen'][prs[:, 0]] | Db['degen'][prs[:, 1]]) if len(prs) else np.zeros(0, bool)
        sh = prs[keep]; o = dict(n_pairs=int(len(prs)), n_shared=int(len(sh)), degen_frac=float(1 - keep.mean()) if len(prs) else None)
        if len(sh) >= MINSH:
            va, vb = gpa_profile(Fa, sh[:, 0]), gpa_profile(Fb, sh[:, 1])
            o.update(r_dyn=spearman(va, vb), control=rot_share(Da['Pm'][sh[:, 0]], Db['Pm'][sh[:, 1]]),
                     ident=float(np.mean([Da['seq'][i] == Db['seq'][j] for i, j in sh])), shared_frac=float(len(sh) / min(Da['ne'], Db['ne'])),
                     control_meanpt=rot_share(Da['m'][sh[:, 0]], Db['m'][sh[:, 1]]),
                     r_dyn_fullsup=spearman(((Fa[:, sh[:, 0]] - Fa[:, sh[:, 0]].mean(0)) ** 2).sum((-1, -2)).mean(0), ((Fb[:, sh[:, 1]] - Fb[:, sh[:, 1]].mean(0)) ** 2).sum((-1, -2)).mean(0)))
            if tag == 'main':
                pa = prs; va2, vb2 = gpa_profile(Fa, pa[:, 0]), gpa_profile(Fb, pa[:, 1])
                o.update(r_dyn_withdegen=spearman(va2, vb2), control_withdegen=rot_share(Da['Pm'][pa[:, 0]], Db['Pm'][pa[:, 1]]))
        out[tag] = o
    d2 = float((Da['m'] ** 2).sum() + (Db['m'] ** 2).sum() - 2 * rec['S']); out['meanpt_D_scaled'] = float(np.sqrt(max(d2, 0)) / np.sqrt(0.75 * (Da['ne'] + Db['ne'])))
    return out


def part_pair():
    P = json.load(open(OUT + '/pairs.json'))['pairs']; fo = OUT + '/pairs_out.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    PB = {}
    for l in open(OUT + '/pibar.jsonl'):
        r = json.loads(l); PB[tuple(r['ab'])] = r
    todo = sorted([p for p in P if (p[0], p[1]) in PB and (p[0], p[1]) not in done], key=lambda p: (p[0], p[1])); log('pair', len(todo), '쌍 남음')
    small, big = {}, {}
    def get(d):
        if d not in small: small[d] = dom_small(d)[0]
        if d not in big:
            if len(big) > 40: big.pop(next(iter(big)))
            z = dom_get(d); big[d] = None if z is None else z['A']
        return small[d], big[d]
    t0 = time.time()
    with open(fo, 'a') as fh:
        for k, p in enumerate(todo):
            try:
                Da, Fa = get(p[0]); Db, Fb = get(p[1])
                o = one_pair(p[0], p[1], PB[(p[0], p[1])], Da, Db, Fa, Fb); o.update(H=p[2], ident_table=p[3], tm_long=p[6])
                fh.write(json.dumps(o) + '\n')
            except Exception as e:
                log('ERR pair', p[:2], repr(e)[:150])
            if k % 500 == 0: log('pair', k, '/', len(todo), round(time.time() - t0), 's')
    log('PAIR_DONE')


# ---------- stats ----------
def ranks(v): return np.argsort(np.argsort(v)).astype(float)
def wls_resid(y, X, w):
    Xd = np.column_stack([np.ones(len(y))] + list(X)); W = np.sqrt(w); beta = np.linalg.lstsq(Xd * W[:, None], y * W, rcond=None)[0]; return y - Xd @ beta
def wcorr(a, b, w):
    w = w / w.sum(); ma, mb = (w * a).sum(), (w * b).sum(); return float((w * (a - ma) * (b - mb)).sum() / np.sqrt((w * (a - ma) ** 2).sum() * (w * (b - mb) ** 2).sum()))
def pcorr(ry, rx, Z, w): return wcorr(wls_resid(ry, Z, w), wls_resid(rx, Z, w), w)


def analyse(rows, key='main', y='r_dyn', x='ident', z=('control',), seed=20261028, boot=True, sf_fixed=False):
    rows = [r for r in rows if r[key].get(y) is not None]
    if len(rows) < 10: return dict(n=len(rows))
    Y = ranks(np.array([r[key][y] for r in rows])); X_ = ranks(np.array([r[key][x] for r in rows])); Zs = [ranks(np.array([r[key][c] if c in r[key] else r[c] for r in rows])) for c in z]
    H = np.array([r['H'] for r in rows]); cnt = {h: (H == h).sum() for h in set(H)}; base = np.array([1.0 / cnt[h] for h in H])
    if sf_fixed:
        for arr in [Y, X_] + Zs:
            for h in set(H): s_ = H == h; arr[s_] -= arr[s_].mean()
    est = pcorr(Y, X_, Zs, base); out = dict(n=len(rows), n_sf=len(cnt), est=est)
    if boot:
        rng = np.random.default_rng(seed); doms = sorted({r['ab'][0] for r in rows} | {r['ab'][1] for r in rows}); di = {d: i for i, d in enumerate(doms)}
        ia = np.array([di[r['ab'][0]] for r in rows]); ib = np.array([di[r['ab'][1]] for r in rows]); sfl = sorted(cnt); si = np.array([sfl.index(h) for h in H]); bs = []
        for _ in range(NB):
            w_ = rng.exponential(1, len(doms)); v_ = rng.exponential(1, len(sfl)); bs.append(pcorr(Y, X_, Zs, base * v_[si] * w_[ia] * w_[ib]))
        bs = np.array(bs); lo, hi = np.percentile(bs, [2.5, 97.5])
        zb = np.arctanh(np.clip(bs, -0.999, 0.999)); deff = float(zb.var() / (1.06 / (len(rows) - 3)))
        out.update(lo=float(lo), hi=float(hi), design_effect=deff, n_eff=float((len(rows) - 3) / deff + 3))
    return out


def verdict(m, gate):
    if m.get('lo') is None: return '판정 불가(쌍 부족)'
    if m['lo'] > 0 and m['hi'] >= 0.15: return '효과 있음' if gate.get('lo') is not None and gate['lo'] > 0 else '판정 불가(정렬 교란을 배제하지 못함: 교집합 판이 서지 않음)'
    if m['lo'] > 0: return '있지만 최소 효과(0.15)보다 작다'
    if m['lo'] > -0.15 and m['hi'] < 0.15: return '동등(±0.15)'
    return '판정 불가(검출되지 않음)'


def part_stats():
    rows = [json.loads(l) for l in open(OUT + '/pairs_out.jsonl')]; Pj = json.load(open(OUT + '/pairs.json'))
    capset = {tuple(Pj['pairs'][i][:2]) for i in Pj['cap50']}
    main = analyse(rows, 'main'); both = [r for r in rows if r['main'].get('r_dyn') is not None and r['inter'].get('r_dyn') is not None]
    gate_main = analyse(both, 'main', seed=20261029); gate_inter = analyse(both, 'inter', seed=20261030)
    v = verdict(main, gate_inter)
    # Freedman–Lane null (점검용): 동일성 순위를 통제 순위에 회귀한 잔차를 상과 안에서 순열
    R_ = [r for r in rows if r['main'].get('r_dyn') is not None]; Y = ranks(np.array([r['main']['r_dyn'] for r in R_])); Xr = ranks(np.array([r['main']['ident'] for r in R_])); Zr = ranks(np.array([r['main']['control'] for r in R_]))
    H = np.array([r['H'] for r in R_]); cnt = {h: (H == h).sum() for h in set(H)}; base = np.array([1.0 / cnt[h] for h in H])
    Zd = np.column_stack([np.ones(len(Zr)), Zr]); Wq = np.sqrt(base); beta = np.linalg.lstsq(Zd * Wq[:, None], Xr * Wq, rcond=None)[0]; fit = Zd @ beta; res = Xr - fit
    rng = np.random.default_rng(20261031); nulls = []
    groups = {h: np.where(H == h)[0] for h in cnt}
    for _ in range(NPERM):
        rp = res.copy()
        for h, ix in groups.items(): rp[ix] = res[rng.permutation(ix)]
        nulls.append(pcorr(Y, fit + rp, [Zr], base))
    nulls = np.array(nulls); fl = dict(null_lo=float(np.percentile(nulls, 2.5)), null_hi=float(np.percentile(nulls, 97.5)), p_two=float((np.abs(nulls) >= abs(main['est'])).mean()))
    # 민감도
    sens = dict(sf_fixed=analyse(rows, sf_fixed=True, seed=1), shared_frac_cov=analyse(rows, z=('control', 'shared_frac'), seed=2),
                cap50=analyse([r for r in rows if tuple(r['ab']) in capset], seed=3), control_meanpt=analyse(rows, z=('control_meanpt',), seed=4),
                control_tm=analyse(rows, z=('tm_long',), seed=5), control_meanptD=analyse(rows, z=('meanpt_D_scaled',), seed=6),
                fullsup_profile=analyse(rows, y='r_dyn_fullsup', seed=7), with_degen=analyse(rows, y='r_dyn_withdegen', z=('control_withdegen',), seed=8))
    # 다섯 층(통제의 사분위 경계 + 상위 10%)
    c = np.array([r['main']['control'] for r in R_]); qs = np.percentile(c, [25, 50, 75, 90]); st = np.digitize(c, qs); strata = []
    for k in range(5):
        ix = np.where(st == k)[0]
        if len(ix) > 10: strata.append(dict(stratum=k, n=int(len(ix)), w=float(base[ix].sum()), rho=wcorr(ranks(Y[ix]), ranks(Xr[ix]), base[ix])))
    sens['strata'] = dict(per=strata, combined=float(sum(s_['w'] * s_['rho'] for s_ in strata) / sum(s_['w'] for s_ in strata)) if strata else None)
    drop = dict(n_pairs=len(Pj['pairs']), with_pibar=len(rows), main_ok=len(R_), inter_ok=len(both),
                main_fail_minshared=sum(1 for r in rows if r['main'].get('r_dyn') is None), inter_fail_minshared=sum(1 for r in rows if r['main'].get('r_dyn') is not None and r['inter'].get('r_dyn') is None))
    idq = lambda rs: np.percentile([r['main']['ident'] for r in rs], [5, 25, 50, 75, 95]).round(3).tolist() if rs else None
    rep = dict(ident_dist=idq(R_), r_dyn_median=float(np.median([r['main']['r_dyn'] for r in R_])), r_dyn_iqr=np.percentile([r['main']['r_dyn'] for r in R_], [25, 75]).round(3).tolist(),
               ident_dist_dropped_by_inter=idq([r for r in R_ if r['inter'].get('r_dyn') is None]), shared_edges_median=float(np.median([r['main']['n_shared'] for r in R_])),
               tost_010=bool(main.get('lo', -1) > -0.10 and main.get('hi', 1) < 0.10), degen_frac_median=float(np.median([r['main']['degen_frac'] for r in R_])))
    out = dict(written=time.strftime('%F %T'), verdict=v, main=main, joint_gate=dict(pair_set='교집합 ≥ 30 변인 쌍', main=gate_main, inter=gate_inter),
               freedman_lane=fl, sensitivity=sens, dropouts=drop, report=rep,
               pending_sensitivity=['메도이드·결정 구조 통제', 'Q 0.7·0.8', '반복 나눔 신뢰도 보정', '로그 퍼짐 비 반응', '281 쌍 판'],
               rule='prereg 1b: 편 Spearman(r_dyn, 동일성 | 공유 핵심부 거리), 쌍 무게 1/상과 쌍 수, 교차 군집 부트스트랩 2,000. 효과 있음 = 하한 > 0 이고 상한 ≥ 0.15 이고 같은 쌍 집합의 교집합 판 하한 > 0')
    jdump(out, RES + '/r38_p3.json', bk=True); log('P3', v, {k: round(main[k], 4) for k in ('est', 'lo', 'hi') if k in main}, '교집합', {k: round(gate_inter[k], 4) for k in ('est', 'lo', 'hi') if k in gate_inter}); log('P3_STATS_DONE')


ALL = ['pairs', 'dom', 'pibar', 'pair', 'stats']
if __name__ == '__main__':
    for part in (ALL if PARTS == ['all'] else PARTS):
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('P3_DONE', ','.join(PARTS))
