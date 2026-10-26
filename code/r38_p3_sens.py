"""(prereg round36 1b, 1차 P3 의 민감도 다섯 가운데 넷; 281 쌍 판은 r38_b6.py) r38_p3.py 의 함수를 불러 쓴다(그 파일은 바꾸지 않음).
대응은 주 판의 π̄ 공유 변(퇴화 변 뺀 것)으로 고정하고, 바뀌는 것만 바꿔 낸다:
  (1) 통제: 메도이드 프레임의 π̄ 짝 몫, 결정(시작) 구조의 π̄ 짝 몫 (h5 pdbProteinAtoms 좌표), 메도이드 정렬 자체의 짝 몫(메도이드끼리 128 출발 정렬, GPU).
      메도이드 = 앙상블 안 d_fix² 합이 가장 작은 프레임.
  (2) Q 0.7·0.8: 같은 규칙(뒤 절반, 반복마다 ≥ 40 이면 최대 100, 반복 ≥ 2)으로 문턱만 바꾼 앙상블에서 공유 변 GPA 윤곽(r_dyn)과 사영 평균 점 통제.
      공유 변 = 주 판 공유 변에서 그 앙상블의 퇴화 변을 더 뺀 것(≥ 30).
  (3) 반복 나눔 신뢰도 보정: 도메인마다 반복을 정렬해 앞 ⌈k/2⌉ 대 나머지(5 개면 3 대 2)로 나눠 같은 공유 변 GPA 윤곽을 내고,
      r_half = 두 윤곽의 Spearman, 신뢰도 = 2 r_half / (1 + r_half) (Spearman–Brown). r_dyn_corr = r_dyn / √(신뢰도_A 신뢰도_B), [−1, 1] 로 자름; 신뢰도 ≤ 0 이면 뺌.
  (4) 반응 = 로그 퍼짐 비: −|log(Σ 퍼짐_A / Σ 퍼짐_B)| (공유 변 GPA 윤곽; 클수록 닮음, r_dyn 과 같은 방향).
분석은 r38_p3.analyse 그대로(가중 편 Pearson, 교차 군집 부트스트랩 2,000).
부분: dom2 (CPU: 메도이드·시작 구조·Q 0.7·0.8 앙상블) → medalign (GPU: 메도이드 정렬) → pair2 (CPU) → stats2 (CPU)
사용: python r38_p3_sens.py <ROOT> <부분들 쉼표 | all> [작업자 수=3]
출력: derived/r38/p3/{dom2/<d>.npz, medalign.jsonl, pairs_sens.jsonl}, results/r38_p3_sens.json"""
import os, sys
ROOT = sys.argv[1]; PARTS2 = sys.argv[2].split(',') if len(sys.argv) > 2 else ['all']; NW2 = int(sys.argv[3]) if len(sys.argv) > 3 else 3
_src = open(ROOT + '/code/r38_p3.py').read(); _src = _src[:_src.index("ALL = ['pairs'")]
sys.argv = ['r38_p3', ROOT, 'none', str(NW2)]; exec(_src)
OUT2 = OUT + '/dom2'; os.makedirs(OUT2, exist_ok=True)
QV = (0.7, 0.8)


def start_frame(d):
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); txt = h[d]['pdbProteinAtoms'][()].decode(); h.close()
    lines = [l for l in txt.splitlines() if l.startswith(('ATOM', 'HETATM'))]; smd, imd = md_atoms(txt)
    X = np.array([[[float(lines[k][30:38]), float(lines[k][38:46]), float(lines[k][46:54])] for k in row] for row in imd])   # (L, 3, 3) CA, C, O
    return frames(X[None])[0], smd


def build_q(d, q):
    """dom_build 와 같은 규칙, 문턱만 q."""
    global QT
    old = QT; QT = q
    try: return dom_build(d)
    finally: QT = old


def _dom2_worker(d):
    fo = f'{OUT2}/{d}.npz'
    if os.path.exists(fo): return d, 'skip'
    try:
        z = dom_get(d)
        if z is None: np.savez(fo, skip='주 앙상블 없음'); return d, 'skip-main'
        A = z['A']; Dm = np.asarray(fs.fixed_sq_distance(A, A)); med = A[int(Dm.sum(1).argmin())]
        st, smd = start_frame(d); out = dict(medoid=med, start=st)
        for q in QV:
            r, info = build_q(d, q); tag = f'q{int(q * 100):03d}'
            if r is None: out[tag + '_skip'] = np.array(1); continue
            np.savez(f'{CACHE}/{d}.{tag}.npz', **r); out.update({tag + '_Pm': r['Pm'], tag + '_degen': r['degen'], tag + '_n': np.array(len(r['A']))})
        np.savez(fo, **out); return d, 'ok'
    except Exception as e:
        return d, 'ERR ' + repr(e)[:150]


def part_dom2():
    doms = json.load(open(OUT + '/pairs.json'))['domains']; from multiprocessing import Pool
    with Pool(NW2) as p:
        for k, (d, st) in enumerate(p.imap_unordered(_dom2_worker, doms)):
            if st not in ('skip',) and (k % 100 == 0 or st.startswith('ERR')): log('dom2', k, d, st)
    log('DOM2_DONE')


def q_get(d, tag):
    f = f'{CACHE}/{d}.{tag}.npz'
    if os.path.exists(f): z = np.load(f); return {k: z[k] for k in z.files}
    r, info = build_q(d, float(tag[1:]) / 100)
    if r is None: return None
    np.savez(f, **r); return r


def part_medalign():
    """메도이드끼리 정렬(실제 구조, 항등 + 무작위 127 출발) — pibar 와 같은 묶음 계산. 짝 몫 = Σ 짝 ‖x − U y‖² / 짝 수."""
    dev = 'cuda' if torch.cuda.is_available() else 'cpu'; P = json.load(open(OUT + '/pairs.json'))['pairs']; fo = OUT + '/medalign.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    M = {}
    for d in sorted({p[0] for p in P} | {p[1] for p in P}):
        f = f'{OUT2}/{d}.npz'
        if os.path.exists(f):
            z = np.load(f)
            if 'medoid' in z.files: M[d] = z['medoid'].astype(np.float64)
    todo = sorted([p for p in P if (p[0], p[1]) not in done and p[0] in M and p[1] in M], key=lambda p: (len(M[p[0]]) // 16, len(M[p[1]]) // 16)); log('medalign', len(todo), dev)
    BS = 48
    for b0 in range(0, len(todo), BS):
        B = todo[b0:b0 + BS]; n = max(len(M[p[0]]) for p in B); m = max(len(M[p[1]]) for p in B)
        x = fs.pad_stack([M[p[0]] for p in B], n).to(dev, torch.float64); y = fs.pad_stack([M[p[1]] for p in B], m).to(dev, torch.float64)
        al = fs.align(x, y, n_random=127, iters=10, seed=b0); Ms = al.match.cpu().numpy(); Us = al.rotation.cpu().numpy(); Sv = al.score.cpu().numpy()
        with open(fo, 'a') as fh:
            for q, p in enumerate(B):
                na, nb = len(M[p[0]]), len(M[p[1]]); ii, jj = np.nonzero(Ms[q][:na, :nb]); a_, b_ = M[p[0]][ii], M[p[1]][jj]
                share = float(((a_ - np.einsum('xy,iyz->ixz', Us[q], b_)) ** 2).sum() / max(len(ii), 1))
                fh.write(json.dumps(dict(ab=[p[0], p[1]], S=float(Sv[q]), n_match=int(len(ii)), share=share)) + '\n')
    log('MEDALIGN_DONE')


def halves_profile(F, rep, idx):
    u = sorted(set(rep.tolist())); k = (len(u) + 1) // 2; g1 = np.isin(rep, u[:k]); g2 = ~g1
    if g1.sum() < 10 or g2.sum() < 10: return None
    return spearman(gpa_profile(F[g1], idx), gpa_profile(F[g2], idx))


def sb(r): return None if r is None or r <= 0 else 2 * r / (1 + r)


def part_pair2():
    P = json.load(open(OUT + '/pairs.json'))['pairs']; fo = OUT + '/pairs_sens.jsonl'; done = set()
    if os.path.exists(fo):
        for l in open(fo):
            try: done.add(tuple(json.loads(l)['ab']))
            except Exception: pass
    PB = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(OUT + '/pibar.jsonl')}
    MA = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(OUT + '/medalign.jsonl')} if os.path.exists(OUT + '/medalign.jsonl') else {}
    todo = sorted([p for p in P if (p[0], p[1]) in PB and (p[0], p[1]) not in done]); log('pair2', len(todo))
    cache = {}
    def get(d):
        if d not in cache:
            if len(cache) > 30: cache.pop(next(iter(cache)))
            s_, _ = dom_small(d); z = dom_get(d); f = f'{OUT2}/{d}.npz'; e = dict(np.load(f)) if os.path.exists(f) else {}
            cache[d] = (s_, z, e, {t: q_get(d, t) for t in ('q070', 'q080') if f'{t}_Pm' in e})
        return cache[d]
    t0 = time.time()
    with open(fo, 'a') as fh:
        for k, p in enumerate(todo):
            try:
                (Da, za, ea, qa), (Db, zb, eb, qb) = get(p[0]), get(p[1]); rec = PB[(p[0], p[1])]
                if Da is None or Db is None or za is None or zb is None: continue
                prs = np.array(rec['best'], int).reshape(-1, 2); keep = ~(Da['degen'][prs[:, 0]] | Db['degen'][prs[:, 1]]) if len(prs) else np.zeros(0, bool); sh = prs[keep]
                o = dict(ab=[p[0], p[1]], H=p[2], n_shared=int(len(sh)))
                if len(sh) >= MINSH:
                    va, vb = gpa_profile(za['A'], sh[:, 0]), gpa_profile(zb['A'], sh[:, 1]); r0 = spearman(va, vb)
                    o.update(r_dyn=r0, ident=float(np.mean([Da['seq'][i] == Db['seq'][j] for i, j in sh])), control=rot_share(Da['Pm'][sh[:, 0]], Db['Pm'][sh[:, 1]]),
                             logspread=-abs(float(np.log(va.sum() / vb.sum()))))
                    if 'medoid' in ea and 'medoid' in eb: o['control_medoid'] = rot_share(ea['medoid'][sh[:, 0]].astype(np.float64), eb['medoid'][sh[:, 1]].astype(np.float64))
                    if 'start' in ea and 'start' in eb: o['control_start'] = rot_share(ea['start'][sh[:, 0]], eb['start'][sh[:, 1]])
                    if (p[0], p[1]) in MA: o['control_medalign'] = MA[(p[0], p[1])]['share']
                    ha, hb = sb(halves_profile(za['A'], za['rep'], sh[:, 0])), sb(halves_profile(zb['A'], zb['rep'], sh[:, 1]))
                    o.update(rel_a=ha, rel_b=hb, r_dyn_corr=None if ha is None or hb is None else float(np.clip(r0 / np.sqrt(ha * hb), -1, 1)))
                    for t in ('q070', 'q080'):
                        if t in qa and t in qb and qa[t] is not None and qb[t] is not None:
                            kq = ~(ea[f'{t}_degen'][sh[:, 0]] | eb[f'{t}_degen'][sh[:, 1]]); s2 = sh[kq]
                            if len(s2) >= MINSH:
                                o[f'r_dyn_{t}'] = spearman(gpa_profile(qa[t]['A'], s2[:, 0]), gpa_profile(qb[t]['A'], s2[:, 1]))
                                o[f'control_{t}'] = rot_share(ea[f'{t}_Pm'][s2[:, 0]], eb[f'{t}_Pm'][s2[:, 1]]); o[f'ident_{t}'] = float(np.mean([Da['seq'][i] == Db['seq'][j] for i, j in s2]))
                fh.write(json.dumps(o) + '\n')
            except Exception as e:
                log('ERR pair2', p[:2], repr(e)[:150])
            if k % 500 == 0: log('pair2', k, '/', len(todo), round(time.time() - t0), 's')
    log('PAIR2_DONE')


def part_stats2():
    rows = [json.loads(l) for l in open(OUT + '/pairs_sens.jsonl')]
    W = lambda rs: [dict(ab=r['ab'], H=r['H'], main=r) for r in rs]
    ok = lambda key: [r for r in rows if r.get(key) is not None and r.get('r_dyn') is not None]
    out = dict(written=time.strftime('%F %T'),
               control_medoid=analyse(W(ok('control_medoid')), z=('control_medoid',), seed=11),
               control_start=analyse(W(ok('control_start')), z=('control_start',), seed=12),
               control_medalign=analyse(W(ok('control_medalign')), z=('control_medalign',), seed=13),
               reliability_corrected=analyse(W([r for r in rows if r.get('r_dyn_corr') is not None]), y='r_dyn_corr', seed=14),
               reliability=dict(rel_median=float(np.median([x for r in rows for x in (r.get('rel_a'), r.get('rel_b')) if x is not None])) if rows else None),
               logspread_response=analyse(W(ok('logspread')), y='logspread', seed=15))
    for t in ('q070', 'q080'):
        rr = [dict(ab=r['ab'], H=r['H'], main=dict(r_dyn=r[f'r_dyn_{t}'], ident=r[f'ident_{t}'], control=r[f'control_{t}'])) for r in rows if r.get(f'r_dyn_{t}') is not None]
        out[t] = analyse(rr, seed=16 if t == 'q070' else 17)
    out['rule'] = 'r38_p3_sens 머리말. 대응은 주 판 π̄ 공유 변으로 고정. 판정은 주 판(r38_p3.json)이 하고 이것은 민감도'
    jdump(out, RES + '/r38_p3_sens.json', bk=True); log('P3_SENS', {k: (round(v['est'], 3) if isinstance(v, dict) and 'est' in v else None) for k, v in out.items()}); log('P3_SENS_DONE')


ALL2 = ['dom2', 'medalign', 'pair2', 'stats2']
if __name__ == '__main__':
    for part in (ALL2 if PARTS2 == ['all'] else PARTS2):
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('P3_SENS_PARTS_DONE', ','.join(PARTS2))
