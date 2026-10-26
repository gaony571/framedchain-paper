"""(prereg round36 6 절 최소판; 범위 줄임 2026-10-20 "연산 시간 표 + 인증 비율; 인증은 정적 축 BnB 500 으로 갈음"; 등록 = Drive JSON round36_sec6_2026-10-20, 돌리기 전)
부분: time (GPU, 대기열 끝 쪽) → stats (CPU; 정적 축·줄기 336 D 뒤).
time: 구조 쌍 시간. 정적 축 후보 쌍(results/r37_trunk_candidates.tsv 의 level T·A·C·X)에서 두 메도이드가 있는 쌍을 긴 쪽 변 수 구간(≤ 100, 101–300, 301–500)마다 30 (default_rng(20261078)).
  도구: 기본 판 휴리스틱(fs.align, 항등 + 무작위 127 출발, iters 10; GPU 가 되면 GPU, 동기화한 벽시계), BnB 인증(fs.bnb rel_tol 0.005; 인증 여부·시간),
  TM-align 재구현(stage0.rulers.tm_align, 메도이드 Cα, d₀ = 긴 쪽; Python, CPU), C++ TM-align(PyPI tmtools 의 tm_align, 메모리 안 호출),
  US-align(GitHub pylelab/USalign 을 g++ -O3 로 컴파일; Cα 만 쓴 PDB 두 개를 넘기는 프로세스 한 번 = 파일 쓰기 제외, 프로세스 시작 포함). 첫 쌍을 한 번 더 돌려 데우기로 뺀다.
  (사용자 허락 2026-10-20 UTC: tmtools·US-align 내려받기.) 도구가 없으면 그 열은 "설치 안 함".
stats: (1) 인증 비율 = 정적 축 BnB 500 (derived/r38/static/bnb.jsonl) 에서 인증이고 ≤ 60 s 인 비율 (예측 ≥ 95%);
  (2) 휴리스틱 정확도 = (S_ub − S_heur)/S_ub ≤ 0.002 인 비율 (예측 ≥ 95%; D 의 상대 오차 상한 rel_err_ub 도 함께 적음);
  (3) 시간 표 = 구간 × 도구의 중앙값·90% 점; 앙상블 쌍 시간(서술) = 줄기 336 D 의 도메인당 초(모든 칸)와 R10a′ 쌍당 초를 길이 구간별 중앙값으로.
사용: python r39_sec6.py <ROOT> time|stats 출력: derived/r39/sec6_time.jsonl, results/r39_sec6.json"""
import os, sys
ROOT = sys.argv[1]; PARTS6 = sys.argv[2].split(',')
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)
from stage0.rulers import tm_align, d0_tm
import pandas as pd
import subprocess
OUT6 = ROOT + '/derived/r39'; os.makedirs(OUT6, exist_ok=True); BINS6 = [('≤100', 0, 100), ('101-300', 101, 300), ('301-500', 301, 500)]
try:
    import tmtools; HAVE_TMT = True
except Exception: HAVE_TMT = False
USALIGN = '/content/usalign/USalign' if os.path.exists('/content/usalign/USalign') else None; TMPS6 = '/content/sec6_tmp'; os.makedirs(TMPS6, exist_ok=True)


def ca_pdb(P, f):
    with open(f, 'w') as fh:
        sp = lambda k: ' ' * k # PDB 칸 맞춤(공백 수를 글자로 쓰지 않음)
        for i, (x, y, z) in enumerate(P): fh.write('ATOM' + sp(2) + f'{i + 1:5d}' + sp(2) + 'CA' + sp(2) + f'ALA A{i + 1:4d}' + sp(4) + f'{x:8.3f}{y:8.3f}{z:8.3f}' + sp(2) + '1.00' + sp(2) + '0.00' + sp(11) + 'C\n')
        fh.write('END\n')


def medoid_of(d):
    f2 = f'{ROOT}/derived/r38/p3/dom2/{d}.npz'
    if os.path.exists(f2):
        z = np.load(f2)
        if 'medoid' in z.files: return z['medoid'].astype(np.float64)
    f = f'{ROOT}/derived/r38/static/medoid/{d}.npy'
    return np.load(f).astype(np.float64) if os.path.exists(f) else None


def sync():
    if DEV == 'cuda': torch.cuda.synchronize()


def time_pair(x, y):
    o = {}
    X, Y = torch.as_tensor(x).to(DEV, DT), torch.as_tensor(y).to(DEV, DT)
    sync(); t0 = time.time(); al = fs.align(X, Y); sync(); o['heur_sec'] = time.time() - t0; o['S_heur'] = float(al.score)
    t0 = time.time(); c = fs.bnb(x, y, rel_tol=0.005, device=DEV); o['bnb_sec'] = time.time() - t0; o['bnb_sec_rep'] = float(c['sec']); o['certified'] = bool(c['certified'])
    PA, PB = np.asarray(fs.ca_trace(x)), np.asarray(fs.ca_trace(y)); t0 = time.time(); r = tm_align(PA, PB, d0_tm(max(len(PA), len(PB)))); o['tmalign_py_sec'] = time.time() - t0; o['TM'] = float(r['S'])
    if HAVE_TMT:
        PA_, PB_ = np.ascontiguousarray(PA, dtype=np.float64), np.ascontiguousarray(PB, dtype=np.float64); t0 = time.time(); rr = tmtools.tm_align(PA_, PB_, 'A' * len(PA_), 'A' * len(PB_)); o['tmalign_cpp_sec'] = time.time() - t0
        o['TM_cpp'] = float(max(rr.tm_norm_chain1, rr.tm_norm_chain2))
    if USALIGN:
        fa, fb = TMPS6 + '/a.pdb', TMPS6 + '/b.pdb'; ca_pdb(PA, fa); ca_pdb(PB, fb); t0 = time.time(); p = subprocess.run([USALIGN, fa, fb, '-outfmt', '2'], capture_output=True, text=True); o['usalign_sec'] = time.time() - t0
        o['usalign_ok'] = p.returncode == 0
    return o


def part_time():
    C = pd.read_csv(RES + '/r37_trunk_candidates.tsv', sep='\t'); C = C[C.level.isin(['T', 'A', 'C', 'X'])].sort_values(['level', 'a', 'b']).reset_index(drop=True)
    M = {}; rows = []
    for r in C.itertuples():
        for d in (r.a, r.b):
            if d not in M: M[d] = medoid_of(d)
        if M[r.a] is not None and M[r.b] is not None: rows.append((r.a, r.b, max(len(M[r.a]), len(M[r.b]))))
    rng = np.random.default_rng(20261078); pick = []
    for nm, lo, hi in BINS6:
        cand = [x for x in rows if lo <= x[2] <= hi]; sel = [cand[i] for i in sorted(rng.choice(len(cand), min(30, len(cand)), replace=False))] if cand else []
        pick += [(nm, *x) for x in sel]; log('sec6 구간', nm, '후보', len(cand), '고름', len(sel))
    fo = OUT6 + '/sec6_time.jsonl'; done = set()
    if os.path.exists(fo): done = {tuple(json.loads(l)['ab']) for l in open(fo)}
    if pick: time_pair(M[pick[0][1]], M[pick[0][2]])
    for nm, a, b, n in pick:
        if (a, b) in done: continue
        try: o = dict(bin=nm, ab=[a, b], n=len(M[a]), m=len(M[b]), dev=DEV, **time_pair(M[a], M[b]))
        except Exception as e: o = dict(bin=nm, ab=[a, b], err=repr(e)[:200])
        with open(fo, 'a') as fh: fh.write(json.dumps(o) + '\n')
        log('sec6', nm, a, b, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in o.items() if k.endswith('sec') or k == 'certified'})
        if DEV == 'cuda': torch.cuda.empty_cache()
    log('SEC6_TIME_DONE')


def part_stats():
    bn = [json.loads(l) for l in open(ROOT + '/derived/r38/static/bnb.jsonl')]; ok = [r for r in bn if 'err' not in r]
    cert = [r['certified'] and r['sec'] <= 60 for r in ok]; gap = [(r['S_ub'] - r['S_heur']) / r['S_ub'] for r in ok if r['S_ub'] > 0]; rel = [r['rel_err_ub'] for r in ok if r['rel_err_ub'] is not None]
    out = dict(written=time.strftime('%F %T'), bnb=dict(n=len(bn), n_err=len(bn) - len(ok), frac_certified_le60s=float(np.mean(cert)) if cert else None, frac_certified_any=float(np.mean([r['certified'] for r in ok])) if ok else None,
                                                      frac_sim_gap_le_0p002=float(np.mean([g <= 0.002 for g in gap])) if gap else None, sim_gap_median=float(np.median(gap)) if gap else None, sim_gap_max=float(max(gap)) if gap else None,
                                                      rel_err_ub_median=float(np.median(rel)) if rel else None, rel_err_ub_max=float(max(rel)) if rel else None, sec_median=float(np.median([r['sec'] for r in ok])) if ok else None,
                                                      sec_p90=float(np.percentile([r['sec'] for r in ok], 90)) if ok else None))
    b = out['bnb']; out['prediction'] = dict(cert_ge95=(b['frac_certified_le60s'] or 0) >= 0.95, heur_ge95=(b['frac_sim_gap_le_0p002'] or 0) >= 0.95)
    # (고침 1, 2026-10-20 19:1x UTC, 사후 서술) 등록 정의는 S_ub = S_lb + 허용 폭(0.005 sim_long) 이라 정의상 0.002 안일 수 없다. 같은 단위(sim_long = S / (0.75 × 긴 쪽 변 수))로 BnB 가 찾은 최선값과 휴리스틱의 차를 함께 낸다
    sl = [(r['S_lb'] - r['S_heur']) / (0.75 * max(r['n'], r['m'])) for r in ok]; out['bnb_posthoc'] = dict(frac_heur_within_0p002_simlong_of_S_lb=float(np.mean([g <= 0.002 for g in sl])), max_shortfall_simlong=float(max(sl)),
        n_bnb_better=int(sum(g > 1e-12 for g in sl)), note='사후 서술: 등록 정의(S_ub − S_heur)/S_ub ≤ 0.002 는 S_ub 에 허용 폭 0.005 sim_long 이 늘 더해져 설 수 없음(등록 오류). 판정 값은 그대로 보고')
    T = [json.loads(l) for l in open(OUT6 + '/sec6_time.jsonl')] if os.path.exists(OUT6 + '/sec6_time.jsonl') else []; T = [t for t in T if 'err' not in t]; tab = {}
    for nm, _, _ in BINS6:
        v = [t for t in T if t['bin'] == nm]
        ks6 = [k for k in ('heur_sec', 'bnb_sec', 'tmalign_py_sec', 'tmalign_cpp_sec', 'usalign_sec') if all(k in t for t in v)]
        if v: tab[nm] = dict(n=len(v), **{f'{k}_median': float(np.median([t[k] for t in v])) for k in ks6}, **{f'{k}_p90': float(np.percentile([t[k] for t in v], 90)) for k in ks6},
                             frac_certified=float(np.mean([t['certified'] for t in v])), missing_tools=[k for k in ('tmalign_cpp_sec', 'usalign_sec') if k not in ks6])
    out['time_table'] = tab
    E = [json.load(open(f)) for f in glob.glob(ROOT + '/derived/r39/D/*.json')]; ens = {}
    for nm, lo, hi in BINS6:
        v = [j['sec'] for j in E if lo <= j.get('ne', 0) <= hi and 'sec' in j]
        if v: ens[nm] = dict(n=len(v), sec_median=float(np.median(v)), note='줄기 336 D: 도메인당 모든 칸(온도 4 × 창 2 + R0 나눔), 실행마다 N = 100, 실행 쌍마다 프레임 쌍 2,000')
    R10 = [json.load(open(f)) for f in glob.glob(ROOT + '/derived/r38/r10a/*.json')]; ens['R10a′ 쌍당'] = dict(n=len(R10), sec_median=float(np.median([r['sec'] for r in R10 if 'sec' in r])) if R10 else None)
    out['ensemble_time'] = ens
    out['rule'] = 'prereg 6 절 최소판: 인증 ≥ 95% (인증이고 ≤ 60 s), 휴리스틱 (S_ub − S_heur)/S_ub ≤ 0.002 인 쌍 ≥ 95%; 연산량은 보고만. 인증 표본 = 정적 축 BnB 500 (7g).'
    jdump(out, RES + '/r39_sec6.json', bk=True); log('SEC6', out['bnb'], out['prediction'], out['time_table']); log('SEC6_STATS_DONE')


if __name__ == '__main__':
    for part in PARTS6:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
