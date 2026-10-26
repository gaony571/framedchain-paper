"""(round41 g5 측정; Drive JSON round41_2026-10-21 g5_design 07:33:45) 같은 장비 연산량 표. 입력은 r41_compute_prep.py 가 만든 /content/ct (로컬).
부분 (대기열이 차례로 부름; 측정 중 다른 계산·Drive 입출력 없음, 결과는 /content/ct/out 에 쓰고 sum 에서 Drive 로 복사):
  cpu     (taskset 코어 하나, 스레드 1) 500 쌍마다: C++ TM-align(tmtools, 메모리 안), US-align(프로세스 하나, 시작 포함), D 1 출발(fs.align 항등 출발, iters 10, CPU);
          D 128 출발(CPU) 은 무작위 50 쌍. 첫 5 쌍은 데우기로 한 번 먼저 돌려 뺀다. 무작위 50 쌍을 한 번 더 재어 재현성.
  foldseek createdb (스레드 1; 570 구조와 200 구조), 전체 대 전체 200 (--exhaustive-search 1 -e inf; 3Di+AA 와 --alignment-type 1; 스레드 1 은 코어 하나, 12 는 기계 전체),
          단일 쌍 easy-search 20 쌍(끝에서 끝까지).
  usava   US-align 전체 대 전체 200 (-dir, 한 프로세스, 코어 하나).
  gpu     D 128 출발 묶음(정적 축과 같은 묶음 규칙)·D 1 출발 묶음 500 쌍, BnB rel_tol 0.005 500 쌍(차례로), 전체 대 전체 D 1 출발 200 (19,900 쌍)·D 128 출발 100 (4,950 쌍). 최대 GPU 메모리.
  ens     앙상블 칸 하나: 도메인 셋(ne ~80·150·250)의 320 K·379 K 실행 쌍 하나, 프레임 쌍 2,000 의 줄인 계산(MeanStart + red_pairs) 3 회.
  sum     요약 → results/r41_compute.json, 원자료 사본 derived/r41/compute/.
사용: python r41_compute.py <ROOT> <부분>"""
import os, sys
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '1'
ROOT = sys.argv[1]; PART = sys.argv[2]
if PART in ('cpu', 'foldseek', 'usava', 'sum'): os.environ['CUDA_VISIBLE_DEVICES'] = ''
import json, time, glob, shutil, subprocess, resource
import numpy as np
T = '/content/ct'; OUTC = T + '/out'; os.makedirs(OUTC, exist_ok=True); FS = T + '/foldseek/bin/foldseek'; US = T + '/USalign/USalign'; SEED = 20261089


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def loadx(d): return np.load(f'{T}/x/{d}.npy')


def wr(name, rows):
    with open(f'{OUTC}/{name}.jsonl', 'w') as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')


def part_cpu():
    import torch; torch.set_num_threads(1); sys.path[:0] = [ROOT + '/framespace_pkg']; import framespace as fs; from tmtools import tm_align
    from framespace.frames import ca_trace
    P = json.load(open(f'{T}/pairs.json')); ds = sorted({p['a'] for p in P} | {p['b'] for p in P}); X = {d: loadx(d) for d in ds}
    CA = {d: np.ascontiguousarray(ca_trace(X[d])) for d in ds}; XT = {d: torch.as_tensor(X[d], dtype=torch.float64) for d in ds}
    SQ = {}
    for d in ds:
        s = [l[17:20] for l in open(f'{T}/pdb/{d}.pdb') if l.startswith('ATOM') and l[12:16] == ' CA ']; SQ[d] = 'A' * len(s)
    def one(p, which):
        a, b = p['a'], p['b']; o = {}
        if 'tm' in which: t0 = time.perf_counter(); r = tm_align(CA[a], CA[b], SQ[a], SQ[b]); o['tmalign_sec'] = time.perf_counter() - t0; o['tm'] = float(min(r.tm_norm_chain1, r.tm_norm_chain2))
        if 'us' in which: t0 = time.perf_counter(); q = subprocess.run([US, f'{T}/pdb/{a}.pdb', f'{T}/pdb/{b}.pdb', '-outfmt', '2'], capture_output=True, text=True); o['usalign_sec'] = time.perf_counter() - t0; o['us_ok'] = q.returncode == 0
        if 'd1' in which: t0 = time.perf_counter(); al = fs.align(XT[a], XT[b], n_random=0, iters=10); o['d1_cpu_sec'] = time.perf_counter() - t0; o['d1_S'] = float(al.score)
        if 'd128' in which: t0 = time.perf_counter(); al = fs.align(XT[a], XT[b], n_random=127, iters=10, seed=0); o['d128_cpu_sec'] = time.perf_counter() - t0; o['d128_S'] = float(al.score)
        return o
    for p in P[:5]: one(p, ('tm', 'us', 'd1'))                      # 데우기(버림)
    rows = []; t0 = time.time()
    for k, p in enumerate(P):
        rows.append(dict(a=p['a'], b=p['b'], n=p['n'], m=p['m'], **one(p, ('tm', 'us', 'd1'))))
        if k % 100 == 0: log('cpu', k, round(time.time() - t0), 's')
    rss_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss; rss_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    rng = np.random.default_rng(SEED); sub = sorted(rng.choice(len(P), 50, replace=False).tolist())
    rep = [dict(i=i, **one(P[i], ('tm', 'us', 'd1'))) for i in sub]
    d128 = [dict(i=i, a=P[i]['a'], b=P[i]['b'], n=P[i]['n'], m=P[i]['m'], **one(P[i], ('d128',))) for i in sub]
    wr('cpu', rows); wr('cpu_repeat', rep); wr('cpu_d128', d128)
    json.dump(dict(maxrss_self_kb=rss_self, maxrss_child_kb=rss_child, affinity=sorted(os.sched_getaffinity(0)), torch_threads=torch.get_num_threads(), sub=sub), open(f'{OUTC}/cpu_meta.json', 'w'))
    log('CT_CPU_DONE', len(rows))


def runt(cmd, pin=None):
    full = (['taskset', '-c', str(pin)] if pin is not None else []) + cmd; t0 = time.perf_counter(); q = subprocess.run(full, capture_output=True, text=True); dt = time.perf_counter() - t0
    return dt, q.returncode, q.stderr[-400:]


def part_foldseek():
    W = T + '/fsw'; shutil.rmtree(W, ignore_errors=True); os.makedirs(W); o = dict(version=subprocess.run([FS, 'version'], capture_output=True, text=True).stdout.strip())
    os.makedirs(W + '/p200'); ava = open(f'{T}/ava200.txt').read().split()
    for d in ava: os.symlink(f'{T}/pdb/{d}.pdb', f'{W}/p200/{d}.pdb')
    o['createdb_570'] = runt([FS, 'createdb', f'{T}/pdb', f'{W}/db570', '--threads', '1'], pin=1)[:2]
    o['createdb_200'] = runt([FS, 'createdb', f'{W}/p200', f'{W}/db200', '--threads', '1'], pin=1)[:2]
    o['search'] = {}
    for at, nm in (('2', '3Di+AA'), ('1', 'TM-align')):
        for th in (1, 12):
            tmp = f'{W}/tmp_{at}_{th}'; aln = f'{W}/aln_{at}_{th}'
            ev = 'inf'; dt, rc, err = runt([FS, 'search', f'{W}/db200', f'{W}/db200', aln, tmp, '--exhaustive-search', '1', '-e', ev, '--alignment-type', at, '--threads', str(th)], pin=1 if th == 1 else None)
            if rc:
                for f in glob.glob(aln + '*'): os.remove(f)
                ev = '1e10'; dt, rc, err = runt([FS, 'search', f'{W}/db200', f'{W}/db200', aln, tmp, '--exhaustive-search', '1', '-e', ev, '--alignment-type', at, '--threads', str(th)], pin=1 if th == 1 else None)
            subprocess.run([FS, 'convertalis', f'{W}/db200', f'{W}/db200', aln, aln + '.m8', '--threads', '12'], capture_output=True)
            nh = sum(1 for _ in open(aln + '.m8')) if os.path.exists(aln + '.m8') else None
            o['search'][f'{nm}|threads{th}'] = dict(sec=dt, rc=rc, evalue=ev, n_hits=nh, n_pairs_all=len(ava) ** 2, err=err if rc else None); log('foldseek', nm, th, round(dt, 1), 's', nh, ev)
    P = json.load(open(f'{T}/pairs.json')); rng = np.random.default_rng(SEED + 1); sub = sorted(rng.choice(len(P), 20, replace=False).tolist()); e2e = []
    for i in sub:
        p = P[i]; tmp = f'{W}/e2e_tmp'; shutil.rmtree(tmp, ignore_errors=True)
        dt, rc, _ = runt([FS, 'easy-search', f'{T}/pdb/{p["a"]}.pdb', f'{T}/pdb/{p["b"]}.pdb', f'{W}/e2e.m8', tmp, '--exhaustive-search', '1', '-e', 'inf', '--threads', '1'], pin=1)
        e2e.append(dict(i=i, n=p['n'], m=p['m'], sec=dt, rc=rc))
    o['easy_search_pair'] = e2e; o['maxrss_child_kb'] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    json.dump(o, open(f'{OUTC}/foldseek.json', 'w'), indent=1); log('CT_FOLDSEEK_DONE')


def part_usava():
    ava = open(f'{T}/ava200.txt').read().split(); open(f'{T}/ava200_files.txt', 'w').write('\n'.join(f'{d}.pdb' for d in ava) + '\n')
    t0 = time.perf_counter(); q = subprocess.run(['taskset', '-c', '1', US, '-dir', f'{T}/pdb/', f'{T}/ava200_files.txt', '-outfmt', '2'], capture_output=True, text=True); dt = time.perf_counter() - t0
    nl = sum(1 for l in q.stdout.splitlines() if l and not l.startswith('#'))
    json.dump(dict(sec=dt, rc=q.returncode, n_lines=nl, n_pairs_expected=len(ava) * (len(ava) - 1) // 2, maxrss_child_kb=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss), open(f'{OUTC}/usava.json', 'w'), indent=1)
    log('CT_USAVA_DONE', round(dt, 1), nl)


def part_gpu():
    import torch; sys.path[:0] = [ROOT + '/framespace_pkg']; import framespace as fs
    dev = 'cuda'; P = json.load(open(f'{T}/pairs.json')); ds = sorted({p['a'] for p in P} | {p['b'] for p in P}); X = {d: loadx(d) for d in ds}
    cap = torch.cuda.get_device_properties(0).total_memory / (128 * 8 * 16); sync = torch.cuda.synchronize
    def batches(prs):
        todo = sorted(prs, key=lambda p: (len(X[p[0]]) // 16, len(X[p[1]]) // 16)); b0 = 0; out = []
        while b0 < len(todo):
            k_ = 1
            while k_ < 48 and b0 + k_ < len(todo) and (k_ + 1) * max(len(X[p[0]]) for p in todo[b0:b0 + k_ + 1]) * max(len(X[p[1]]) for p in todo[b0:b0 + k_ + 1]) <= cap: k_ += 1
            out.append(todo[b0:b0 + k_]); b0 += k_
        return out
    def run_batches(prs, nr, tag):
        B = batches(prs); rows = []; torch.cuda.reset_peak_memory_stats(); sync(); T0 = time.perf_counter()
        for bi, Bt in enumerate(B):
            n = max(len(X[p[0]]) for p in Bt); m = max(len(X[p[1]]) for p in Bt)
            xx = fs.pad_stack([X[p[0]] for p in Bt], n).to(dev, torch.float64); yy = fs.pad_stack([X[p[1]] for p in Bt], m).to(dev, torch.float64); sync(); t0 = time.perf_counter()
            al = fs.align(xx, yy, n_random=nr, iters=10, seed=bi); sync(); dt = time.perf_counter() - t0; Sv = al.score.cpu().numpy()
            for q, p in enumerate(Bt): rows.append(dict(a=p[0], b=p[1], n=len(X[p[0]]), m=len(X[p[1]]), batch=bi, batch_size=len(Bt), batch_sec=dt, amort_sec=dt / len(Bt), S=float(Sv[q])))
            del al, xx, yy
        sync(); tot = time.perf_counter() - T0; peak = torch.cuda.max_memory_allocated(); log('gpu', tag, len(prs), '쌍', len(B), '묶음', round(tot, 1), 's'); return rows, dict(total_sec=tot, n_batches=len(B), peak_mem_bytes=peak)
    meta = {}
    x0 = fs.pad_stack([X[ds[0]]], len(X[ds[0]])).to(dev, torch.float64); fs.align(x0, x0, n_random=127, iters=2); sync()        # 데우기
    pr = [(p['a'], p['b']) for p in P]
    r128, meta['d128_500'] = run_batches(pr, 127, 'D128'); r1, meta['d1_500'] = run_batches(pr, 0, 'D1'); wr('gpu_d128', r128); wr('gpu_d1', r1)
    rows = []; torch.cuda.reset_peak_memory_stats(); T0 = time.perf_counter()
    for k, p in enumerate(P):
        sync(); t0 = time.perf_counter(); c = fs.bnb(X[p['a']], X[p['b']], rel_tol=0.005, device=dev); sync(); dt = time.perf_counter() - t0
        rows.append(dict(a=p['a'], b=p['b'], n=p['n'], m=p['m'], wall_sec=dt, bnb_sec=float(c['sec']), certified=bool(c['certified'])))
        if k % 100 == 0: log('gpu bnb', k, round(time.time() - T0), 's')
    meta['bnb_500'] = dict(total_sec=time.perf_counter() - T0, peak_mem_bytes=torch.cuda.max_memory_allocated()); wr('gpu_bnb', rows)
    for fn, nr, tag in (('ava200.txt', 0, 'ava200_d1'), ('ava100.txt', 127, 'ava100_d128')):
        a = open(f'{T}/{fn}').read().split(); prs = [(a[i], a[j]) for i in range(len(a)) for j in range(i + 1, len(a))]
        rr, meta[tag] = run_batches(prs, nr, tag); meta[tag]['n_pairs'] = len(prs); wr('gpu_' + tag, rr)
    json.dump(dict(meta, device=torch.cuda.get_device_name(0), torch=torch.__version__), open(f'{OUTC}/gpu_meta.json', 'w'), indent=1); log('CT_GPU_DONE')


def part_ens():
    _src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]; sys.argv = ['r38_s0', ROOT, 'none']; g = {'__name__': 'r38s0'}; exec(_src, g)
    torch = g['torch']; info = json.load(open(f'{T}/prep.json')); rows = []
    for e in info['ens']:
        z = np.load(f"{T}/ens/{e['dom']}.npz"); EA = torch.as_tensor(z['EA'], dtype=torch.float64, device='cuda'); EB = torch.as_tensor(z['EB'], dtype=torch.float64, device='cuda')
        for rep in range(4):
            rng = np.random.default_rng(SEED + rep); torch.cuda.synchronize(); t0 = time.perf_counter()
            ms = g['MeanStart'](EA, EB); ia, ib = g['upairs'](0, 0, rng); d, _, S, ku, wid = g['red_pairs'](ms, ia, ib); torch.cuda.synchronize(); dt = time.perf_counter() - t0
            rows.append(dict(dom=e['dom'], ne=e['ne'], rep=rep, sec=dt, n_pairs=len(ia), mean_d=float(np.mean(d))))
            log('ens', e['dom'], e['ne'], rep, round(dt, 2), 's')
    wr('ens', rows); log('CT_ENS_DONE')


def q(v, p): return float(np.percentile(v, p))


def summ(v):
    v = np.asarray(v, float); return dict(n=int(len(v)), median=float(np.median(v)), p95=q(v, 95), max=float(v.max()), mean=float(v.mean()), total=float(v.sum()))


def slope(rows, k):
    nm = np.array([r['n'] * r['m'] for r in rows], float); t = np.array([r[k] for r in rows], float); ok = t > 0
    return float(np.polyfit(np.log(nm[ok]), np.log(t[ok]), 1)[0])


def part_sum():
    rd = lambda n: [json.loads(l) for l in open(f'{OUTC}/{n}.jsonl')] if os.path.exists(f'{OUTC}/{n}.jsonl') else []
    cpu, rep, d128c = rd('cpu'), rd('cpu_repeat'), rd('cpu_d128'); g128, g1, bnb = rd('gpu_d128'), rd('gpu_d1'), rd('gpu_bnb'); ens = rd('ens')
    S = dict(written=time.strftime('%F %T'), rule='round41 g5 (Drive JSON round41_2026-10-21 g5_design): 같은 장비, 측정 중 다른 계산 없음. 판정 없음.', per_pair={}, slope_log_nm={}, repeat_ratio={})
    for k, rows in (('tmalign_sec', cpu), ('usalign_sec', cpu), ('d1_cpu_sec', cpu), ('d128_cpu_sec', d128c), ('amort_sec|gpu_d128', g128), ('amort_sec|gpu_d1', g1), ('wall_sec|gpu_bnb', bnb)):
        kk = k.split('|')[0]
        if rows and kk in rows[0]: S['per_pair'][k] = summ([r[kk] for r in rows]); S['slope_log_nm'][k] = slope(rows, kk)
    byi = {i: r for i, r in enumerate(cpu)}
    for k in ('tmalign_sec', 'usalign_sec', 'd1_cpu_sec'):
        if rep: S['repeat_ratio'][k] = float(np.median([r[k] / byi[r['i']][k] for r in rep if byi[r['i']][k] > 0]))
    if bnb:
        S['bnb'] = dict(certified=int(sum(r['certified'] for r in bnb)), n=len(bnb), le60=int(sum(r['certified'] and r['wall_sec'] <= 60 for r in bnb)))
        mx = max(bnb, key=lambda r: r['wall_sec']); S['bnb'].update(max_pair=dict(a=mx['a'], b=mx['b'], n=mx['n'], m=mx['m'], sec=mx['wall_sec']))
        S['bnb_vs_static_sec'] = float(np.median([r['wall_sec'] / p['static_bnb_sec'] for r, p in zip(bnb, json.load(open(f'{T}/pairs.json'))) if p.get('static_bnb_sec')]))
    if cpu: S['agree_d1_cpu_vs_gpu_S_max_rel'] = float(max(abs(a['d1_S'] - b['S']) / max(abs(b['S']), 1e-12) for a, b in zip(sorted(cpu, key=lambda r: (r['a'], r['b'])), sorted(g1, key=lambda r: (r['a'], r['b']))))) if g1 else None
    for f in ('cpu_meta', 'foldseek', 'usava', 'gpu_meta'):
        if os.path.exists(f'{OUTC}/{f}.json'): S[f] = json.load(open(f'{OUTC}/{f}.json'))
    if 'foldseek' in S:
        for k, v in S['foldseek']['search'].items(): v['amort_sec_per_pair'] = v['sec'] / v['n_pairs_all']
        S['foldseek']['easy_search_pair_summ'] = summ([r['sec'] for r in S['foldseek']['easy_search_pair']])
    if 'usava' in S: S['usava']['amort_sec_per_pair'] = S['usava']['sec'] / S['usava']['n_pairs_expected']
    if 'gpu_meta' in S:
        for t in ('ava200_d1', 'ava100_d128'):
            m = S['gpu_meta'].get(t)
            if m: m['amort_sec_per_pair'] = m['total_sec'] / m['n_pairs']; m['extrap_1000_structures_sec'] = m['amort_sec_per_pair'] * 1000 * 999 / 2
    if ens: S['ensemble_cell'] = {e: summ([r['sec'] for r in ens if r['dom'] == e and r['rep'] > 0]) | dict(ne=next(r['ne'] for r in ens if r['dom'] == e)) for e in sorted({r['dom'] for r in ens})}
    S['prep'] = {k: v for k, v in json.load(open(f'{T}/prep.json')).items() if k not in ('template',)}
    json.dump(S, open(ROOT + '/results/r41_compute.json', 'w'), ensure_ascii=False, indent=1, default=float)
    dst = ROOT + '/derived/r41/compute'; os.makedirs(dst, exist_ok=True)
    for f in glob.glob(OUTC + '/*'): shutil.copy(f, dst)
    b = ROOT + '/results/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(ROOT + '/results/r41_compute.json', b + '/' + time.strftime('%H%M%S_') + 'r41_compute.json')
    log('CT_SUM_DONE')


if __name__ == '__main__':
    log('부분 시작', PART); dict(cpu=part_cpu, foldseek=part_foldseek, usava=part_usava, gpu=part_gpu, ens=part_ens, sum=part_sum)[PART]()
