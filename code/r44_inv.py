"""(round44, 등록 = Drive results/prereg_2026-10-28_addendum.json "r44_inv", 2026-10-28 UTC) 에너지 통계의 반복 수·프레임 수 불변성 측정과 E0′ 성분 분해.
부분:
  rep   : 줄기 D 실행 쌍 표(derived/r39/D/<dom>.json, r39_D_stem.py)에서 쪽마다 실행 5 → 3 (도메인·칸마다 고정 seed 로 부분집합 하나; 서술로 2 도).
          먼저 5 대 5 를 표에서 다시 내어 저장 칸 값과 같은지 본다. 칸 = 뒤 절반 창 MD T 대 MD 320 (T = 348·379·413·450). CPU.
          (모든 부분집합을 평균하면 항등으로 원래 값과 같아 시험이 되지 않으므로 부분집합은 하나만 뽑는다.)
  frm   : 실행당 프레임 100 → 50 (r38_s0 과 같은 실행·프레임에서 하나 건너 솎기). 줄기에서 고정 seed 로 고른 60 도메인(뒤 절반 창에서 320·379·413 K
          실행이 모두 5), T = 379·413. D 는 줄인 계산(실행 쌍마다 프레임 쌍 2,000), 맞춘 좌표는 전수. 점검: 100 판 맞춘 좌표 전수(aall) 한 실행 쌍을
          같은 코드로 다시 내어 저장값과 같은지. GPU.
  split : E0′ 성분 비 √(Δμ²₂_D·변/w̄² + ((ρ − 1)/(ρ + 1))²) 의 평균 이동 몫 m 과 퍼짐 몫 s (몫 = m/(m + s)). AlphaFlow·aSAMt·R3′·온도 칸. 서술. CPU.
  sum   : 판정(등록 그대로). 칸마다 줄인 판 / 원래 판의 상과 중앙 비 [상과 부트스트랩 2,000 의 95%]. 주 = Δ²_D/√n̄ 의 +59·+93 K; 보조 = 비·Δμ²₂_D·맞춘 좌표 판·+28·+130 K.
          선다 = 비 ∈ [0.9, 1.1] 이고 95% 구간이 1 을 포함. 반복 축·프레임 축마다 +59·+93 둘 다 서야 그 축이 선다.
사용: python r44_inv.py <ROOT> <부분들 쉼표>  출력: derived/r44/{rep.json, frm/<dom>.json, split.json}, results/r44_inv.json (사본 results/_backup/<날짜>/)"""
import os, sys, json, glob, time, itertools, zlib, shutil, traceback
import numpy as np
ROOT = sys.argv[1]; PARTS44 = sys.argv[2].split(',')
RES = ROOT + '/results'; DER = ROOT + '/derived'; OUT44 = DER + '/r44'; os.makedirs(OUT44 + '/frm', exist_ok=True)
TS = ('348', '379', '413', '450'); DTK = {'348': 28, '379': 59, '413': 93, '450': 130}; PRIMARY = ('379', '413')
NB = 2000; SEED = dict(rep=20261090, frm=20261091, boot=20261092, split=20261093); NFRM = 60; LO, HI = 0.9, 1.1


def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)
def seeded(*tag): return np.random.default_rng(zlib.crc32('|'.join(map(str, tag)).encode()))


def jdump(o, p, bk=True):
    json.dump(o, open(p + '.tmp', 'w'), ensure_ascii=False, indent=1, default=float); os.replace(p + '.tmp', p)
    if bk:
        b = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + os.path.basename(p))


# ---------- 상과 통계 (r39_core 와 같은 정의) ----------
def cath_sf():
    h = {}
    for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
        if l.startswith('#') or not l.strip(): continue
        f = l.split(); h[f[0]] = '.'.join(f[1:5])
    return h


def sf_groups(sfs):
    u = sorted(set(sfs)); ix = {s: i for i, s in enumerate(u)}; return u, np.array([ix[s] for s in sfs], int)


def sf_meds(v, g, ns):
    v = np.asarray(v, float); out = np.full(ns, np.nan)
    for s in range(ns):
        x = v[(g == s) & np.isfinite(v)]
        if len(x): out[s] = np.median(x)
    return out


def boot_joint(vs, g, ns, rng, B=NB):
    M = np.stack([sf_meds(v, g, ns) for v in vs]); pt = np.nanmedian(M, axis=1); idx = rng.integers(0, ns, (B, ns))
    return pt, np.nanmedian(M[:, idx], axis=2)


def ci(bs):
    bs = np.asarray(bs, float); bs = bs[np.isfinite(bs)]
    return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if len(bs) else [None, None]


def stem_D():
    return {j['dom']: j for j in (json.load(open(f)) for f in sorted(glob.glob(DER + '/r39/D/*.json')))}


def estats(C, BP, BQ, ne):
    D2 = 2 * C - BP - BQ; w = (BP + BQ) / 2
    return dict(C=float(C), B_P=float(BP), B_Q=float(BQ), d2n=float(D2 / np.sqrt(ne)), ratio=float(np.sqrt(max((C / w) ** 2 - 1, 0))) if w > 0 else float('nan'))


def cell_from(cr, fp, fq, ne):
    """cr·fp·fq: 실행 쌍마다 [d, d2, a2000, aall, a2_2000, a2all, ...] (r39_D_stem 의 KEYS 순서)."""
    m = lambda rows, k: float(np.mean([r[k] for r in rows]))
    return dict(D=estats(m(cr, 0), m(fp, 0), m(fq, 0), ne), A=estats(m(cr, 3), m(fp, 3), m(fq, 3), ne),
                dmu2_D=(m(cr, 1) - (m(fp, 1) + m(fq, 1)) / 2) / ne, dmu2_A=(m(cr, 5) - (m(fp, 5) + m(fq, 5)) / 2) / ne)


STATS = (('D.d2n', lambda c: c['D']['d2n']), ('D.ratio', lambda c: c['D']['ratio']), ('dmu2_D', lambda c: c['dmu2_D']),
         ('A.d2n', lambda c: c['A']['d2n']), ('A.ratio', lambda c: c['A']['ratio']), ('dmu2_A', lambda c: c['dmu2_A']))


def compare(rows, ref_key, sub_key, rng):
    """rows: [dict(dom, sf, <ref_key>=cell, <sub_key>=cell)]. 통계마다 상과 중앙(원래·줄인)과 비 [95%]."""
    u, g = sf_groups([r['sf'] for r in rows]); out = dict(n_dom=len(rows), n_sf=len(u))
    for nm, f in STATS:
        a = np.array([f(r[ref_key]) for r in rows]); b = np.array([f(r[sub_key]) for r in rows])
        pt, bs = boot_joint([a, b], g, len(u), rng); rr = bs[1] / bs[0]
        est = float(pt[1] / pt[0]) if pt[0] else None; c95 = ci(rr)
        out[nm] = dict(ref=float(pt[0]), sub=float(pt[1]), ratio=est, ci95=c95,
                       passed=bool(est is not None and LO <= est <= HI and c95[0] is not None and c95[0] <= 1 <= c95[1]))
    return out


# ---------- rep ----------
def part_rep():
    sfmap = cath_sf(); J = stem_D(); rows = {T: [] for T in TS}; chk = []
    for d, j in J.items():
        rp = j['rp'].get('half', {}); ne = j['ne']
        for T in TS:
            c = j['cells'].get(f'MD{T}-MD320|half')
            if not c or c.get('nP') != 5 or c.get('nQ') != 5: continue
            tab = lambda key: {(int(r[0]), int(r[1])): r[2:] for r in rp[key]}
            X, P, Q = tab(f'320-{T}'), tab('320-320'), tab(f'{T}-{T}')
            def cell(SP, SQ):
                SP, SQ = sorted(SP), sorted(SQ)
                return cell_from([X[(a, b)] for a in SP for b in SQ], [P[p] for p in itertools.combinations(SP, 2)], [Q[p] for p in itertools.combinations(SQ, 2)], ne)
            full = cell(range(5), range(5))
            chk.append(max(abs(full['D']['d2n'] - c['D']['d2n']), abs(full['dmu2_D'] - c['dmu2_D']), abs(full['A']['d2n'] - c['Aall']['d2n'])))
            rng = seeded('r44_rep', d, T)
            k3 = cell(rng.choice(5, 3, replace=False), rng.choice(5, 3, replace=False)); k2 = cell(rng.choice(5, 2, replace=False), rng.choice(5, 2, replace=False))
            rows[T].append(dict(dom=d, sf=sfmap.get(d, d), full=full, k3=k3, k2=k2))
    o = dict(written=time.strftime('%F %T'), check_max_abs_diff=float(max(chk)) if chk else None, n_check=len(chk), rows=rows)
    jdump(o, OUT44 + '/rep.json', bk=False); log('R44_REP', {T: len(v) for T, v in rows.items()}, '점검 최대 차', o['check_max_abs_diff']); log('R44_REP_DONE')


# ---------- frm (GPU) ----------
def part_frm():
    _e = open(ROOT + '/code/r38_e7_30.py').read(); _e = _e[:_e.index('\nfor d in SEL:')]
    G = dict(__name__='r38_e7_30_lib'); argv0 = sys.argv; sys.argv = ['r38_e7_30', ROOT]; exec(_e, G); sys.argv = argv0
    J = stem_D(); ok = [d for d, j in sorted(J.items()) if all(j['nruns'].get(f'{T}|half') == 5 for T in ('320', '379', '413')) and all(f'MD{T}-MD320|half' in j['cells'] for T in PRIMARY)]
    sel = sorted(np.random.default_rng(SEED['frm']).choice(ok, min(NFRM, len(ok)), replace=False).tolist())
    jdump(dict(written=time.strftime('%F %T'), n_eligible=len(ok), domains=sel), OUT44 + '/frm_domains.json', bk=False); log('frm 도메인', len(sel), '/', len(ok))
    np_, fs, cost, NP = G['np'], G['fs'], G['cost'], G['NPAIR']; nn = 50; TT = G['TT']
    for d in sel:
        fo = f'{OUT44}/frm/{d}.json'
        if os.path.exists(fo): continue
        t0 = time.time(); rng = seeded('r44_frm', d); G['_runs'].clear(); G['_md'].clear()
        ref = G['ref_of'](d); R = {T: G['dom_runs'](d, T, 'half', TT)[0] for T in ('320', '379', '413')}; L = G['dom_runs'](d, '320', 'half', TT)[1]; ne = L - 1
        res = dict(dom=d, ne=ne, nn=nn, cells={}, check={})
        a100 = {T: fs.align_to_reference(G['frames'](np.concatenate([X for _, X in R[T]])), ref).numpy().reshape(len(R[T]) * 100, -1) for T in ('320', '379')}
        stored = J[d]['rp']['half']['320-379'][0]
        res['check'] = dict(aall_00_recomputed=float(np.sqrt(cost(a100['320'][:100], a100['379'][:100])).mean()), aall_00_stored=float(stored[2 + 3]), run_pair=stored[:2])
        E = {T: G['coords'](np.concatenate([X[::2] for _, X in R[T]])) for T in R}
        A = {T: fs.align_to_reference(G['frames'](np.concatenate([X[::2] for _, X in R[T]])), ref).numpy().reshape(len(R[T]) * nn, -1) for T in R}
        def go(Ta, Tb, rr):
            ms = G['MeanStart'](E[Ta], E[Tb]); IA, IB, GR = [], [], []
            for k, (a, b) in enumerate(rr):
                f = rng.choice(nn * nn, NP, replace=False); IA.append(a * nn + f // nn); IB.append(b * nn + f % nn); GR.append(np.full(NP, k))
            IA, IB, GR = map(np.concatenate, (IA, IB, GR)); dv = G['red_pairs'](ms, IA, IB)[0]
            out = {}
            for k, (a, b) in enumerate(rr):
                cc = cost(A[Ta][a * nn:(a + 1) * nn], A[Tb][b * nn:(b + 1) * nn]); x = dv[GR == k]
                out[(a, b)] = [float(x.mean()), float((x ** 2).mean()), None, float(np.sqrt(cc).mean()), None, float(cc.mean())]
            return out
        P = go('320', '320', list(itertools.combinations(range(5), 2)))
        for T in PRIMARY:
            Q = go(T, T, list(itertools.combinations(range(5), 2))); X = go('320', T, list(itertools.product(range(5), range(5))))
            res['cells'][T] = cell_from(list(X.values()), list(P.values()), list(Q.values()), ne)
        res['sec'] = time.time() - t0; jdump(res, fo, bk=False); log('frm', d, round(res['sec']), 's', '점검', round(res['check']['aall_00_recomputed'] - res['check']['aall_00_stored'], 12))
    log('R44_FRM_DONE')


# ---------- split ----------
def part_split():
    sfmap = cath_sf(); rng = np.random.default_rng(SEED['split']); out = dict(written=time.strftime('%F %T'), cells={})
    def dom_split(C, BP, BQ, dmu2, ne):
        w = (BP + BQ) / 2; m = max(dmu2 * ne, 0.0) / w ** 2; rho = BQ / BP; s = ((rho - 1) / (rho + 1)) ** 2
        return dict(m=m, s=s, share_mean=m / (m + s) if m + s > 0 else np.nan, r_comp=float(np.sqrt(m + s)), r_dir=float(np.sqrt(max((C / w) ** 2 - 1, 0))), log_rho=float(np.log(rho)))
    def cell(rows):
        M = [(d, dom_split(**x)) for d, x in rows]
        if not M: return None
        u, g = sf_groups([sfmap.get(d, d) for d, _ in M]); ks = ('m', 's', 'share_mean', 'r_comp', 'r_dir', 'log_rho')
        pt, bs = boot_joint([np.array([v[k] for _, v in M]) for k in ks], g, len(u), rng)
        return dict(n_dom=len(M), n_sf=len(u), **{k: dict(est=float(pt[i]), ci95=ci(bs[i])) for i, k in enumerate(ks)})
    J = stem_D()
    for W in ('half', 'full'):
        for T in TS:
            rows = [(d, dict(C=c['D']['C'], BP=c['D']['B_P'], BQ=c['D']['B_Q'], dmu2=c['dmu2_D'], ne=j['ne'])) for d, j in J.items() for c in [j['cells'].get(f'MD{T}-MD320|{W}')] if c]
            out['cells'][f'온도 +{DTK[T]} K|줄기|{W}'] = cell(rows)
    jl = lambda sub: {os.path.basename(f)[:-5]: json.load(open(f)) for f in sorted(glob.glob(f'{DER}/r38/{sub}/*.json'))}
    F = jl('dff_D'); out['cells']['R3′ 역장 + 규약|full'] = cell([(j['dom'], dict(C=j['ff']['C'], BP=j['ff']['B_P'], BQ=j['ff']['B_Q'], dmu2=j['ff']['dmu'], ne=j['ne'])) for j in F.values() if 'ff' in j])
    for nm, sub, key in (('7c aSAMt320 대 MD 320', 'compare_D', 'aSAMt320-MD320|full|0.75'), ('7c AlphaFlow 대 MD 320', 'compare_D2', 'AlphaFlow-MD320|full|0.75')):
        S = jl(sub); out['cells'][nm + '|full'] = cell([(d, dict(C=c['C'], BP=c['B_P'], BQ=c['B_Q'], dmu2=c['dmu'], ne=j['ne'])) for d, j in S.items() for c in [j['cells'].get(key)] if c])
    out['rule'] = '서술(등록 r44_inv split). 도메인마다 m = max(Δμ²₂_D·변, 0)/w̄², s = ((ρ − 1)/(ρ + 1))², ρ = B_Q/B_P (Q = 비교 쪽), 몫 = m/(m + s); 상과 중앙 [상과 부트스트랩 2,000]. r39_E0 의 성분 비와 같은 식.'
    jdump(out, OUT44 + '/split.json', bk=False); log('R44_SPLIT', {k: (v['n_dom'], round(v['share_mean']['est'], 3), round(v['r_dir']['est'], 3)) for k, v in out['cells'].items() if v}); log('R44_SPLIT_DONE')


# ---------- sum ----------
def part_sum():
    sfmap = cath_sf(); J = stem_D(); out = dict(written=time.strftime('%F %T'), rule=('등록 r44_inv: 줄인 판 / 원래 판의 상과 중앙 비 [상과 부트스트랩 2,000 의 95%]. '
           '선다 = 비 ∈ [0.9, 1.1] 이고 구간이 1 을 포함. 주 = Δ²_D/√n̄ 의 +59·+93 K (뒤 절반). 축마다 둘 다 서야 그 축이 선다.'))
    if os.path.exists(OUT44 + '/rep.json'):
        rp = json.load(open(OUT44 + '/rep.json')); out['rep_check_max_abs_diff'] = rp['check_max_abs_diff']; out['rep'] = {}
        for T in TS:
            rows = rp['rows'][T]
            if len(rows) < 10: continue
            out['rep'][f'+{DTK[T]} K'] = dict(k3=compare(rows, 'full', 'k3', np.random.default_rng(SEED['boot'] + int(T))), k2=compare(rows, 'full', 'k2', np.random.default_rng(SEED['boot'] + 1000 + int(T))))
        pr = [out['rep'].get(f'+{DTK[T]} K', {}).get('k3', {}).get('D.d2n', {}).get('passed') for T in PRIMARY]
        out['rep_axis'] = dict(passed=bool(all(pr)), primary=pr)
    fr = [json.load(open(f)) for f in sorted(glob.glob(OUT44 + '/frm/*.json'))]
    if fr:
        out['frm_check_max_abs_diff'] = float(max(abs(r['check']['aall_00_recomputed'] - r['check']['aall_00_stored']) for r in fr)); out['frm'] = {}
        for T in PRIMARY:
            rows = []
            for r in fr:
                c = J[r['dom']]['cells'][f'MD{T}-MD320|half']
                ref = dict(D=dict(d2n=c['D']['d2n'], ratio=c['D']['ratio']), A=dict(d2n=c['Aall']['d2n'], ratio=c['Aall']['ratio']), dmu2_D=c['dmu2_D'], dmu2_A=c['dmu2_A'])
                rows.append(dict(dom=r['dom'], sf=sfmap.get(r['dom'], r['dom']), n100=ref, n50=r['cells'][T]))
            out['frm'][f'+{DTK[T]} K'] = compare(rows, 'n100', 'n50', np.random.default_rng(SEED['boot'] + 2000 + int(T)))
        pf = [out['frm'].get(f'+{DTK[T]} K', {}).get('D.d2n', {}).get('passed') for T in PRIMARY]
        out['frm_axis'] = dict(passed=bool(all(pf)), primary=pf, n_dom=len(fr))
    ra, fa = out.get('rep_axis', {}).get('passed'), out.get('frm_axis', {}).get('passed')
    out['sentence'] = ('반복 수(3 대 5)와 실행당 프레임 수(50 대 100)를 바꿔도 통계의 상과 중앙이 등록 여유(±10%) 안에서 같다' if ra and fa else
                       f'불변성 측정: 반복 축 {"섬" if ra else "안 섬"}, 프레임 축 {"섬" if fa else ("안 섬" if fa is False else "미완")}')
    jdump(out, RES + '/r44_inv.json'); log('R44_SUM', out['sentence']); log('R44_SUM_DONE')


if __name__ == '__main__':
    for part in PARTS44:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
    log('R44_INV_DONE', ','.join(PARTS44))
