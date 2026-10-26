"""(round39 CPU 분석; D 눈금 결정 2026-10-20, 범위 줄임 2026-10-20, 등록 = Drive JSON round39_2026-10-20) 시연 설계 2 절 E0′·E1·E5 (꼭) 와 E2 주 창·E3·E4 가벼운 판 (있으면 좋음).
자료: derived/r39/D/<dom>.json (줄기 336, D; r39_D_stem.py), derived/r38/dff_D (R3′ D), derived/r38/compare_D·compare_D2 (7c D), derived/r38/compare/<dom>.npz (7c 맞춘 좌표 실행 평균·변별 퍼짐),
      derived/r39/e5/<tag>.npz (R3′·R3″ 맞춘 좌표; r39_e5prep.py), derived/r39/pr.json (320 K 유효 차원; 이 파일의 pr 부분).
부분:
  check : 합성 자료 점검(E0′ 두 비, E1 IRLS, E3 귀무 검출률, E4 무잡음 역환산, E5 코사인 추정량). 실제 자료를 읽지 않는다.
  pr : 줄기 336 의 320 K 뒤 절반 앙상블(r39_D_stem 과 같은 실행·프레임) 참여 비 (Σλ)²/Σλ² ÷ 9 변 (E1 귀속 열). CPU.
  e0 : 두 영역 경계. 칸마다 도메인별 직접 비 √((C/w̄)² − 1), 성분 비 √(Δμ²₂_D·변/w̄² + ((ρ − 1)/(ρ + 1))²) (ρ = B_Q/B_P), 두 비의 차, δ, b − a.
  e1 : 크기 지수 α (Δ²_D = c·n^α, n = 변 수; IRLS 주 판, 민감도 넷, 생존 판, 귀속).
  e2 : (가벼운 판) 주 창 둘(전체·뒤 절반)의 D 값과 전체/뒤 절반 비.
  e3 : (가벼운 판) 검출력 곡선(D, Δ²_D/√변; 창 둘, 뽑는 법 둘).
  e4 : (가벼운 판) 빼고 맞히기, 공유 · 표적 3 판, K_μ (Δμ²₂_D), 곡선 셋, 창 둘.
  e5 : 방향 검정(맞춘 좌표): 평균 이동 코사인, 퍼짐 윤곽 순위 상관, 퍼짐 전체 부호; 천장·귀무; 판정 넷.
사용: python r39_core.py <ROOT> <부분들 쉼표> 출력: results/r39_<부분>.json (사본 results/_backup/<날짜>/)"""
import os
os.environ['CUDA_VISIBLE_DEVICES'] = ''
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '2'
import sys, json, glob, time, itertools, zlib, shutil, traceback, warnings
import numpy as np
from scipy.stats import spearmanr, theilslopes
from scipy.optimize import minimize_scalar
warnings.filterwarnings('ignore', category=RuntimeWarning)
ROOT = sys.argv[1]; PARTS39 = sys.argv[2].split(',')
RES = ROOT + '/results'; DER = ROOT + '/derived'; OUT39 = DER + '/r39'; os.makedirs(OUT39, exist_ok=True)
TS = ('348', '379', '413', '450'); DTK = {'348': 28, '379': 59, '413': 93, '450': 130}
SEED = dict(e0=20261070, e1=20261071, e2=20261072, e3=20261073, e4=20261074, e5=20261075, check=20261076)
NB = 2000


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def jdump(o, p, bk=True):
    json.dump(o, open(p + '.tmp', 'w'), ensure_ascii=False, indent=1, default=float); os.replace(p + '.tmp', p)
    if bk:
        b = RES + '/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + os.path.basename(p))


def cath_sf():
    h = {}
    for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
        if l.startswith('#') or not l.strip(): continue
        f = l.split(); h[f[0]] = '.'.join(f[1:5])
    return h


# ---------- 상과 통계 ----------
def sf_groups(sfs):
    u = sorted(set(sfs)); ix = {s: i for i, s in enumerate(u)}; return u, np.array([ix[s] for s in sfs], int)


def sf_meds(v, g, ns):
    v = np.asarray(v, float); out = np.full(ns, np.nan)
    for s in range(ns):
        x = v[(g == s) & np.isfinite(v)]
        if len(x): out[s] = np.median(x)
    return out


def boot_joint(vs, g, ns, rng, B=NB):
    """vs: 같은 도메인 순서의 배열들. 상과 중앙값(상과 안 중앙값의 중앙값)과 상과 복원 추출 B 번(모든 배열에 같은 추출)."""
    M = np.stack([sf_meds(v, g, ns) for v in vs]); pt = np.nanmedian(M, axis=1); idx = rng.integers(0, ns, (B, ns))
    return pt, np.nanmedian(M[:, idx], axis=2)


def ci(bs):
    bs = np.asarray(bs, float); bs = bs[np.isfinite(bs)]
    return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if len(bs) else [None, None]


def summ(pt, bs):
    return dict(est=float(pt) if np.isfinite(pt) else None, ci95=ci(bs))


# ---------- 자료 ----------
def stem_D():
    J = {}
    for f in sorted(glob.glob(OUT39 + '/D/*.json')):
        j = json.load(open(f)); J[j['dom']] = j
    return J


def jload_dir(sub):
    return {os.path.basename(f)[:-5]: json.load(open(f)) for f in sorted(glob.glob(f'{DER}/r38/{sub}/*.json'))}


# ---------- E0′ ----------
def e0_dom(C, BP, BQ, dmu2, ne):
    w = (BP + BQ) / 2; rd = float(np.sqrt(max((C / w) ** 2 - 1, 0))); dm = max(dmu2 * ne, 0.0); rho = BQ / BP
    rc = float(np.sqrt(dm / w ** 2 + ((rho - 1) / (rho + 1)) ** 2)); sbar = w / (np.sqrt(18) * np.sqrt(ne))
    return dict(r_dir=rd, r_comp=rc, diff=rd - rc, rel=(rc / rd - 1) if rd > 0 else np.nan, delta=float(np.sqrt(max(dmu2, 0))), b_minus_a=float(2 * sbar * (rho - 1) / (rho + 1)),
                wn=float(w / np.sqrt(ne)), d2n=float((2 * C - BP - BQ) / np.sqrt(ne)), bound_d2n=float(2 * (np.sqrt(2) - 1) * w / np.sqrt(ne)))


def region(r):
    return None if r is None else ('가까움' if r < 0.5 else ('멂' if r > 2 else '교차'))


def e0_cell(rows, sfmap, rng):
    """rows: [(dom, dict(C, B_P, B_Q, dmu2, ne))]."""
    M = [(d, e0_dom(x['C'], x['B_P'], x['B_Q'], x['dmu2'], x['ne'])) for d, x in rows]
    if not M: return None
    u, g = sf_groups([sfmap.get(d, d) for d, _ in M]); ks = ('r_dir', 'r_comp', 'diff', 'rel', 'delta', 'b_minus_a', 'wn', 'd2n', 'bound_d2n')
    pt, bs = boot_joint([np.array([m[k] for _, m in M]) for k in ks], g, len(u), rng)
    o = dict(n_dom=len(M), n_sf=len(u), **{k: summ(pt[i], bs[i]) for i, k in enumerate(ks)})
    rd = np.array([m['r_dir'] for _, m in M]); o['region'] = region(o['r_dir']['est']); o['side_of_1'] = '비 < 1' if o['r_dir']['est'] < 1 else '비 ≥ 1'
    o['frac_dom'] = {'가까움': float(np.mean(rd < 0.5)), '교차': float(np.mean((rd >= 0.5) & (rd <= 2))), '멂': float(np.mean(rd > 2))}
    o['self_check'] = '맞음' if o['rel']['est'] is not None and abs(o['rel']['est']) <= 0.2 else '어긋남'
    return o


def part_e0():
    rng = np.random.default_rng(SEED['e0']); sfmap = cath_sf(); J = stem_D(); out = dict(written=time.strftime('%F %T'), cells={})
    for W in ('half', 'full'):
        for T in TS:
            rows = [(d, dict(C=c['D']['C'], B_P=c['D']['B_P'], B_Q=c['D']['B_Q'], dmu2=c['dmu2_D'], ne=j['ne'])) for d, j in J.items() for c in [j['cells'].get(f'MD{T}-MD320|{W}')] if c]
            out['cells'][f'온도 +{DTK[T]} K|줄기|{W}'] = e0_cell(rows, sfmap, rng)
    F = jload_dir('dff_D'); rows = [(j['dom'], dict(C=j['ff']['C'], B_P=j['ff']['B_P'], B_Q=j['ff']['B_Q'], dmu2=j['ff']['dmu'], ne=j['ne'])) for j in F.values() if 'ff' in j]
    out['cells']['R3′ 역장 + 규약 + ΔT 20 K|full'] = e0_cell(rows, sfmap, rng)
    CD = jload_dir('compare_D'); CD2 = jload_dir('compare_D2')
    for nm, src, key in (('7c aSAMt320 대 MD 320', CD, 'aSAMt320-MD320|full|0.75'), ('7c AlphaFlow 대 MD 320', CD2, 'AlphaFlow-MD320|full|0.75')):
        rows = [(d, dict(C=c['C'], B_P=c['B_P'], B_Q=c['B_Q'], dmu2=c['dmu'], ne=j['ne'])) for d, j in src.items() for c in [j['cells'].get(key)] if c]
        out['cells'][nm + '|full'] = e0_cell(rows, sfmap, rng)
    for T in TS:
        rows = [(d, dict(C=c['C'], B_P=c['B_P'], B_Q=c['B_Q'], dmu2=c['dmu'], ne=j['ne'])) for d, j in CD.items() for c in [j['cells'].get(f'MD{T}-MD320|full|0.75')] if c]
        out['cells'][f'(서술) 시험 67 온도 +{DTK[T]} K|full'] = e0_cell(rows, sfmap, rng)
    out['not_in_D'] = ['R3″ 규약만', '치환']; out['not_in_D_note'] ='D 판 계산이 없다(범위 줄임). 길이가 다른 행은 prereg 1 절 예측 3 으로 정한다.'
    out['rule'] = ('round39 E0′ (D): 칸마다 도메인별 직접 비·성분 비, 상과 중앙값 [상과 부트스트랩 2,000 의 95%]. 영역은 직접 비의 상과 중앙값으로 가까움 < 0.5, 교차 0.5–2, 멂 > 2 (경계 비 = 1). '
                   '자가 점검 = 성분 비/직접 비 − 1 의 상과 중앙값이 ±0.2 안이면 맞음. Δμ²₂_D 가 음수면 성분 비에서 0 으로 자른다.')
    jdump(out, RES + '/r39_E0.json'); log('R39_E0', {k: (v['n_dom'], round(v['r_dir']['est'], 3), v['region'], v['self_check']) for k, v in out['cells'].items() if v}); log('R39_E0_DONE')


# ---------- E1 ----------
def _c(n, y, w, a):
    x = n ** a; return (w * y * x).sum() / (w * x * x).sum()


def _rss(a, n, y, w):
    return (w * (y - _c(n, y, w, a) * n ** a) ** 2).sum()


def _rss_fe(a, n, y, w, g, ns):
    x = n ** a; c = np.bincount(g, w * y * x, ns) / np.maximum(np.bincount(g, w * x * x, ns), 1e-300); return (w * (y - c[g] * x) ** 2).sum()


def nls(n, y, w, g=None, ns=None):
    f = (lambda a: _rss(a, n, y, w)) if g is None else (lambda a: _rss_fe(a, n, y, w, g, ns))
    return float(minimize_scalar(f, bounds=(-2, 4), method='bounded', options=dict(xatol=1e-6)).x)


def irls(n, y, g=None, ns=None):
    """원 척도 비선형 적합 y = c·n^α. 무가중으로 시작 → log|잔차| ~ log n 기울기 × 2 = γ → 무게 n^(−γ) 로 다시 적합; 세 번 또는 α 변화 < 0.01 까지."""
    w = np.ones_like(y); a = nls(n, y, w, g, ns); gam = 0.0; it = 0
    for it in range(1, 4):
        if g is None: pred = _c(n, y, w, a) * n ** a
        else:
            x = n ** a; c = np.bincount(g, w * y * x, ns) / np.maximum(np.bincount(g, w * x * x, ns), 1e-300); pred = c[g] * x
        r = np.abs(y - pred); ok = r > 0
        gam = 2 * float(np.polyfit(np.log(n[ok]), np.log(r[ok]), 1)[0]); w = n ** (-gam); w = w / w.mean(); a2 = nls(n, y, w, g, ns); ch = abs(a2 - a); a = a2
        if ch < 0.01: break
    return a, gam, it


def wmedian(v, w):
    o = np.argsort(v); cw = np.cumsum(w[o]); return v[o][np.searchsorted(cw, cw[-1] / 2)]


def lad(n, y):
    def obj(a):
        x = n ** a; c = wmedian(y / x, x); return np.abs(y - c * x).sum()
    grid = np.arange(-1, 3.0001, 0.01); a0 = grid[int(np.argmin([obj(a) for a in grid]))]
    return float(minimize_scalar(obj, bounds=(a0 - 0.02, a0 + 0.02), method='bounded', options=dict(xatol=1e-6)).x)


def e1_data(J, T, W):
    D = [(d, j['ne'], j['cells'][f'MD{T}-MD320|{W}']) for d, j in J.items() if f'MD{T}-MD320|{W}' in j['cells']]
    return [x[0] for x in D], np.array([x[1] for x in D], float), np.array([x[2]['D']['D2'] for x in D]), D


def part_e1():
    rng = np.random.default_rng(SEED['e1']); sfmap = cath_sf(); J = stem_D(); W = 'half'; out = dict(written=time.strftime('%F %T'), window=W, steps={})
    alld = sorted(J); u, g = sf_groups([sfmap.get(d, d) for d in alld]); ns = len(u); mem = [np.where(g == s)[0] for s in range(ns)]; pos = {d: i for i, d in enumerate(alld)}
    surv = [d for d in alld if all(f'MD{T}-MD320|{W}' in J[d]['cells'] for T in TS)]; out['n_survival'] = len(surv); out['survival_note'] = '뒤 절반 창에서 네 온도 칸이 모두 있는 도메인(D 판의 생존; round37 의 149 와 같은 정의)'
    pr = json.load(open(OUT39 + '/pr.json')) if os.path.exists(OUT39 + '/pr.json') else {}
    BS = {}; bidx = [np.concatenate([mem[s] for s in rng.integers(0, ns, ns)]) for _ in range(NB)]; ssurv = set(surv)
    for T in TS:
        doms, n, y, D = e1_data(J, T, W); o = dict(dT=DTK[T], n_dom=len(doms))
        q1, q4 = np.percentile(n, [25, 75]); lo, hi = n <= q1, n >= q4; rn = float(np.median(n[hi]) / np.median(n[lo])); d2n = y / np.sqrt(n)
        o['r_n'] = rn; o['quartile_alpha'] = float(0.5 + np.log(np.median(d2n[hi]) / np.median(d2n[lo])) / np.log(rn)) if np.median(d2n[lo]) > 0 and np.median(d2n[hi]) > 0 else None
        a, gam, it = irls(n, y); o['irls'] = dict(alpha=a, gamma=gam, iters=it)
        o['lad_alpha'] = lad(n, y); p = y > 0; o['log_alpha'] = float(np.polyfit(np.log(n[p]), np.log(y[p]), 1)[0]); o['log_excluded'] = int((~p).sum())
        o['theilsen_alpha'] = float(theilslopes(np.log(y[p]), np.log(n[p]))[0])
        gg = np.array([g[pos[d]] for d in doms]); o['sf_fixed_alpha'] = irls(n, y, gg, ns)[0]
        sv = np.array([d in ssurv for d in doms]); o['survival_alpha'] = irls(n[sv], y[sv])[0]
        dm = np.array([x[2]['dmu2_D'] for x in D]); pp = dm > 0; o['dmu_n_slope'] = float(np.polyfit(np.log(n[pp]), np.log(dm[pp]), 1)[0]); o['dmu_n_excluded'] = int((~pp).sum())
        ls = np.array([2 * np.log(x[2]['D']['B_Q'] / x[2]['D']['B_P']) for x in D]); o['logspread_n_slope'] = float(np.polyfit(np.log(n), ls, 1)[0])
        out['steps'][T] = o
        # 상과 부트스트랩: 주 판(IRLS)과 생존 판
        dpos = {d: i for i, d in enumerate(doms)}; amap = np.array([dpos.get(d, -1) for d in alld]); smap = np.array([dpos.get(d, -1) if d in ssurv else -1 for d in alld])
        ba, bsv = [], []
        for bi in bidx:
            k = amap[bi]; k = k[k >= 0]; ks = smap[bi]; ks = ks[ks >= 0]
            ba.append(irls(n[k], y[k])[0] if len(k) > 5 else np.nan); bsv.append(irls(n[ks], y[ks])[0] if len(ks) > 5 else np.nan)
        BS[T] = (np.array(ba), np.array(bsv)); o['irls']['ci95'] = ci(BS[T][0]); o['survival_ci95'] = ci(BS[T][1])
        log('E1', T, o['n_dom'], 'α IRLS', round(a, 3), o['irls']['ci95'], 'γ', round(gam, 2), '생존', round(o['survival_alpha'], 3))
    dd = BS['450'][0] - BS['379'][0]; ds = BS['450'][1] - BS['379'][1]; st = out['steps']
    out['diff_130_59'] = dict(est=st['450']['irls']['alpha'] - st['379']['irls']['alpha'], ci95=ci(dd), lower_2p5=ci(dd)[0])
    out['survival_diff_130_59'] = dict(est=st['450']['survival_alpha'] - st['379']['survival_alpha'], ci95=ci(ds))
    c1 = st['379']['irls']['alpha'] < st['413']['irls']['alpha'] < st['450']['irls']['alpha']; c2 = ci(dd)[0] is not None and ci(dd)[0] > 0
    c3 = {T: not (st[T]['irls']['ci95'][0] <= 0.5 <= st[T]['irls']['ci95'][1]) for T in ('379', '413', '450')}
    out['prediction'] = dict(order_59_93_130=bool(c1), diff_lower_gt0=bool(c2), ci_excludes_half=c3, verdict='맞음' if (c1 and c2 and all(c3.values())) else '틀림')
    out['survival_reading'] = '생존 판에서 ΔT 의존이 사라짐(차의 95% 구간이 0 을 포함) → 주된 출처는 생존 거르기' if (ci(ds)[0] is not None and ci(ds)[0] <= 0 <= ci(ds)[1]) else '생존 판에서도 ΔT 의존이 남음 → 귀속 열(δ·σ·유효 차원)로 출처를 쓴다'
    if pr:
        pk = [d for d in pr if 'pr' in pr[d]]; nn = np.array([pr[d]['ne'] for d in pk], float); f = np.array([pr[d]['pr_frac'] for d in pk]); out['pr_n_slope'] = float(np.polyfit(np.log(nn), np.log(f), 1)[0]); out['pr_n'] = len(pk); out['pr_skipped'] = len(pr) - len(pk)
    fl = [(j['ne'], np.std([j['cells'][f'null320|{W}|{i}']['D']['D2'] for i in range(10)])) for j in J.values() if f'null320|{W}|9' in j['cells']]
    if fl: a_ = np.array(fl, float); ok = a_[:, 1] > 0; out['floor_sd_n_slope'] = float(np.polyfit(np.log(a_[ok, 0]), np.log(a_[ok, 1]), 1)[0]); out['floor_n'] = int(ok.sum())
    # 자연 점검: 시험 67 (+130 K, 전체 창, D) 대 줄기 전체 창
    CD = jload_dir('compare_D'); _, nS, yS, _ = e1_data(J, '450', 'full'); aF = irls(nS, yS)[0]
    tv = [(j['ne'], j['cells']['MD450-MD320|full|0.75']['d2n']) for j in CD.values() if 'MD450-MD320|full|0.75' in j['cells']]
    if tv:
        t_ = np.array(tv, float); obs = np.median(t_[:, 1]) / np.median(yS / np.sqrt(nS)); prd = np.median(t_[:, 0] ** (aF - 0.5)) / np.median(nS ** (aF - 0.5))
        out['natural_check'] = dict(n_test=len(tv), alpha_full_450=aF, observed_ratio=float(obs), predicted_ratio_from_n=float(prd), note='도메인 중앙값 비(서술). 시험 87 → D 판이 있는 67 (P2 집합)')
    out['rule'] = ('round39 E1 (D, 뒤 절반 창 주): Δ²_D = c·n^α (n = 변 수), IRLS(Var ∝ n^γ, γ 를 함께 추정), 상과 부트스트랩 2,000 (모든 걸음에 같은 상과 추출). 예측: α(+59) < α(+93) < α(+130), '
                   'α(+130) − α(+59) 의 부트스트랩 2.5% 점 > 0, 세 걸음 모두 95% 구간이 0.5 를 뺌. D 에서 미리 적은 어림은 없음(4분위 α 는 서술). 생존 판에서 차의 95% 구간이 0 을 포함하면 "생존 거르기" 로 읽는다.')
    jdump(out, RES + '/r39_E1.json'); log('R39_E1', out['prediction'], out['survival_reading']); log('R39_E1_DONE')


# ---------- E2 (가벼운 판) ----------
def r0_band(r0, mem, rng, draws, scheme):
    """R0 띠: 도메인마다 10 나눔 가운데 하나를 무작위로, 뽑는 법(scheme)으로 상과를 고른 뒤 상과 중앙값. 2.5·97.5% 점."""
    v = []
    for _ in range(draws):
        S = scheme(); v.append(np.median([np.median(r0[mem[s], rng.integers(0, r0.shape[1], len(mem[s]))]) for s in S]))
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def part_e2():
    rng = np.random.default_rng(SEED['e2']); sfmap = cath_sf(); J = stem_D(); out = dict(written=time.strftime('%F %T'), steps={})
    for T in TS:
        doms = [d for d, j in J.items() if f'MD{T}-MD320|half' in j['cells'] and f'MD{T}-MD320|full' in j['cells']]
        u, g = sf_groups([sfmap.get(d, d) for d in doms]); ns = len(u)
        h = np.array([J[d]['cells'][f'MD{T}-MD320|half']['D']['d2n'] for d in doms]); f = np.array([J[d]['cells'][f'MD{T}-MD320|full']['D']['d2n'] for d in doms])
        mh = np.array([J[d]['cells'][f'MD{T}-MD320|half']['dmu2_D'] for d in doms]); mf = np.array([J[d]['cells'][f'MD{T}-MD320|full']['dmu2_D'] for d in doms])
        pt, bs = boot_joint([h, f, mh, mf], g, ns, rng)
        out['steps'][T] = dict(dT=DTK[T], n_dom=len(doms), n_sf=ns, d2n_half=summ(pt[0], bs[0]), d2n_full=summ(pt[1], bs[1]), ratio_full_over_half=summ(pt[1] / pt[0], bs[1] / bs[0]),
                               dmu2_half=summ(pt[2], bs[2]), dmu2_full=summ(pt[3], bs[3]), dmu2_ratio_full_over_half=summ(pt[3] / pt[2], bs[3] / bs[2]))
    for W in ('half', 'full'):
        doms = [d for d, j in J.items() if f'null320|{W}|9' in j['cells']]; u, g = sf_groups([sfmap.get(d, d) for d in doms]); ns = len(u); mem = [np.where(g == s)[0] for s in range(ns)]
        r0 = np.array([[J[d]['cells'][f'null320|{W}|{i}']['D']['d2n'] for i in range(10)] for d in doms])
        out[f'R0_band_{W}'] = dict(n_dom=len(doms), band=r0_band(r0, mem, rng, NB, lambda: rng.integers(0, ns, ns)))
    out['rule'] = 'round39 E2 가벼운 판 (D): 주 창 둘(전체·뒤 절반, N = 100)의 Δ²_D/√변 과 Δμ²₂_D/변, 같은 도메인 집합의 전체/뒤 절반 비(상과 중앙값의 비, 같은 상과 추출). 서술. 보충 창·Q 문턱은 맞춘 좌표 판(이 범위 밖).'
    jdump(out, RES + '/r39_E2.json'); log('R39_E2', {T: round(v['ratio_full_over_half']['est'], 3) for T, v in out['steps'].items()}); log('R39_E2_DONE')


# ---------- E3 (가벼운 판) ----------
NS3 = (5, 10, 20, 50, 100, 200, 'all')


def scheme_of(meth, n, ns, mem, rng):
    tot = sum(len(m) for m in mem)
    if meth == '가':
        k = ns if n == 'all' else n; return lambda: rng.integers(0, ns, k)
    tgt = tot if n == 'all' else n
    def f():
        S, c = [], 0
        while c < tgt:
            s = int(rng.integers(0, ns)); S.append(s); c += len(mem[s])
        return np.array(S)
    return f


def sign_perm(dif, rng):
    k = len(dif); obs = np.median(dif)
    if k <= 10: sg = np.array(list(itertools.product((1, -1), repeat=k)))
    else: sg = np.vstack([np.ones(k), rng.choice((1, -1), (999, k))])
    return float(np.mean(np.median(sg * dif[None], axis=1) >= obs))


def part_e3():
    rng = np.random.default_rng(SEED['e3']); sfmap = cath_sf(); J = stem_D(); out = dict(written=time.strftime('%F %T'), rows=[], n08={}); R3 = 500
    for W in ('half', 'full'):
        for T in TS:
            doms = [d for d, j in J.items() if f'MD{T}-MD320|{W}' in j['cells'] and f'null320|{W}|9' in j['cells']]
            u, g = sf_groups([sfmap.get(d, d) for d in doms]); ns = len(u); mem = [np.where(g == s)[0] for s in range(ns)]
            v = np.array([J[d]['cells'][f'MD{T}-MD320|{W}']['D']['d2n'] for d in doms]); r0 = np.array([[J[d]['cells'][f'null320|{W}|{i}']['D']['d2n'] for i in range(10)] for d in doms])
            mT = np.array([np.median(v[m]) for m in mem])
            for meth in ('가', '나'):
                curve = []
                for n in NS3:
                    sch = scheme_of(meth, n, ns, mem, rng); band = r0_band(r0, mem, rng, NB, sch); det, nsf, ndm, kinds = [], [], [], []
                    for _ in range(R3):
                        S = sch(); nsf.append(len(S)); ndm.append(int(sum(len(mem[s]) for s in S)))
                        if len(S) >= 15:
                            lb = np.percentile(np.median(mT[S][rng.integers(0, len(S), (200, len(S)))], axis=1), 2.5); det.append(lb > band[1]); kinds.append('부트스트랩')
                        else:
                            dif = np.array([np.median(v[mem[s]] - r0[mem[s], rng.integers(0, 10, len(mem[s]))]) for s in S]); det.append(sign_perm(dif, rng) < 0.05); kinds.append('부호 순열')
                    p = float(np.mean(det)); row = dict(window=W, T=T, dT=DTK[T], method=meth, n=n, p_detect=p, mc_se=float(np.sqrt(p * (1 - p) / R3)), n_sf_median=float(np.median(nsf)), n_dom_median=float(np.median(ndm)),
                                                        test='부호 순열(참고)' if kinds.count('부호 순열') > R3 / 2 else '부트스트랩', R0_band=list(band))
                    out['rows'].append(row); curve.append((n, p))
                n08 = next((n for n, p in curve if p >= 0.8), None); out['n08'][f'{W}|{T}|{meth}'] = n08
                log('E3', W, T, meth, [(n, round(p, 2)) for n, p in curve], 'n0.8', n08)
    h = out['n08']; p20 = {T: next(r['p_detect'] for r in out['rows'] if r['window'] == 'half' and r['T'] == T and r['method'] == '나' and r['n'] == 20) for T in ('379', '413', '450')}
    n28 = h['half|348|나']
    out['prediction'] = dict(p_n20_half_na=p20, pred1=all(p >= 0.8 for p in p20.values()), n08_28K=n28,
                             pred2=('맞음' if n28 in (100, 200, 'all') else ('도달 못 함' if n28 is None else '맞지 않음(더 적은 n 에서 검출)')))
    out['rule'] = ('round39 E3 가벼운 판 (D, Δ²_D/√변): 걸음 × 창(뒤 절반·전체) × 뽑는 법(가 = 상과 n 개 복원, 나 = 도메인 n 개가 될 때까지 상과 복원) × n, 500 회. '
                   '상과 수 ≥ 15 이면 상과 중앙값의 부트스트랩(200) 2.5% 점이 같은 뽑는 법·n 으로 다시 낸 R0 띠(도메인마다 나눔 하나 무작위, 2,000 회)의 97.5% 점을 넘으면 검출; '
                   '< 15 이면 상과 단위 부호 순열(도메인 칸 − 같은 도메인 R0 나눔 하나) 일방 p < 0.05. 예측(나, 뒤 절반): +59 K 이상 n = 20 에서 ≥ 0.8, +28 K 의 n0.8 ≥ 100.')
    jdump(out, RES + '/r39_E3.json'); log('R39_E3', out['prediction']); log('R39_E3_DONE')


# ---------- E4 (가벼운 판) ----------
def rp_d2(j, W, key):
    rows = j['rp'][W].get(key); k = j['rp_keys'].index('d2') + 2
    return None if rows is None else {(int(r[0]), int(r[1])): r[k] for r in rows}


def dmu2_sub(j, W, T, S):
    cr = rp_d2(j, W, f'320-{T}'); ff = rp_d2(j, W, '320-320'); tt = rp_d2(j, W, f'{T}-{T}')
    cross = np.mean([v for (k, t), v in cr.items() if t in S]); fT = np.mean([tt[(a, b)] for a, b in itertools.combinations(sorted(S), 2)])
    return (cross - (np.mean(list(ff.values())) + fT) / 2) / j['ne']


def pw(xk, yk):
    """단조(누적 최대)로 만든 조각 선형 곡선의 역함수. 마지막 구간을 선형으로 늘린다. 반환 inv(y) -> (x, 표시)."""
    xk = np.asarray(xk, float); yk = np.maximum.accumulate(np.asarray(yk, float))
    def inv(y):
        if y <= 0: return 0.0, '0 이하'
        for i in range(1, len(xk)):
            if y <= yk[i] and yk[i] > yk[i - 1]: return float(xk[i - 1] + (y - yk[i - 1]) * (xk[i] - xk[i - 1]) / (yk[i] - yk[i - 1])), ''
        s = (yk[-1] - yk[-2]) / (xk[-1] - xk[-2]) if len(xk) > 1 else 0
        return (float(xk[-1] + (y - yk[-1]) / s), '외삽') if s > 0 else (float(xk[-1]), '평평')
    return inv


def part_e4():
    rng = np.random.default_rng(SEED['e4']); sfmap = cath_sf(); J = stem_D(); out = dict(written=time.strftime('%F %T'), rows=[], checks={})
    for W in ('half', 'full'):
        Y = {d: {T: j['cells'][f'MD{T}-MD320|{W}']['dmu2_D'] for T in TS if f'MD{T}-MD320|{W}' in j['cells']} for d, j in J.items()}
        worst = 0.0
        for d, j in list(J.items())[:40]:
            for T in Y[d]:
                nT = len(j['rp'][W]['runs'][T]); v = dmu2_sub(j, W, T, set(range(nT))); worst = max(worst, abs(v - Y[d][T]) / max(abs(Y[d][T]), 1e-12))
        out['checks'][f'subset_all_equals_cell|{W}'] = worst; assert worst < 1e-8, f'실행 쌍 표 재구성 불일치 {worst}'
        for Tx in TS:
            K = DTK[Tx]; rest = [T for T in TS if T != Tx]
            if Tx == '450': rest = ['348', '379', '413']
            dl = list(J); uA, gA = sf_groups([sfmap.get(d, d) for d in dl])
            gk = {T: float(np.nanmedian(sf_meds(np.array([Y[d].get(T, np.nan) for d in dl]), gA, len(uA)))) for T in rest}
            xs = [0] + [DTK[T] for T in rest]; gs = [0] + [gk[T] for T in rest]; ginv = pw(xs, gs); gfun = lambda x: np.interp(x, xs, np.maximum.accumulate(gs))
            res = {c: [] for c in ('모양×진폭', '선형 보간', '모집단')}; doms = []
            for d, j in J.items():
                if Tx not in Y[d] or len(j['rp'][W]['runs'][Tx]) < 3: continue
                pts = [T for T in rest if T in Y[d]]
                if not pts: continue
                amp = sum(Y[d][T] for T in pts) / sum(gfun(DTK[T]) for T in pts)
                linv = pw([0] + [DTK[T] for T in pts], [0] + [Y[d][T] for T in pts])
                nT = len(j['rp'][W]['runs'][Tx]); ks = {c: [] for c in res}; fl = {c: [] for c in res}
                for S in itertools.combinations(range(nT), 3):
                    y = dmu2_sub(j, W, Tx, set(S))
                    for c, (k, f) in (('모양×진폭', ginv(y / amp) if amp > 0 else (np.nan, '진폭 ≤ 0')), ('선형 보간', linv(y)), ('모집단', ginv(y))):
                        ks[c].append(k); fl[c].append(f)
                doms.append(d)
                for c in res: res[c].append((d, np.array(ks[c], float), fl[c]))
            sfs = [sfmap.get(d, d) for d in doms]; u, g = sf_groups(sfs); ns = len(u)
            for c, rr in res.items():
                Kd = np.array([np.nanmean(k) for _, k, _ in rr]); ae = np.array([np.nanmean(np.abs(k - K)) for _, k, _ in rr]); sd = np.array([np.nanstd(k) for _, k, _ in rr])
                pt, bs = boot_joint([Kd], g, ns, rng); fls = [f for _, _, ff in rr for f in ff]
                out['rows'].append(dict(window=W, T=Tx, K=K, curve=c, n_dom=len(rr), n_sf=ns, n_splits_mean=float(np.mean([len(k) for _, k, _ in rr])) if rr else None, split_sd=float(np.nanmean(sd)) if rr else None,
                                        K_sf=summ(pt[0], bs[0]), bias=float(pt[0] - K), dom_MAE=float(np.nanmean(ae)), dom_err_90w=float(np.nanpercentile(Kd - K, 95) - np.nanpercentile(Kd - K, 5)) if rr else None,
                                        frac_zero=float(np.mean([f == '0 이하' for f in fls])) if fls else None, frac_extrap=float(np.mean([f == '외삽' for f in fls])) if fls else None, frac_flat=float(np.mean([f == '평평' for f in fls])) if fls else None,
                                        doms=[d for d, _, _ in rr], ae=ae.tolist()))
    def row(W, T, c): return next(r for r in out['rows'] if r['window'] == W and r['T'] == T and r['curve'] == c)
    V = {}
    for W in ('half', 'full'):
        mae = {c: np.mean([row(W, T, c)['dom_MAE'] for T in ('379', '413')]) for c in ('모양×진폭', '선형 보간', '모집단')}
        main = '선형 보간' if mae['모양×진폭'] > mae['선형 보간'] else '모양×진폭'
        gate = {T: abs(row(W, T, main)['bias']) <= 10 for T in ('379', '413')}
        r450, r413 = row(W, '450', main), row(W, '413', main); s450 = set(r450['doms']); m413 = np.mean([a for d, a in zip(r413['doms'], r413['ae']) if d in s450]) if s450 else np.nan
        p450 = bool(abs(r450['bias']) <= 10 and r450['dom_MAE'] <= 1.5 * m413)
        V[W] = dict(MAE_379_413=mae, main_curve=main, replaced=main != '모양×진폭', gate_le10K=gate, gate_pass=all(gate.values()), pass_450=p450, MAE_413_same_set=float(m413), last_trusted='+130' if p450 else '+93',
                    pred_order=bool(mae['모양×진폭'] <= mae['선형 보간'] < mae['모집단']))
    out['verdict'] = V
    for r in out['rows']: r.pop('ae', None)
    out['rule'] = ('round39 E4 가벼운 판 (D, Δμ²₂_D): 공유 · 표적 3 판(뺀 온도의 반복 셋 = 표적, 320 K 는 곡선과 같은 다섯; 실행 쌍 표에서 다시 냄) × 곡선 셋 × 창 둘(뒤 절반 주·전체). '
                   '모양 g = 뺀 온도 없이 만든 상과 중앙 곡선(누적 최대로 단조, 마지막 구간 선형 연장); 진폭 = Σy/Σg (무게 1/g 의 가중 최소제곱). 예측(뒤 절반): 379·413 K 상과 중앙 |K̂ − K| ≤ 10 K; '
                   'MAE 모양×진폭 ≤ 선형 보간 < 모집단. 대체 규칙: 모양×진폭 MAE(379·413 평균) > 선형 보간이면 주 곡선 = 선형 보간. 450 K 판: |편향| ≤ 10 K 이고 MAE ≤ 1.5 × 같은 도메인 413 K MAE. '
                   '뺀 것(가벼운 판): 표적 5 판들·반복 나눔·같은 반복, 절대 창 50–100 ns, K_σ, 바닥 폭 점검, Q 문턱.')
    jdump(out, RES + '/r39_E4.json'); log('R39_E4', V); log('R39_E4_DONE')


# ---------- E5 ----------
def dmu2_runs(A, Bm):
    cr = ((A[:, None] - Bm[None]) ** 2).sum(-1).mean(); fa = np.mean([((A[i] - A[k]) ** 2).sum() for i, k in itertools.combinations(range(len(A)), 2)])
    fb = np.mean([((Bm[i] - Bm[k]) ** 2).sum() for i, k in itertools.combinations(range(len(Bm)), 2)]); return float(cr - (fa + fb) / 2)


def splits(n):
    h = n // 2; return [(list(a), [k for k in range(n) if k not in a]) for a in itertools.combinations(range(n), h)] if n >= 4 else []


def null_max(P):
    sp = splits(len(P)); return max(dmu2_runs(P[a], P[b]) for a, b in sp) if sp else None


def num_cross(xb, tb, P, shared):
    if not shared: return float(np.dot(xb, tb))
    n = len(P); S = P.sum(0); full = np.dot(n * xb - S, n * tb - S); diag = sum(np.dot(xb - P[k], tb - P[k]) for k in range(n)); return float((full - diag) / (n * (n - 1)))


def sp_(a, b):
    return float(spearmanr(a, b).correlation)


def dir_pair(Xm, Xs, Rm, Rs, Tm, Ts, Pm, Ps, shared, nm):
    """X 의 이동(기준 R)과 T 의 이동(기준 P = MD 320)의 코사인·퍼짐 윤곽 상관. shared 면 R = P 이고 두 쪽에 서로 다른 320 K 반복을 쓴다."""
    o = {}; dX = dmu2_runs(Xm, Rm); dT = dmu2_runs(Tm, Pm); o['dmu_X'], o['dmu_T'] = dX, dT
    ok = nm is not None and dX > nm and dT > nm; o['protected'] = bool(ok)
    xb = Xm.mean(0) - (0 if shared else Rm.mean(0)); tb = Tm.mean(0) - (0 if shared else Pm.mean(0))
    o['cos'] = num_cross(xb, tb, Pm, shared) / np.sqrt(dX * dT) if ok else np.nan
    if shared:
        v = []
        for a, b in splits(len(Ps)):
            v.append(sp_(np.log(Xs.mean(0) / Ps[a].mean(0)), np.log(Ts.mean(0) / Ps[b].mean(0)))); v.append(sp_(np.log(Xs.mean(0) / Ps[b].mean(0)), np.log(Ts.mean(0) / Ps[a].mean(0))))
        o['prof'] = float(np.mean(v)) if v else np.nan
    else: o['prof'] = sp_(np.log(Xs.mean(0) / Rs.mean(0)), np.log(Ts.mean(0) / Ps.mean(0)))
    o['sign'] = float(np.mean(np.log(Xs.mean(0) / Rs.mean(0))))
    return o


def null_dir(Pm, Ps):
    if len(Pm) < 4: return np.nan, np.nan
    c, p = [], []
    for q in itertools.combinations(range(len(Pm)), 4):
        for (a, b), (e, f) in (((q[0], q[1]), (q[2], q[3])), ((q[0], q[2]), (q[1], q[3])), ((q[0], q[3]), (q[1], q[2]))):
            u, w = Pm[a] - Pm[b], Pm[e] - Pm[f]; c.append(np.dot(u, w) / np.linalg.norm(u) / np.linalg.norm(w)); p.append(sp_(np.log(Ps[a] / Ps[b]), np.log(Ps[e] / Ps[f])))
    return float(np.mean(c)), float(np.mean(p))


def e5_domain(D, Tver):
    """D: dict(X=(m, s), R=(m, s), shared, P=(m, s), MD={T: (m, s)}). Tver = 'top' | '379'."""
    MDT = {T: v for T, v in D['MD'].items() if len(v[0]) >= 2}
    T = ([T for T in ('450', '413', '379', '348') if T in MDT] or [None])[0] if Tver == 'top' else ('379' if '379' in MDT else None)
    if T is None or len(D['P'][0]) < 2: return None
    Pm, Ps = D['P']; nm = null_max(Pm); o = dict(T=T, null_max=nm)
    o.update(dir_pair(D['X'][0], D['X'][1], D['R'][0], D['R'][1], MDT[T][0], MDT[T][1], Pm, Ps, D['shared'], nm))
    o['heat_sign'] = float(np.mean(np.log(MDT[T][1].mean(0) / Ps.mean(0))))
    cx, ct = ('379', '413') if ('379' in MDT and '413' in MDT) else (('348', '379') if ('348' in MDT and '379' in MDT) else (None, None))
    if cx:
        c = dir_pair(MDT[cx][0], MDT[cx][1], Pm, Ps, MDT[ct][0], MDT[ct][1], Pm, Ps, True, nm); o['ceil_cos'] = c['cos']; o['ceil_prof'] = c['prof']; o['ceil_pair'] = f'{cx}·{ct}'
    else: o['ceil_cos'] = o['ceil_prof'] = np.nan
    o['null_cos'], o['null_prof'] = null_dir(Pm, Ps)
    return o


def e5_load_7c(model, doms):
    out = {}
    for d in doms:
        f = f'{DER}/r38/compare/{d}.npz'
        if not os.path.exists(f): continue
        z = np.load(f); k = f'{model}|0.75'
        if f'{k}|mean' not in z or 'MD320|full|0.75|mean' not in z: continue
        P = (z['MD320|full|0.75|mean'], z['MD320|full|0.75|spr'])
        out[d] = dict(dom=d, X=(z[f'{k}|mean'], z[f'{k}|spr']), R=P, shared=True, P=P, MD={T: (z[f'MD{T}|full|0.75|mean'], z[f'MD{T}|full|0.75|spr']) for T in TS if f'MD{T}|full|0.75|mean' in z})
    return out


def e5_load_prep(cell):
    out = {}
    for f in sorted(glob.glob(OUT39 + '/e5/*.json')):
        m = json.load(open(f))
        if m['cell'] != cell: continue
        z = np.load(f[:-5] + '.npz'); P = (z['MD320|mean'], z['MD320|spr'])
        out[m['tag']] = dict(dom=m['dom'], X=(z['X|mean'], z['X|spr']), R=(z['Xref|mean'], z['Xref|spr']), shared=bool(m['xref_is_md320']), P=P,
                             MD={T: (z[f'MD{T}|mean'], z[f'MD{T}|spr']) for T in TS if f'MD{T}|mean' in z}, meta=m)
    return out


def overlap(a, b):
    return not (a[1] < b[0] or b[1] < a[0])


def verdict5(x, ce, nu):
    if None in (x[0], ce[0], nu[0]): return '자료 부족'
    oc, on = overlap(x, ce), overlap(x, nu)
    if oc and not on: return '온도와 같은 방향'
    if oc and on: return '판정 불가'
    if on: return '직교에 가까움'
    if x[1] < ce[0]: return '온도와 반대 방향' if x[1] < 0 else '온도와 다른 방향'
    return '온도와 같은 방향(천장 위)'


def sign_verdict(c):
    return '0 과 구별 안 됨' if c[0] is None or c[0] <= 0 <= c[1] else ('가열과 같은 부호' if c[0] > 0 else '가열과 반대 부호')


def e5_cell(DD, sfmap, rng, Tver):
    R = [(t, D['dom'], e5_domain(D, Tver)) for t, D in DD.items()]; R = [(t, d, o) for t, d, o in R if o]
    if not R: return None
    o = dict(n_dom=len(R), T_used={T: sum(1 for *_, x in R if x['T'] == T) for T in TS})
    Pc = [x for x in R if x[2]['protected']]; o['n_dom_cos'] = len(Pc); o['n_excluded_cos'] = len(R) - len(Pc)
    if len(Pc) >= 3:
        u, g = sf_groups([sfmap.get(d, d) for _, d, _ in Pc]); pt, bs = boot_joint([np.array([x[k] for *_, x in Pc]) for k in ('cos', 'ceil_cos', 'null_cos')], g, len(u), rng)
        o['cos'], o['ceil_cos'], o['null_cos'] = summ(pt[0], bs[0]), summ(pt[1], bs[1]), summ(pt[2], bs[2]); o['n_sf_cos'] = len(u)
        o['verdict_cos'] = verdict5(o['cos']['ci95'], o['ceil_cos']['ci95'], o['null_cos']['ci95'])
    else: o['verdict_cos'] = '자료 부족(보호를 통과한 도메인 < 3)'
    u, g = sf_groups([sfmap.get(d, d) for _, d, _ in R]); ks = ('prof', 'ceil_prof', 'null_prof', 'sign', 'heat_sign'); pt, bs = boot_joint([np.array([x[k] for *_, x in R]) for k in ks], g, len(u), rng)
    for i, k in enumerate(ks): o[k] = summ(pt[i], bs[i])
    o['n_sf'] = len(u); o['verdict_prof'] = verdict5(o['prof']['ci95'], o['ceil_prof']['ci95'], o['null_prof']['ci95']); o['verdict_sign'] = sign_verdict(o['sign']['ci95'])
    o['domains'] = {t: {k: (v if not isinstance(v, float) or np.isfinite(v) else None) for k, v in x.items()} for t, _, x in R}
    return o


def part_e5():
    rng = np.random.default_rng(SEED['e5']); sfmap = cath_sf(); out = dict(written=time.strftime('%F %T'), cells={})
    p67 = open(DER + '/r38/p2_67.txt').read().split(); all7c = sorted(os.path.basename(f)[:-4] for f in glob.glob(DER + '/r38/compare/*.npz'))
    cells = {'R3′ 역장 + 규약 + ΔT 20 K (ATLAS 대 mdCATH 320 K)': e5_load_prep('R3p'), 'R3″ 규약만 (ONE 대 MSR)': e5_load_prep('R3pp'),
             '7c aSAMt320 대 MD 320 (P2 67)': e5_load_7c('aSAMt320', p67), '7c AlphaFlow 대 MD 320 (P2 67)': e5_load_7c('AlphaFlow', p67),
             '(민감도) 7c aSAMt320 (시험 전체)': e5_load_7c('aSAMt320', all7c), '(민감도) 7c AlphaFlow (시험 전체)': e5_load_7c('AlphaFlow', all7c)}
    for nm, DD in cells.items():
        for Tver in ('top', '379'):
            o = e5_cell(DD, sfmap, rng, Tver); out['cells'][f'{nm}|T={Tver}'] = o
            if o: log('E5', nm, Tver, o['n_dom'], o.get('n_dom_cos'), 'cos', (o.get('cos') or {}).get('est'), o['verdict_cos'], '윤곽', o['prof']['est'], o['verdict_prof'], '부호', o['sign']['est'], o['verdict_sign'])
    r = out['cells'].get('R3′ 역장 + 규약 + ΔT 20 K (ATLAS 대 mdCATH 320 K)|T=top')
    out['prediction'] = dict(R3p_sign_opposite=(r['verdict_sign'] == '가열과 반대 부호') if r else None, note='정상 점검(round37 의 퍼짐 비 −0.47 의 형식화). 코사인·윤곽 상관은 미리 정하지 않음')
    out['rule'] = ('round39 E5 주 판 (맞춘 좌표): 도메인마다 cos = mean_{i,j,k≠l}⟨m̄_X,i − m̄_320,k, m̄_T,j − m̄_320,l⟩ / √(Δμ²_X Δμ²_T) (기준이 MD 320 이 아니면 k ≠ l 조건 없음), '
                   'Δμ² 가 그 도메인 320 K 반복 나눔 귀무의 최댓값 안이면 뺌. 퍼짐 윤곽 = 변별 로그 퍼짐 변화의 Spearman (두 윤곽은 서로 다른 320 K 반복; 2:3 나눔 양방향 평균), 전체 부호 = 변 평균 로그 퍼짐 변화. '
                   'T = 접힌 프레임으로 N 을 채운 맨 위 온도(실행 ≥ 2) 와 379 K 판. 천장 = 379 대 413 (없으면 348 대 379), 귀무 = 320 K 반복 차 벡터끼리(서로 다른 넷). 상과 중앙값, 상과 부트스트랩 2,000 의 95%. '
                   '판정: 천장과 겹치고 귀무와 안 겹침 = 같은 방향; 상한 < 천장 하한이고 귀무와 안 겹침 = 다른 방향(상한 < 0 이면 반대); 귀무와 겹치고 천장과 안 겹침 = 직교에 가까움; 둘 다 겹침 = 판정 불가. '
                   '창: 전체(7c 모델 윤곽은 전체 창만, R3′ 는 d_ff 와 같은 프레임). R3″ 의 방향(ONE 대 MSR)은 물리적 부호가 없다.')
    jdump(out, RES + '/r39_E5.json'); log('R39_E5_DONE')


# ---------- pr ----------
def part_pr():
    src = open(ROOT + '/code/r38_s0.py').read(); src = src[:src.index("ALL = ['pick'")]; argv0 = sys.argv; sys.argv = ['r38_s0', ROOT, 'none']; S0 = dict(__name__='r38_s0_ns'); exec(src, S0); sys.argv = argv0
    dom_runs, md_load, fs, frames, QT, N, torch = (S0[k] for k in ('dom_runs', 'md_load', 'fs', 'frames', 'QT', 'N', 'torch'))
    fo = OUT39 + '/pr.json'; P = json.load(open(fo)) if os.path.exists(fo) else {}
    for d in sorted(os.path.basename(f)[:-4] for f in glob.glob(DER + '/mdcath_5T/*.npz')):
        if d in P: continue
        try:
            R, L = dom_runs(d, '320', 'half')
            if len(R) < 2: P[d] = dict(skip='실행 < 2'); continue
            MD, Qf, L = md_load(d); f = [X[Qf(X) >= QT] for X in MD['320']]; allf = np.concatenate([x for x in f if len(x)])[:600]; ref = torch.as_tensor(fs.mean_point(frames(allf)))
            A = fs.align_to_reference(frames(np.concatenate([X for _, X in R])), ref).numpy().reshape(len(R) * N, -1); A = A - A.mean(0)
            lam = np.linalg.eigvalsh(A @ A.T); lam = lam[lam > 1e-12 * lam.max()]; pr = float(lam.sum() ** 2 / (lam ** 2).sum())
            P[d] = dict(ne=L - 1, pr=pr, pr_frac=pr / (9 * (L - 1)), n_frames=len(A))
        except Exception as e:
            traceback.print_exc(); P[d] = dict(skip=repr(e)[:200])
        if len(P) % 20 == 0: jdump(P, fo, bk=False); log('pr', len(P))
    P = {d: v for d, v in P.items()}; jdump(P, fo, bk=False); log('R39_PR_DONE', sum('pr' in v for v in P.values()))


# ---------- check (합성) ----------
def part_check():
    rng = np.random.default_rng(SEED['check']); out = {}
    # 1 E0′: 가우스 앙상블에서 직접 비와 성분 비(유클리드, 반복 사이 평균 거리)
    def gauss(n, dl, a, b, R=5, N=100):
        d = 9 * n; mu = np.zeros(d); mu[0] = dl * np.sqrt(n)
        P = [rng.normal(0, a, (N, d)) for _ in range(R)]; Q = [rng.normal(0, b, (N, d)) + mu for _ in range(R)]
        md = lambda A, B: np.sqrt(np.maximum(((A[:, None] - B[None]) ** 2).sum(-1), 0)).mean()
        C = np.mean([md(x, y) for x in P for y in Q]); BP = np.mean([md(P[i], P[k]) for i, k in itertools.combinations(range(R), 2)]); BQ = np.mean([md(Q[i], Q[k]) for i, k in itertools.combinations(range(R), 2)])
        mP = np.stack([x.mean(0) for x in P]); mQ = np.stack([y.mean(0) for y in Q]); dm = dmu2_runs(mQ, mP) / n
        return e0_dom(C, BP, BQ, dm, n)
    t1 = [(cfg, gauss(*cfg)) for cfg in ((60, 0.05, 0.3, 0.3), (60, 0.6, 0.3, 0.3), (60, 0.0, 0.3, 0.45), (60, 0.2, 0.3, 0.33), (60, 3.0, 0.3, 0.3))]
    out['E0_two_ratios'] = [dict(cfg=c, r_dir=m['r_dir'], r_comp=m['r_comp'], rel=m['rel']) for c, m in t1]; out['E0_pass'] = bool(all(abs(m['rel']) < 0.1 for c, m in t1 if m['r_dir'] > 0.1))
    # 2 E1: 알려진 α·곱셈 잡음 + 덧셈 잡음(음수 꼬리)
    n = rng.integers(40, 500, 320).astype(float); ok1 = []
    for al in (0.6, 0.9, 1.2):
        y = 0.002 * n ** al * (1 + 0.3 * rng.normal(size=len(n))) + 0.002 * rng.normal(size=len(n)); a, gam, it = irls(n, y); ok1.append(dict(alpha_true=al, alpha_hat=a, gamma=gam, iters=it))
    out['E1_irls'] = ok1; out['E1_pass'] = bool(all(abs(x['alpha_hat'] - x['alpha_true']) < 0.08 for x in ok1))
    # 3 E3: 효과 없음(T = R0 나눔 하나)에서 검출률, 큰 효과에서 검출률
    ns = 60; mem = [np.arange(i * 3, i * 3 + 3) for i in range(ns)]; r0 = rng.normal(0, 0.002, (ns * 3, 10)); v0 = r0[:, 0] + 0; v1 = r0[:, 0] + 0.01; R = []
    for v, nm in ((v0, '효과 없음'), (v1, '큰 효과')):
        mT = np.array([np.median(v[m]) for m in mem]); sch = lambda: rng.integers(0, ns, 30); band = r0_band(r0, mem, rng, 500, sch); det = []
        for _ in range(200):
            S = sch(); lb = np.percentile(np.median(mT[S][rng.integers(0, len(S), (200, len(S)))], axis=1), 2.5); det.append(lb > band[1])
        R.append((nm, float(np.mean(det))))
    out['E3_detect'] = R; out['E3_pass'] = bool(R[0][1] <= 0.05 and R[1][1] >= 0.95)
    # 4 E4: 무잡음 곡선에서 역환산
    xs = [0, 28, 93, 130]; gs = [0, 0.002, 0.012, 0.03]; inv = pw(xs, gs); k59 = inv(float(np.interp(59, xs, gs)))[0]; k150 = inv(0.03 + (0.03 - 0.012) / 37 * 20)[0]
    out['E4_inverse'] = dict(K59=k59, K150=k150); out['E4_pass'] = bool(abs(k59 - 59) < 1e-9 and abs(k150 - 150) < 1e-6)
    # 5 E5: 알려진 각도의 코사인(공유 320 K 잡음, 실행 평균 수준)
    d = 900; u = rng.normal(size=d); u /= np.linalg.norm(u); w_ = rng.normal(size=d); w_ -= w_.dot(u) * u; w_ /= np.linalg.norm(w_); R5 = []
    for th in (0, 60, 90, 180):
        cs = []
        for _ in range(60):
            mu0 = rng.normal(size=d); sig = 0.02; P =mu0 + sig * rng.normal(size=(5, d)); T = mu0 + 1.0 * u + sig * rng.normal(size=(5, d))
            X = mu0 + 0.8 * (np.cos(np.radians(th)) * u + np.sin(np.radians(th)) * w_) + sig * rng.normal(size=(2, d))
            dX, dT = dmu2_runs(X, P), dmu2_runs(T, P); cs.append(num_cross(X.mean(0), T.mean(0), P, True) / np.sqrt(dX * dT))
        R5.append(dict(theta=th, cos_true=float(np.cos(np.radians(th))), cos_mean=float(np.mean(cs)), cos_sd=float(np.std(cs))))
    out['E5_cos'] = R5; out['E5_pass'] = bool(all(abs(x['cos_mean'] - x['cos_true']) < 0.05 for x in R5))
    out['all_pass'] = bool(all(out[k] for k in out if k.endswith('_pass'))); out['written'] = time.strftime('%F %T')
    jdump(out, RES + '/r39_check.json'); log('R39_CHECK', {k: out[k] for k in out if k.endswith('_pass')}, '모두', out['all_pass']); log('R39_CHECK_DONE')


if __name__ == '__main__':
    for part in PARTS39:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
