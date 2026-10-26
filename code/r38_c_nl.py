"""(결과 뒤 등록한 강건성 점검 — 확인 실험 C 판정(2026-10-20 14:43 UTC) 뒤, 계산 전에 Drive JSON round36_c_nl_2026-10-20 에 등록; 사용자 "그렇게 해봐")
물음: C-1·C-2 의 편 Spearman 은 통제를 순위에서 선형으로 뺀다. 통제와 결과·동일성 사이가 비선형이면 남은 효과의 일부가 비선형 교란일 수 있다.
(가) 주: 유연한 통제(이중 기계학습 꼴 편상관). Y = 순위(r_dyn), X = 순위(동일성), Z = 구조·크기·정렬 변수 8 개
     (control, control_meanpt, r_anm, meanpt_D_scaled, log n_shared, shared_frac, degen_frac, log L_long).
     E[Y|Z]·E[X|Z] 를 HistGradientBoosting(초매개변수 고정, 조기 멈춤 없음, 쌍 무게 = 1/상과 쌍 수)으로 상과 단위 5 겹 교차 적합 → 잔차의 쌍 무게 상관.
     구간 = 상과 Exp(1) 군집 부트스트랩 2,000 (잔차 고정). 같은 행에서 견줌: 등록 C-2 꼴(선형, control + r_anm), 선형 + Z 8 개.
     고침 1 (실제 자료 전, 합성 점검이 실패해서): 처음 판(300 회, 잎 15, Z 8 개)은 합성 자료의 곱 상호작용을 못 잡아 −0.20 치우침.
     주 판 = 특성 44 개(표준화 Z, 제곱, 두 변수 곱) + 1000 회·잎 31·최소 잎 20. 합성 치우침 −0.04 (치우침 예산으로 보고). 처음 판은 nl_plain 으로 함께 낸다.
(나) 모양: 비선형 잔차 rY 의 동일성 십분위별 무게 평균(상과 부트스트랩 95%) — 단조·문턱·포화를 본다(서술).
(다) 일반 의존: 잔차 (rY, rX) 의 거리상관, 상과 안 순열(rX 를 상과 안에서 섞음) 1,000 회의 p (쌍 ≤ 4,000 은 전부, 넘으면 고정 seed 표본 4,000).
자료: C 주 판(판정), C 교집합 판·H-풀 주 판(서술).
판정 문장(계산 전 고정; C 주 판의 비선형 판으로):
  하한 > 0 → "구조·ANM·크기의 비선형 통제 뒤에도 서열 동일성의 효과가 남는다(선형 판 대비 [비])"
  하한 ≤ 0 이고 상한 < 0.15 → "동역학 닮음의 서열 효과는 구조 변수의 비선형 함수로 설명된다"
  그 밖 → "판정 불가(비선형 통제 뒤 구간이 넓음)"
check: 합성 자료(실제 값 안 봄; 상과 구조만 C 와 같은 꼴로 흉내) — (1) Y, X 가 Z 의 비선형 함수 + 잡음이고 직접 효과 0 이면 선형 판은 0 에서 벗어나고 비선형 판은 0 근처,
       (2) 직접 효과를 심으면 비선형 판이 되찾음.
사용: python r38_c_nl.py <ROOT> check,run 출력: results/r38_c_nl.json"""
import sys, os, json, time
import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold
ROOT = sys.argv[1]; PARTS = sys.argv[2].split(','); RES = ROOT + '/results'; FO = RES + '/r38_c_nl.json'
ZK = ('control', 'control_meanpt', 'r_anm', 'meanpt_D_scaled', 'log_n_shared', 'shared_frac', 'degen_frac', 'log_L_long')
HP = dict(max_iter=1000, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=20, early_stopping=False, random_state=20261051) # 주 판 (합성 점검으로 고름, 고침 1)
HP_PLAIN = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=40, early_stopping=False, random_state=20261051) # 처음 등록 판(용량 민감도)


def augment(Z):
    """표준화한 Z, 그 제곱, 두 변수 곱(8 + 8 + 28 = 44). 곱 상호작용은 나무가 한 변수씩 나눠서는 찾기 어렵다(합성 점검)."""
    S = (Z - Z.mean(0)) / Z.std(0); return np.column_stack([S, S ** 2] + [S[:, i] * S[:, j] for i in range(S.shape[1]) for j in range(i + 1, S.shape[1])])
OUT = json.load(open(FO)) if os.path.exists(FO) else {}
def log(*a): print(time.strftime('%H:%M:%S'), *a, flush=True)


# r38_p3 와 같은 정의
def ranks(v): return np.argsort(np.argsort(v)).astype(float)
def wls_resid(y, X, w):
    Xd = np.column_stack([np.ones(len(y))] + list(X)); W = np.sqrt(w); beta = np.linalg.lstsq(Xd * W[:, None], y * W, rcond=None)[0]; return y - Xd @ beta
def wcorr(a, b, w):
    w = w / w.sum(); ma, mb = (w * a).sum(), (w * b).sum(); return float((w * (a - ma) * (b - mb)).sum() / np.sqrt((w * (a - ma) ** 2).sum() * (w * (b - mb) ** 2).sum()))
def pcorr(ry, rx, Z, w): return wcorr(wls_resid(ry, Z, w), wls_resid(rx, Z, w), w)


def load(outdir, anmf, tag='main'):
    an = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(anmf)}; rows = []
    for l in open(outdir + '/pairs_out.jsonl'):
        r = json.loads(l); m = r.get(tag, {}); a = an.get(tuple(r['ab']), {}); am = a.get(tag, {})
        if m.get('r_dyn') is None or m.get('ident') is None or m.get('control') is None or am.get('r_anm') is None: continue
        z = dict(control=m['control'], control_meanpt=m.get('control_meanpt'), r_anm=am['r_anm'], meanpt_D_scaled=r.get('meanpt_D_scaled'),
                 log_n_shared=float(np.log(m['n_shared'])), shared_frac=m.get('shared_frac'), degen_frac=m.get('degen_frac'), log_L_long=float(np.log(a['L_long'])) if a.get('L_long') else None)
        if any(z[k] is None or not np.isfinite(z[k]) for k in ZK): continue
        rows.append(dict(H=r['H'], y=m['r_dyn'], x=m['ident'], **z))
    return rows


def nl_resid(Y, X, Z, si, base, hp=None, aug=True):
    rY, rX = np.zeros(len(Y)), np.zeros(len(X)); F = augment(Z) if aug else Z
    for tr, te in GroupKFold(5).split(F, Y, groups=si):
        for tgt, res in ((Y, rY), (X, rX)):
            m = HistGradientBoostingRegressor(**(hp or HP)).fit(F[tr], tgt[tr], sample_weight=base[tr]); res[te] = tgt[te] - m.predict(F[te])
    return rY, rX


def dcor_perm(a, b, groups, rng, n_perm=1000, nmax=4000):
    if len(a) > nmax:
        ix = np.sort(rng.choice(len(a), nmax, replace=False)); a, b, groups = a[ix], b[ix], groups[ix]
    def dc(v):
        D = np.abs(v[:, None] - v[None]); return D - D.mean(0)[None] - D.mean(1)[:, None] + D.mean()
    A, Bm = dc(a), dc(b); vA, vB = (A * A).mean(), (Bm * Bm).mean(); obs = float(np.sqrt(max((A * Bm).mean(), 0) / np.sqrt(vA * vB)))
    gl = {}
    for i, g in enumerate(groups): gl.setdefault(g, []).append(i)
    gl = [np.array(v) for v in gl.values()]; ge = 0
    for _ in range(n_perm):
        p = np.arange(len(b))
        for v in gl: p[v] = rng.permutation(v)
        s = float(np.sqrt(max((A * Bm[p][:, p]).mean(), 0) / np.sqrt(vA * vB))); ge += s >= obs
    return dict(n=len(a), dcor=obs, p_within_sf=(ge + 1) / (n_perm + 1))


def analyse(rows, seed, B=2000, shape=True, dcor=True):
    H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); cnt = np.bincount(si); base = 1.0 / cnt[si]
    Y = ranks(np.array([r['y'] for r in rows])); X = ranks(np.array([r['x'] for r in rows])); Z = np.array([[r[k] for k in ZK] for r in rows], float); Zr = [ranks(Z[:, j]) for j in range(Z.shape[1])]
    Z2 = [Zr[ZK.index('control')], Zr[ZK.index('r_anm')]]
    est = dict(lin2=pcorr(Y, X, Z2, base), linall=pcorr(Y, X, Zr, base)); rY, rX = nl_resid(Y, X, Z, si, base); est['nl'] = wcorr(rY, rX, base)
    pY, pX = nl_resid(Y, X, Z, si, base, HP_PLAIN, aug=False); est['nl_plain'] = wcorr(pY, pX, base)
    rng = np.random.default_rng(seed); bs = {k: [] for k in est}; xr = np.array([r['x'] for r in rows]); qs = np.quantile(xr, np.linspace(0, 1, 11)); bins = np.clip(np.searchsorted(qs, xr, side='right') - 1, 0, 9); sh = [[] for _ in range(10)]
    for _ in range(B):
        w = base * rng.exponential(1, len(sfl))[si]
        bs['lin2'].append(pcorr(Y, X, Z2, w)); bs['linall'].append(pcorr(Y, X, Zr, w)); bs['nl'].append(wcorr(rY, rX, w)); bs['nl_plain'].append(wcorr(pY, pX, w))
        if shape:
            for k in range(10): sh[k].append(float((w[bins == k] * rY[bins == k]).sum() / w[bins == k].sum()))
    out = dict(n=len(rows), n_sf=len(sfl), est=est, ci95={k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] for k, v in bs.items()},
               r2_Y_on_Z=float(1 - np.average(rY ** 2, weights=base) / np.average((Y - np.average(Y, weights=base)) ** 2, weights=base)),
               r2_X_on_Z=float(1 - np.average(rX ** 2, weights=base) / np.average((X - np.average(X, weights=base)) ** 2, weights=base)))
    if shape:
        out['shape'] = [dict(ident_lo=float(qs[k]), ident_hi=float(qs[k + 1]), n=int((bins == k).sum()), mean_rY=float((base[bins == k] * rY[bins == k]).sum() / base[bins == k].sum()),
                             ci95=[float(np.percentile(sh[k], 2.5)), float(np.percentile(sh[k], 97.5))]) for k in range(10)]
    if dcor: out['dcor'] = dcor_perm(rY, rX, si, np.random.default_rng(seed + 1))
    return out


def verdict(a):
    lo, hi = a['ci95']['nl']
    if lo > 0: return f"구조·ANM·크기의 비선형 통제 뒤에도 서열 동일성의 효과가 남는다(선형 판 대비 {a['est']['nl'] / a['est']['lin2']:.2f})"
    if hi < 0.15: return '동역학 닮음의 서열 효과는 구조 변수의 비선형 함수로 설명된다'
    return '판정 불가(비선형 통제 뒤 구간이 넓음)'


def part_check():
    """합성 자료: 상과 구조만 C 와 같은 꼴(상과 242, 쌍 4,152 의 크기 분포)로 흉내, 값은 모두 합성."""
    Hc = [json.loads(l)['H'] for l in open(ROOT + '/derived/r38/c/pairs_out.jsonl')]; rng = np.random.default_rng(20261054); res = {}
    sf = sorted(set(Hc)); u = {h: rng.normal(size=8) for h in sf}
    for case, beta in (('no_direct', 0.0), ('direct_0.3', 0.3)):
        rows = []
        for h in Hc:
            z = 0.7 * u[h] + 0.7 * rng.normal(size=8)
            fz = np.sin(2 * z[0]) + z[1] ** 2 + z[2] * z[3] + np.abs(z[4]); gz = np.cos(2 * z[0]) + z[1] ** 2 - z[2] * z[3] + 0.5 * np.abs(z[4])
            x = gz + 0.8 * rng.normal(); y = fz + beta * x + 0.8 * rng.normal()
            rows.append(dict(H=h, y=y, x=x, **{k: float(z[i]) for i, k in enumerate(ZK)}))
        a = analyse(rows, 20261055, B=300, shape=False, dcor=False); res[case] = dict(est=a['est'], ci95=a['ci95']); log('check', case, a['est'], a['ci95'])
    ok1 = res['no_direct']['ci95']['nl'][0] <= 0 <= res['no_direct']['ci95']['nl'][1] and not (res['no_direct']['ci95']['lin2'][0] <= 0 <= res['no_direct']['ci95']['lin2'][1])
    ok2 = res['direct_0.3']['ci95']['nl'][0] > 0
    OUT['check'] = dict(res, pass_no_direct=bool(ok1), pass_direct=bool(ok2), note='합성: 직접 효과 0 이면 선형 판(control·r_anm 꼴 두 변수)이 거짓 양성이고 비선형 판은 0 을 포함해야 함; 직접 효과 0.3 이면 비선형 판 하한 > 0')
    log('NL_CHECK', ok1, ok2)


def part_run():
    C = ROOT + '/derived/r38/c'; P3 = ROOT + '/derived/r38/p3'
    OUT['C_main'] = analyse(load(C, C + '/anm_pairs.jsonl', 'main'), 20261052); OUT['C_main']['verdict'] = verdict(OUT['C_main']); log('C_main', OUT['C_main']['est'], OUT['C_main']['ci95'], OUT['C_main']['verdict'])
    OUT['C_inter'] = analyse(load(C, C + '/anm_pairs.jsonl', 'inter'), 20261056, dcor=False); log('C_inter', OUT['C_inter']['est'], OUT['C_inter']['ci95'])
    OUT['H_main'] = analyse(load(P3, C + '/anm_pairs_hpool.jsonl', 'main'), 20261053); log('H_main', OUT['H_main']['est'], OUT['H_main']['ci95'])
    OUT['rule'] = '결과 뒤 등록(계산 전 Drive JSON round36_c_nl_2026-10-20). 판정 = C 주 판 비선형 판: 하한 > 0 → 남는다; 하한 ≤ 0 이고 상한 < 0.15 → 비선형 구조로 설명됨; 그 밖 판정 불가'


if __name__ == '__main__':
    for part in PARTS:
        log('부분 시작', part); globals()['part_' + part](); OUT['written'] = time.strftime('%F %T')
        json.dump(OUT, open(FO + '.tmp', 'w'), ensure_ascii=False, indent=1); os.replace(FO + '.tmp', FO); log('부분 끝', part)
    log('C_NL_DONE')
