"""(탐색 — 결과를 본 뒤; 등록 판정 아님) 1차 P3 의 등록 구간 방법(상과 × 두 도메인 Exp(1) 무게의 곱)이 내포 구조에서 맞는 크기인지 보정 시뮬레이션으로 본다.
results/r38_p3_explore2.json.
(1) 실제 자료: 공동 관문 쌍 집합(주 판·교집합 판)과 민감도 판들을 상과만 Exp(1) 부트스트랩으로 다시 낸다(점추정은 같음).
(2) 보정 시뮬레이션: 실제 쌍·도메인·상과 구조와 무게(1/상과 쌍 수)를 그대로 두고, 세 변수(r_dyn·동일성·통제의 순위, 표준화)를
    V_p = s_상과 + d_a + d_b + e_p 꼴로 만든다. 분산 성분은 실제 자료에서 어림: 같은 상과·도메인 안 겹침 쌍끼리 공분산 = σ_s²,
    한 도메인을 나눈 쌍끼리 공분산 = σ_s² + σ_d². 귀무: Y 와 X 는 Z 를 통해서만 상관(실제 순위 상관 크기), Y ⟂ X | Z.
    반복마다 점추정과 두 부트스트랩(등록 판, 상과만; 각 400)의 SE·95% 구간. 정답 = 반복 사이 점추정의 SD, 피복률 = 구간에 0.
사용: CUDA_VISIBLE_DEVICES= python r38_p3_explore2.py <ROOT>"""
import sys, json, time
ROOT = sys.argv[1]
P3 = {'__name__': 'p3x2'}; _s = open(ROOT + '/code/r38_p3.py').read(); _s = _s[:_s.index("ALL = ['pairs'")]
sys.argv = ['r38_p3', ROOT, 'none', '1']; exec(_s, P3)
np = P3['np']; ranks, pcorr, wcorr = P3['ranks'], P3['pcorr'], P3['wcorr']
out = dict(written=time.strftime('%F %T'), note='탐색(결과를 본 뒤). 등록 판정(r38_p3.json)을 바꾸지 않는다')
rows_all = [json.loads(l) for l in open(P3['OUT'] + '/pairs_out.jsonl')]


def setup(rows, key='main', y='r_dyn', x='ident', z=('control',)):
    rows = [r for r in rows if r[key].get(y) is not None]
    Y = ranks(np.array([r[key][y] for r in rows])); X = ranks(np.array([r[key][x] for r in rows])); Z = [ranks(np.array([r[key][c] if c in r[key] else r[c] for r in rows])) for c in z]
    H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); cnt = np.bincount(si); base = 1.0 / cnt[si]
    doms = sorted({r['ab'][0] for r in rows} | {r['ab'][1] for r in rows}); di = {d: i for i, d in enumerate(doms)}
    ia = np.array([di[r['ab'][0]] for r in rows]); ib = np.array([di[r['ab'][1]] for r in rows])
    return Y, X, Z, si, base, ia, ib, len(sfl), len(doms)


def boots(Y, X, Z, si, base, ia, ib, nsf, ndom, rng, B):
    reg, sf = [], []
    for _ in range(B):
        v = rng.exponential(1, nsf); w = rng.exponential(1, ndom)
        reg.append(pcorr(Y, X, Z, base * v[si] * w[ia] * w[ib])); sf.append(pcorr(Y, X, Z, base * v[si]))
    return np.array(reg), np.array(sf)


ci = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
# (1) 실제 자료
rng = np.random.default_rng(20261051); real = {}
both = [r for r in rows_all if r['main'].get('r_dyn') is not None and r['inter'].get('r_dyn') is not None]
for name, rows, kw in (('main', rows_all, {}), ('gate_main', both, {}), ('gate_inter', both, dict(key='inter')), ('control_meanpt', rows_all, dict(z=('control_meanpt',))),
                       ('control_tm', rows_all, dict(z=('tm_long',))), ('with_degen', rows_all, dict(y='r_dyn_withdegen', z=('control_withdegen',))), ('shared_frac_cov', rows_all, dict(z=('control', 'shared_frac')))):
    S = setup(rows, **kw); est = pcorr(S[0], S[1], S[2], S[4]); reg, sf = boots(*S, rng, 1000)
    real[name] = dict(n=len(S[0]), est=float(est), ci_reg=ci(reg), ci_sf=ci(sf), se_reg=float(reg.std()), se_sf=float(sf.std()))
    print(name, real[name], flush=True)
out['real_sf_only'] = real
# (2) 보정 시뮬레이션
Y, X, Z, si, base, ia, ib, nsf, ndom = setup(rows_all); z0 = Z[0]
std = lambda v: (v - v.mean()) / v.std()
Ys, Xs, Zs = std(Y), std(X), std(z0)
n = len(Ys); rng = np.random.default_rng(20261052)
# 겹침 쌍 표본: 같은 상과·도메인 안 겹침 / 한 도메인 공유
by_sf = {}
for p in range(n): by_sf.setdefault(si[p], []).append(p)
dis, shr = [], []
for s_, ps in by_sf.items():
    if len(ps) < 2: continue
    ps = np.array(ps); k = min(len(ps) * 4, 4000)
    a_ = rng.choice(ps, k); b_ = rng.choice(ps, k); m_ = a_ != b_; a_, b_ = a_[m_], b_[m_]
    sh_ = (ia[a_] == ia[b_]) | (ia[a_] == ib[b_]) | (ib[a_] == ia[b_]) | (ib[a_] == ib[b_]); dis += list(zip(a_[~sh_], b_[~sh_])); shr += list(zip(a_[sh_], b_[sh_]))
dis, shr = np.array(dis), np.array(shr)
def comps(V):
    cs = float(np.mean(V[dis[:, 0]] * V[dis[:, 1]])); cd = float(np.mean(V[shr[:, 0]] * V[shr[:, 1]])); s2 = max(cs, 0.0); d2 = max(cd - cs, 0.0); e2 = max(1 - s2 - 2 * d2, 0.05)
    return s2, d2, e2
cY, cX, cZ = comps(Ys), comps(Xs), comps(Zs); rxz = float(wcorr(Xs, Zs, base)); ryz = float(wcorr(Ys, Zs, base))
out['sim_setup'] = dict(n_pairs=n, n_sf=nsf, n_dom=ndom, n_disjoint_pairs=len(dis), n_share_pairs=len(shr), comps_Y=cY, comps_X=cX, comps_Z=cZ, r_XZ=rxz, r_YZ=ryz)
print(out['sim_setup'], flush=True)
def field(c, rng):
    s2, d2, e2 = c; return np.sqrt(s2) * rng.normal(size=nsf)[si] + np.sqrt(d2) * (rng.normal(size=ndom)[ia] + rng.normal(size=ndom)[ib]) + np.sqrt(e2) * rng.normal(size=n)
REP, BB = 150, 400; ests, se_reg, se_sf, cov_reg, cov_sf = [], [], [], [], []
for rep in range(REP):
    Zf = field(cZ, rng); Zf = (Zf - Zf.mean()) / Zf.std()
    Xf = rxz * Zf + np.sqrt(max(1 - rxz ** 2, 0)) * field(cX, rng); Yf = ryz * Zf + np.sqrt(max(1 - ryz ** 2, 0)) * field(cY, rng)
    Yr, Xr, Zr = ranks(Yf), ranks(Xf), [ranks(Zf)]; e = pcorr(Yr, Xr, Zr, base); reg, sf = boots(Yr, Xr, Zr, si, base, ia, ib, nsf, ndom, rng, BB)
    ests.append(e); se_reg.append(reg.std()); se_sf.append(sf.std()); a1, b1 = ci(reg); a2, b2 = ci(sf); cov_reg.append(a1 <= 0 <= b1); cov_sf.append(a2 <= 0 <= b2)
    if rep % 25 == 0: print('sim', rep, round(float(np.std(ests)) if len(ests) > 1 else 0, 4), round(float(np.mean(se_reg)), 4), round(float(np.mean(se_sf)), 4), flush=True)
out['sim'] = dict(reps=REP, boot=BB, true_sd_est=float(np.std(ests)), mean_est=float(np.mean(ests)), mean_se_reg=float(np.mean(se_reg)), mean_se_sf=float(np.mean(se_sf)),
                  coverage_reg=float(np.mean(cov_reg)), coverage_sf=float(np.mean(cov_sf)),
                  read='정답(반복 사이 점추정 SD)에 가까운 SE 를 주고 피복률이 0.95 근처인 방법이 맞는 크기. 등록 판 SE ≫ 정답 SD 면 등록 판이 과하게 보수적')
json.dump(out, open(ROOT + '/results/r38_p3_explore2.json', 'w'), ensure_ascii=False, indent=1, default=float)
print(json.dumps(out['sim'], ensure_ascii=False)); print('P3X2_DONE')
