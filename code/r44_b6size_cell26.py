"""(round44, 등록 = results/prereg_2026-10-28_addendum3.json "r44_b6size", 2026-10-28 06:04:10 UTC) 연속 축(B6) 크기 점검.
이 계산은 Colab 노트북 셀 26 에서 돌았다(2026-10-28 06:04:37 UTC). 아래 "셀 26 그대로" 줄부터 끝까지가 그 셀의 글자 그대로이고,
그 위의 머리(모듈 불러오기, R = 프로젝트 루트)는 이 파일만으로 돌리려고 저장소를 꾸릴 때 더한 것이다. 앞 셀에서 jdump 가 정의되지 않았으면 json.dump 로 쓴다.
사용: python r44_b6size_cell26.py <ROOT>   입력: derived/r38/b6/*.json   출력: results/r44_b6size.json"""
import os, sys, json, time
import numpy as np
R = sys.argv[1]
# ---- 셀 26 그대로 ----
# 26. r44_b6size (등록 셋째, 06:04:10 UTC): 편 Spearman(D.d2n, TM | log nbar), 도메인 Exp(1) 곱 무게 2,000
import glob
rows = []
for f in sorted(glob.glob(f'{R}/derived/r38/b6/*.json')):
    z = json.load(open(f)); Dz = z.get('D') or {}
    if Dz.get('d2n') is not None and z.get('TM') is not None and z.get('nbar'):
        rows.append((z['a'], z['b'], float(Dz['d2n']), float(z['TM']), float(z['nbar'])))
print('쌍', len(rows))
def rk(v): return np.argsort(np.argsort(v)).astype(float)
Y = rk(np.array([r[2] for r in rows])); X = rk(np.array([r[3] for r in rows])); Z = rk(np.log([r[4] for r in rows]))
doms = sorted({x for r in rows for x in r[:2]}); ia = np.array([doms.index(r[0]) for r in rows]); ib = np.array([doms.index(r[1]) for r in rows])
def wcor(u, v, w):
    w = w / w.sum(); u = u - (w * u).sum(); v = v - (w * v).sum(); return float((w * u * v).sum() / np.sqrt((w * u * u).sum() * (w * v * v).sum()))
def wres(u, z, w):
    A = np.c_[np.ones_like(z), z]; W = np.sqrt(w); b = np.linalg.lstsq(A * W[:, None], u * W, rcond=None)[0]; return u - A @ b
def pc(w): return wcor(wres(Y, Z, w), wres(X, Z, w), w)
def uc(w): return wcor(Y, X, w)
w1 = np.ones(len(rows)); est = dict(partial=pc(w1), plain=uc(w1)); rng = np.random.default_rng(20261098); bs = {'partial': [], 'plain': []}
for _ in range(2000):
    v = rng.exponential(1, len(doms)); w = v[ia] * v[ib]; bs['partial'].append(pc(w)); bs['plain'].append(uc(w))
out44 = dict(written=time.strftime('%F %T UTC', time.gmtime()), n_pairs=len(rows), n_dom=len(doms),
             partial=dict(est=est['partial'], lo=float(np.percentile(bs['partial'], 2.5)), hi=float(np.percentile(bs['partial'], 97.5))),
             plain=dict(est=est['plain'], lo=float(np.percentile(bs['plain'], 2.5)), hi=float(np.percentile(bs['plain'], 97.5))),
             rho_logn_d2n=uc(w1) if False else wcor(Y, Z, w1), rho_logn_TM=wcor(X, Z, w1))
out44['prediction_holds'] = bool(out44['partial']['est'] <= -0.6 and out44['partial']['hi'] < 0)
out44['rule'] = '등록 r44_b6size (prereg_2026-10-28_addendum3.json): 예측 = 편상관 ≤ −0.6 이고 97.5% 상한 < 0.'
jdump(out44, f'{R}/results/r44_b6size.json', bk=True) if 'jdump' in globals() else json.dump(out44, open(f'{R}/results/r44_b6size.json', 'w'), ensure_ascii=False, indent=1)
print(json.dumps(out44, ensure_ascii=False, indent=1))
