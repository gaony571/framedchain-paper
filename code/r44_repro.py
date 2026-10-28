"""(2026-10-28 UTC, 원고 방법 절의 "0.4.1 이 분석에 쓴 코드를 재현하는가" 점검; 외부 검토 2026-10-28 의 10 번) 결과 자료는 쓰지 않는다.
분석은 framespace 0.2.0 + 연구 모듈 framespace_wip/stage0 + 스크립트(code/r38_s0.py 의 red_pairs·estats)로 돌았다. 셋 모두 재현 저장소 framedchain-paper 의 실행본(등록 md5)을 쓴다.
같은 합성 입력에서 framedchain 0.4.1 이 이들과 같은 값을 내는지 본다(같은 seed, CPU, float64):
  거리·정렬(distance, align), 인증(bnb), 평균 점(mean_point), 앙상블 W₂(ensemble_w2) — 0.2.0 대 0.4.1
  줄인 계산: stage0.band 의 MeanStart·band_windows·align_band 와 r38_s0.red_pairs 대 framedchain.reduced 의 같은 이름·pair_scores
  에너지 통계: r38_s0.estats 대 framedchain.ensemble.energy_stats
사용: python r44_repro.py <atlas 경로>   출력: 표준 출력(최대 절대 차)"""
import sys, os, math
A = sys.argv[1]
sys.path[:0] = [os.path.join(A, 'framedchain-paper', 'framespace_pkg'), os.path.join(A, 'framedchain-paper', 'framespace_wip'), os.path.join(A, 'framedchain')]
import numpy as np, torch
import framespace as f0
import framedchain as fc
from framedchain import reduced as rd, ensemble as en
import stage0.band as sb
torch.set_default_dtype(torch.float64)
print('framespace', getattr(f0, '__version__', '?'), f0.__file__.split(os.sep)[-3], '| framedchain', fc.__version__)

g = torch.Generator().manual_seed(20261095)
def chain(n):
    step = fc.expmap(torch.tensor([0.3, 1.2, 0.4]))
    R = [torch.eye(3)]
    for _ in range(n - 1): R.append(R[-1] @ step @ fc.expmap(0.2 * torch.randn(3, generator=g)))
    return torch.stack(R) * (1 + 0.02 * torch.randn(n, 1, 1, generator=g)) / 2
def ens(base, F, sigma):
    W = fc.expmap(sigma * torch.randn(F, len(base), 3, generator=g))
    return torch.einsum('nxy,fnyz->fnxz', base, W)
worst = {}
def rec(k, v): worst[k] = max(worst.get(k, 0.0), float(v))

for _ in range(15):
    x, y = chain(int(torch.randint(15, 40, (1,), generator=g))), chain(int(torch.randint(15, 40, (1,), generator=g)))
    d0, a0 = f0.distance(x, y); d1, a1 = fc.distance(x, y)
    rec('distance', abs(float(d0) - float(d1))); rec('align score', abs(float(a0.score) - float(a1.score)))
    rec('align match', float((a0.match != a1.match).sum()))
for _ in range(3):
    x, y = chain(20), chain(24)
    b0, b1 = f0.bnb(x, y), fc.bnb(x, y)
    rec('bnb S_lb/S_ub', max(abs(b0['S_lb'] - b1['S_lb']), abs(b0['S_ub'] - b1['S_ub'])))
E = ens(chain(25), 30, 0.08); rec('mean_point', float((f0.mean_point(E) - fc.mean_point(E)).abs().max()))
EA, EB = ens(chain(18), 8, 0.08), ens(chain(21), 8, 0.08)
w0, w1 = f0.ensemble_w2(EA, EB), fc.ensemble_w2(EA, EB); rec('ensemble_w2', abs(float(w0['W2']) - float(w1['W2'])))

# 줄인 계산: 실행 3 × 프레임 20, 두 앙상블
EA = torch.cat([ens(chain(30), 20, 0.07) for _ in range(3)]); EB = torch.cat([ens(chain(34), 20, 0.07) for _ in range(3)])
m0, m1 = sb.MeanStart(EA, EB), rd.MeanStart(EA, EB)
rec('MeanStart center', float(np.abs(np.asarray(m0.center) - np.asarray(m1.center)).max()))
ia = torch.randint(0, len(EA), (500,), generator=g); ib = torch.randint(0, len(EB), (500,), generator=g)
rec('MeanStart starts', float((m0.starts(ia, ib, n_local=1) - m1.starts(ia, ib, n_local=1)).abs().max()))
src = open(os.path.join(A, 'framedchain-paper', 'code', 'r38_s0.py'), encoding='utf-8').read()
def grab(name):
    i = src.index(f'def {name}('); j = src.index('\n\n\n', i); return src[i:j]
ns = dict(np=np, torch=torch, band_windows=sb.band_windows, align_band=sb.align_band, sq_norm=f0.distance.sq_norm if hasattr(f0.distance, 'sq_norm') else __import__('framespace.distance', fromlist=['sq_norm']).sq_norm, K0=5)
exec(grab('red_pairs'), ns); exec(grab('estats'), ns)
r0 = ns['red_pairs'](m0, ia.numpy(), ib.numpy(), want_gap=True); r1 = rd.pair_scores(m1, ia.numpy(), ib.numpy(), want_gap=True)
d1 = np.asarray(r1['d']); g1 = np.asarray(r1['gap'])
rec('red_pairs d', np.abs(r0[0] - d1).max()); rec('red_pairs gap', np.abs(r0[1] - g1).max())
for C, BP, BQ, ne in ((5.1, 4.8, 4.9, 30), (2.0, 1.9, 2.2, 55), (7.3, 7.31, 7.28, 120)):
    e0, e1 = ns['estats'](C, BP, BQ, ne), en.energy_stats(C, BP, BQ, n_edges=ne)
    rec('estats D2', abs(e0['D2'] - e1['D2'])); rec('estats ratio', abs(e0['ratio'] - e1['ratio'])); rec('estats D2/sqrt(n)', abs(e0['d2n'] - e1['D2_sqrt_n']))
for k, v in worst.items(): print(f'{k:22s} max |diff| = {v:.3g}')
print('ALL_EQUAL' if all(v == 0 for v in worst.values()) else ('MAX', max(worst.values())))
