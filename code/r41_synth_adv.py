"""(round41 g6 추가, 2026-10-21 UTC; Drive JSON round41_2026-10-21 g6_adv — 결과 뒤 추가: 무작위 편집 판에서 대조 (iv)·(v) 도 위반 0 이라 대비가 없음을 본 뒤)
적대적 세 짝 두 생성기, 각 500 (seed 20261088). 함수는 r41_synth.py 를 그대로 불러 쓴다(그 파일은 바꾸지 않음).
  A 줄이고 지우기: x = 기준 열(ℝ: AR(1), 평면: 꺾인 사슬), 항 1–5 개 S 를 골라 y = x 에서 S 를 ε ~ U(0, 0.1) 배로 줄임, z = x 에서 S 를 지움 (평면은 y, z 에 무작위 전역 회전)
  B 늘이기: x = 값 c 를 n 번(n 5–30, 작은 잡음), y = [c, c′], z = c′ 를 n 번 (ℝ 만; 평면은 같은 꼴의 변 벡터)
예측: (i)·(ii)·(iii) 위반 0 (iii 은 인증 구간으로 확실한 위반), (iv) c = 2|a| 는 A 에서 위반 > 0, (v) DTW 는 B 에서 위반 > 0.
사용: python r41_synth_adv.py <ROOT>   출력: results/r41_synth_adv.json"""
import os, sys
ROOT = sys.argv[1]
_s = open(ROOT + '/code/r41_synth.py').read(); exec(_s[:_s.index("\nif __name__ == '__main__':")])
FO2 = RES + '/r41_synth_adv.json'; NA = 500; SEED2 = 20261088


def gen_A_series(rng):
    x = series(rng); k = int(rng.integers(1, 6)); S = np.sort(rng.choice(len(x), min(k, len(x) - 1), replace=False)); y = x.copy(); y[S] *= rng.uniform(0, 0.1, len(S)); z = np.delete(x, S); return x, y, z


def gen_A_curve(rng):
    x = curve(rng); k = int(rng.integers(1, 6)); S = np.sort(rng.choice(len(x), min(k, len(x) - 1), replace=False)); y = x.copy(); y[S] *= rng.uniform(0, 0.1, len(S))[:, None]; z = np.delete(x, S, 0)
    return x, rot(y, rng.uniform(0, 2 * np.pi)), rot(z, rng.uniform(0, 2 * np.pi))


def gen_B_series(rng):
    n = int(rng.integers(5, 31)); c, c2 = rng.normal(0, 1.4, 2); x = c + 0.01 * rng.normal(size=n); y = np.array([c, c2]) + 0.01 * rng.normal(size=2); z = c2 + 0.01 * rng.normal(size=n); return x, y, z


def gen_B_curve(rng):
    n = int(rng.integers(5, 31)); t1, t2 = rng.uniform(0, 2 * np.pi, 2); u1, u2 = np.array([np.cos(t1), np.sin(t1)]), np.array([np.cos(t2), np.sin(t2)])
    x = u1 + 0.01 * rng.normal(size=(n, 2)); y = np.stack([u1, u2]) + 0.01 * rng.normal(size=(2, 2)); z = u2 + 0.01 * rng.normal(size=(n, 2)); return x, rot(y, rng.uniform(0, 2 * np.pi)), rot(z, rng.uniform(0, 2 * np.pi))


def tri_exact(trips, f):
    viol = 0; slack = []
    for x, y, z in trips:
        dxy, dyz, dxz = f(x, y)[0], f(y, z)[0], f(x, z)[0]
        for a, b, c in ((dxz, dxy, dyz), (dxy, dxz, dyz), (dyz, dxy, dxz)):
            slack.append((b + c - a) / max(b + c, 1e-12))
            if a > b + c + 1e-9 * (1 + b + c): viol += 1
    return dict(n_ineq=3 * len(trips), violations=int(viol), slack_rel_min=float(np.min(slack)), slack_rel_q01=float(np.percentile(slack, 1)))


def tri_plane(trips):
    viol = und = 0; slack = []; cap = 0
    for x, y, z in trips:
        Dxy, Dyz, Dxz = plane_D(x, y), plane_D(y, z), plane_D(x, z); cap += sum(D[2]['capped'] for D in (Dxy, Dyz, Dxz))
        for (alo, ahi, _), (blo, bhi, _), (clo, chi, _) in ((Dxz, Dxy, Dyz), (Dxy, Dxz, Dyz), (Dyz, Dxy, Dxz)):
            if alo > bhi + chi: viol += 1
            elif ahi > blo + clo: und += 1
            slack.append((blo + clo - ahi) / max(bhi + chi, 1e-12))
    return dict(n_ineq=3 * len(trips), violations_certain=int(viol), undecided=int(und), slack_rel_min=float(np.min(slack)), slack_rel_q01=float(np.percentile(slack, 1)), bnb_capped=int(cap))


if __name__ == '__main__':
    rng = np.random.default_rng(SEED2); t0 = time.time(); out = dict(n_per_generator=NA, seed=SEED2, generators={})
    TA = [gen_A_series(rng) for _ in range(NA)]; TB = [gen_B_series(rng) for _ in range(NA)]; CA = [gen_A_curve(rng) for _ in range(NA)]; CB = [gen_B_curve(rng) for _ in range(NA)]
    for gname, T_ in (('A', TA), ('B', TB)):
        out['generators'][gname] = {fam: tri_exact(T_, FAM[fam]) for fam in ('i', 'ii', 'iv', 'v')}; log('생성기', gname, out['generators'][gname])
    for gname, T_ in (('A', CA), ('B', CB)):
        out['generators'][gname]['iii'] = tri_plane(T_); log('평면', gname, out['generators'][gname]['iii'], round(time.time() - t0), 's')
    G = out['generators']
    out['prediction'] = dict(i_ii_iii_zero=bool(all(G[g][f]['violations'] == 0 for g in 'AB' for f in ('i', 'ii')) and all(G[g]['iii']['violations_certain'] == 0 for g in 'AB')),
                             iv_on_A_positive=bool(G['A']['iv']['violations'] > 0), v_on_B_positive=bool(G['B']['v']['violations'] > 0))
    out['written'] = time.strftime('%F %T'); out['sec'] = round(time.time() - t0, 1)
    out['rule'] = 'round41 g6_adv (결과 뒤 추가): 예측 (i)·(ii)·(iii) 위반 0, (iv) 는 A 에서 > 0, (v) 는 B 에서 > 0.'
    jdump(out, FO2); log('G6_ADV_DONE', out['prediction'])
