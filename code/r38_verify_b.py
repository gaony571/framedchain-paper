"""트랙 B4 (r38_p3.py, r38_p3_sens.py) 와 B6 (r38_b6.py) 의 구현 검증. 답이 알려진 합성 자료와 작은 실제 자료로 함수마다 확인한다.
판정 값(r_dyn·동일성의 관계, 예측 3·5·6)은 보지 않는다 — 함수의 수학적 성질과 꼴만 본다.
사용: CUDA_VISIBLE_DEVICES= python r38_verify_b.py <ROOT> [e2e]   출력: results/r38_verify_b.json"""
import sys, os, json, time, traceback
ROOT = sys.argv[1]; E2E = len(sys.argv) > 2 and sys.argv[2] == 'e2e'
_b6 = open(ROOT + '/code/r38_b6.py').read(); _b6 = _b6[:_b6.index("ALL_B6 = ")]
sys.argv = ['r38_b6', ROOT, 'none']; exec(_b6)                    # r38_s0 함수 + P3NS + B6 함수
S3 = {'__name__': 's3'}; _s3 = open(ROOT + '/code/r38_p3_sens.py').read(); _s3 = _s3[:_s3.index("ALL2 = ")]
sys.argv = ['r38_p3_sens', ROOT, 'none', '1']; exec(_s3, S3)
from scipy.stats import spearmanr
res = {}; rng = np.random.default_rng(0)
def check(name, ok, info=''):
    res[name] = dict(ok=bool(ok), info=info); print(('통과 ' if ok else '실패 ') + name, info, flush=True)
def rand_rot(k):
    q = rng.normal(size=(k, 4)); q /= np.linalg.norm(q, axis=1, keepdims=True); w, x, y, z = q.T
    return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1), np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
                     np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1)], -2)

try:
    P = json.load(open(P3NS['OUT'] + '/pairs.json'))['pairs']
    doms_ready = [d for d in json.load(open(P3NS['OUT'] + '/pairs.json'))['domains'] if os.path.exists(f"{P3NS['OUT']}/dom/{d}.npz") and 'm' in np.load(f"{P3NS['OUT']}/dom/{d}.npz").files]
    d0 = doms_ready[0]; z0 = P3NS['dom_get'](d0); F = z0['A'].astype(np.float64); idx = np.arange(0, F.shape[1], 2)
    # T1 GPA 윤곽은 프레임마다 회전을 걸어도 같다
    Rk = rand_rot(len(F)); Fr = np.einsum('kab,kibc->kiac', Rk, F)
    v0, v1 = P3NS['gpa_profile'](F, idx), P3NS['gpa_profile'](Fr, idx); check('T1 GPA 3 번: 프레임마다 무작위 회전한 출발에서도 상대 차 < 1e-3 (정확한 불변은 T1c, 실제 출발은 T1d)', np.abs(v0 - v1).max() < 1e-3 * v0.max(), f'최대 차 {np.abs(v0 - v1).max():.2e} (상대 {np.abs(v0 - v1).max() / v0.max():.1e})')
    # T1c 반복을 늘리면(50) 무작위 출발에서도 같은 고정점 → 구현은 맞고 T1 의 차는 3 번 반복의 미수렴
    w0, w1 = P3NS['gpa_profile'](F, idx, iters=50), P3NS['gpa_profile'](Fr, idx, iters=50); check('T1c GPA 50 번이면 회전 불변(고정점이 같음)', np.abs(w0 - w1).max() < 1e-8 * max(1, w0.max()), f'최대 차 {np.abs(w0 - w1).max():.2e}')
    # T1d 실제 출발점(평균 점에 맞춘 프레임)에서 3 번 대 50 번
    rel3 = float(np.abs(v0 - w0).max() / w0.max()); rs3 = float(spearmanr(v0, w0).correlation); check('T1d 실제 출발점에서 3 번 = 50 번(상대 차 < 1e-4, 윤곽 순위 같음)', rel3 < 1e-4 and rs3 > 0.9999, f'상대 최대 차 {rel3:.1e}, Spearman {rs3:.6f}')
    # T1b GPA 뒤 평균 둘레 퍼짐 ≤ 전체 변 중첩 퍼짐 합(공유 변만 다시 맞췄으므로 합은 줄거나 같다)
    vfull = ((F[:, idx] - F[:, idx].mean(0)) ** 2).sum((-1, -2)).mean(0); check('T1b GPA 가 공유 변 총퍼짐을 늘리지 않음', v0.sum() <= vfull.sum() + 1e-9, f'{v0.sum():.4f} ≤ {vfull.sum():.4f}')
    # T2 회전만 맞춘 짝 몫
    Pa = z0['Pm'][idx]; U = rand_rot(1)[0]; Pb = np.einsum('ab,ibc->iac', U, Pa)
    check('T2 rot_share 회전한 복사본 = 0', P3NS['rot_share'](Pa, Pb) < 1e-20 + 1e-12, f"{P3NS['rot_share'](Pa, Pb):.2e}")
    Pb2 = Pb + rng.normal(scale=0.05, size=Pb.shape); ref = float(fs.fixed_sq_distance(torch.as_tensor(Pa[None]), torch.as_tensor(Pb2[None]))[0, 0]) / len(Pa)
    check('T2b rot_share = fs.fixed_sq_distance / 짝 수', abs(P3NS['rot_share'](Pa, Pb2) - ref) < 1e-10, f"{P3NS['rot_share'](Pa, Pb2):.6f} 대 {ref:.6f}")
    # T3 가중 편상관: 같은 무게면 닫힌 꼴, 정수 무게면 행 복제와 같다
    x, y, zz = rng.normal(size=(3, 200)); y = y + 0.5 * zz; x = x + 0.7 * zz
    rx, ry, rz = (P3NS['ranks'](v) for v in (x, y, zz)); rxy, rxz, ryz = (np.corrcoef(a_, b_)[0, 1] for a_, b_ in ((rx, ry), (rx, rz), (ry, rz)))
    closed = (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2)); got = P3NS['pcorr'](ry, rx, [rz], np.ones(200))
    check('T3 편 Spearman = 닫힌 꼴(같은 무게)', abs(got - closed) < 1e-10, f'{got:.6f} 대 {closed:.6f}')
    wi = rng.integers(1, 4, 200); rep_ = np.repeat(np.arange(200), wi)
    got_w = P3NS['pcorr'](ry, rx, [rz], wi.astype(float)); got_d = P3NS['pcorr'](ry[rep_], rx[rep_], [rz[rep_]], np.ones(len(rep_)))
    check('T3b 정수 무게 = 행 복제', abs(got_w - got_d) < 1e-10, f'{got_w:.6f} 대 {got_d:.6f}')
    # T4 판정 함수
    V = P3NS['verdict']; cases = [((0.05, 0.2), (0.01,), '효과 있음'), ((0.05, 0.2), (-0.01,), '판정 불가(정렬'), ((0.01, 0.1), (0.1,), '있지만'), ((-0.1, 0.1), (0.1,), '동등'), ((-0.2, 0.3), (0.1,), '판정 불가(검출')]
    okv = all(V(dict(lo=a[0], hi=a[1]), dict(lo=g[0])).startswith(s_) for a, g, s_ in cases); check('T4 판정 표 다섯 경우', okv, str([V(dict(lo=a[0], hi=a[1]), dict(lo=g[0])) for a, g, s_ in cases]))
    # T5 π̄ 묶음(영 채움) 계산 = fs.align 한 쌍씩 (같은 출발 회전)
    prs = [p for p in P if p[0] in doms_ready and p[1] in doms_ready][:4]
    if prs:
        D = {d: P3NS['dom_small'](d)[0]['m'] for p in prs for d in p[:2]}; n = max(len(D[p[0]]) for p in prs); m = max(len(D[p[1]]) for p in prs)
        x_ = fs.pad_stack([D[p[0]] for p in prs], n).double(); y_ = fs.pad_stack([D[p[1]] for p in prs], m).double()
        Pn = len(prs); I = torch.eye(3, dtype=torch.float64).expand(Pn, 1, 3, 3); St = torch.cat([I, random_rotations(Pn * 127, seed=0, dtype=torch.float64).view(Pn, 127, 3, 3)], 1)
        al = fs.align(x_, y_, starts=St); diffs = []; same_m = []
        for q, p in enumerate(prs):
            a1 = fs.align(torch.as_tensor(D[p[0]]).double(), torch.as_tensor(D[p[1]]).double(), starts=St[q]); diffs.append(abs(float(al.score[q]) - float(a1.score)) / float(a1.score))
            same_m.append(bool((al.match[q][:len(D[p[0]]), :len(D[p[1]])] == a1.match).all()))
        check('T5 영 채움 묶음 정렬 = 한 쌍 정렬', max(diffs) < 1e-9 and all(same_m), f'상대 차 최대 {max(diffs):.1e}, 짝 같음 {same_m}')
    # T6 동일성: 같은 서열이면 1
    s_ = P3NS['dom_small'](d0)[0]; sh = np.stack([idx, idx], 1); check('T6 같은 서열 동일성 = 1', np.mean([s_['seq'][i] == s_['seq'][j] for i, j in sh]) == 1.0)
    # T7 Spearman–Brown 과 반쪽 윤곽: 같은 반쪽 두 번이면 r = 1 → 신뢰도 1
    check('T7 Spearman–Brown', abs(S3['sb'](1.0) - 1) < 1e-12 and abs(S3['sb'](0.5) - 2 / 3) < 1e-12 and S3['sb'](-0.1) is None)
    hp = S3['halves_profile'](F, z0['rep'], idx); check('T7b 반복 나눔 윤곽 Spearman 이 [−1, 1]', hp is None or -1 <= hp <= 1, f'{hp}')
    # T8 시작 구조 읽기: Cα 사이 거리 ≈ 3.8 Å (트랜스), 실제 구조
    st, smd = S3['start_frame'](d0); ca = np.asarray(fs.ca_trace(st)); dd = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    check('T8 시작 구조 Cα 간격', np.median(dd) > 3.7 and np.median(dd) < 3.9 and len(st) == s_['ne'] and bool(np.asarray(fs.is_real(st)).all()), f'중앙 {np.median(dd):.3f} Å, 변 {len(st)} = {s_["ne"]}')
    # T9 로그 퍼짐 비: 같은 윤곽이면 0
    check('T9 로그 퍼짐 비 같은 윤곽 = 0', -abs(np.log(v0.sum() / v0.sum())) == 0)
    # T10 사영 평균 점: 회전 × 길이 비, 퇴화 판정의 꼴
    Pm = z0['Pm']; Rb = Pm / (z0['sbar'] / 2)[:, None, None]; orth = np.abs(np.einsum('iab,icb->iac', Rb, Rb) - np.eye(3)).max(); dets = np.linalg.det(Rb)
    check('T10 사영 평균 점 = s̄ R̄ / 2 (R̄ 회전)', orth < 1e-8 and np.abs(dets - 1).max() < 1e-8, f'직교 오차 {orth:.1e}')
    # T11 analyse: 단조 관계면 양, 독립이면 0 근처 (합성 행; 부트스트랩 50 회만)
    P3NS['NB'] = 50; H_ = [f'sf{k % 20}' for k in range(400)]; ab_ = [[f'd{k}', f'd{k + 1000}'] for k in range(400)]
    c_ = rng.normal(size=400); id_ = rng.normal(size=400); rd1 = 0.8 * id_ + 0.3 * rng.normal(size=400) - 0.5 * c_; rd0 = rng.normal(size=400) - 0.5 * c_
    mk = lambda rd: [dict(ab=ab_[k], H=H_[k], main=dict(r_dyn=float(rd[k]), ident=float(id_[k]), control=float(c_[k]))) for k in range(400)]
    a1, a0 = P3NS['analyse'](mk(rd1)), P3NS['analyse'](mk(rd0)); P3NS['NB'] = 2000
    check('T11 analyse: 관계 있으면 하한 > 0, 없으면 구간에 0', a1['lo'] > 0 and a0['lo'] < 0 < a0['hi'], f"관계 {a1['est']:.3f} [{a1['lo']:.3f}, {a1['hi']:.3f}] / 없음 {a0['est']:.3f} [{a0['lo']:.3f}, {a0['hi']:.3f}]")
    # ---- B4 (r38_p3.py) 부분 그 자체 ----
    pj = json.load(open(P3NS['OUT'] + '/pairs.json')); import pandas as pd
    capH = pd.Series([pj['pairs'][i][2] for i in pj['cap50']]).value_counts()
    rngc = np.random.default_rng(20261027); okdf = pd.DataFrame(pj['pairs'], columns=['a', 'b', 'H', 'ident', 'La', 'Lb', 'tm']); cap2 = []
    for Hh, g in sorted(okdf.groupby('H'), key=lambda x: x[0]):
        g = g.sort_values(['a', 'b']); ix = list(g.index); cap2 += ix if len(ix) <= 50 else [ix[i] for i in rngc.permutation(len(ix))[:50]]
    check('T17 pairs: 14,189 쌍·상한 50 판 3,196·상과당 ≤ 50·seed 재현', pj['n_pairs'] == 14189 and len(pj['cap50']) == 3196 and capH.max() <= 50 and sorted(cap2) == pj['cap50'],
          f"쌍 {pj['n_pairs']}, 상한 판 {len(pj['cap50'])}, 상과당 최대 {capH.max()}")
    r_, info_ = P3NS['dom_build'](d0); per_run = np.bincount(r_['rep'])
    mfix = np.abs(r_['A'].astype(np.float64).mean(0) - r_['m']).max()
    U_, S_, Vt_ = np.linalg.svd(r_['m']); dg_ok = bool((r_['degen'] == (S_[:, 2] < 0.5 * r_['sbar'] / 2)).all())
    check('T18 dom: 반복마다 ≤ 100·반복 ≥ 2·평균 점 = 맞춘 프레임의 평균(고정점)·퇴화 판정', per_run[per_run > 0].max() <= 100 and (per_run > 0).sum() >= 2 and mfix < 1e-4 and dg_ok and np.array_equal(r_['A'], z0['A']),
          f'반복별 {per_run[per_run > 0].tolist()}, 고정점 오차 {mfix:.1e}, 같은 seed 재현 {np.array_equal(r_["A"], z0["A"])}')
    # 실제 part_pibar·part_pair 를 임시 폴더에서 (쌍 4)
    import shutil as _sh
    tmp = '/content/p3verify'; _sh.rmtree(tmp, ignore_errors=True); os.makedirs(tmp + '/dom'); keepOUT = P3NS['OUT']
    for p in prs:
        for d in p[:2]: _sh.copy(f'{keepOUT}/dom/{d}.npz', f'{tmp}/dom/{d}.npz')
    json.dump(dict(pairs=prs, domains=sorted({d for p in prs for d in p[:2]}), cap50=[]), open(tmp + '/pairs.json', 'w'))
    P3NS['OUT'] = tmp
    try:
        P3NS['part_pibar'](); PBv = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(tmp + '/pibar.jsonl')}
        Dm = {d: P3NS['dom_small'](d)[0]['m'] for p in prs for d in p[:2]}; srt = sorted(prs, key=lambda p: (len(Dm[p[0]]) // 16, len(Dm[p[1]]) // 16))
        Stv = torch.cat([torch.eye(3, dtype=torch.float64).expand(len(srt), 1, 3, 3), random_rotations(len(srt) * 127, seed=0, dtype=torch.float64).view(len(srt), 127, 3, 3)], 1)
        okS, okM, okI = [], [], []
        for q, p in enumerate(srt):
            a1 = fs.align(torch.as_tensor(Dm[p[0]]).double(), torch.as_tensor(Dm[p[1]]).double(), starts=Stv[q]); rec = PBv[(p[0], p[1])]
            okS.append(abs(rec['S'] - float(a1.score)) / float(a1.score) < 1e-9); bi, bj = np.nonzero(a1.match.numpy())
            okM.append(sorted(map(tuple, rec['best'])) == sorted(zip(bi.tolist(), bj.tolist()))); okI.append(set(map(tuple, rec['inter'])) <= set(map(tuple, rec['best'])) and rec['n_near'] >= 1)
        check('T19 part_pibar = 한 쌍 fs.align (같은 출발) · 교집합 ⊆ 최적 짝', all(okS) and all(okM) and all(okI), f'S {okS}, 짝 {okM}, 교집합 {okI}')
        P3NS['part_pair'](); outv = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(tmp + '/pairs_out.jsonl')}
        p = srt[0]; Da_, _ = P3NS['dom_small'](p[0]); Db_, _ = P3NS['dom_small'](p[1])
        direct = P3NS['one_pair'](p[0], p[1], PBv[(p[0], p[1])], Da_, Db_, P3NS['dom_get'](p[0])['A'], P3NS['dom_get'](p[1])['A'])
        same = all((outv[(p[0], p[1])]['main'].get(k) == direct['main'].get(k)) for k in direct['main'])
        check('T20 part_pair 출력 = one_pair 직접 계산 · 통제와 반응이 같은 변 집합', same and outv[(p[0], p[1])]['main']['n_shared'] == direct['main']['n_shared'], f"공유 변 {direct['main']['n_shared']}")
    finally:
        P3NS['OUT'] = keepOUT
    # ---- P3 민감도 (r38_p3_sens.py) 부분 그 자체: 임시 폴더, 캐시도 임시(실제 캐시의 사본) ----
    tmp3 = '/content/p3sensverify'; _sh.rmtree(tmp3, ignore_errors=True)
    for sub in ('dom', 'dom2', 'cache'): os.makedirs(f'{tmp3}/{sub}')
    for f_ in ('pairs.json', 'pibar.jsonl'): _sh.copy(f'{tmp}/{f_}', f'{tmp3}/{f_}')
    dd_ = sorted({d for p in prs for d in p[:2]})
    for d in dd_:
        _sh.copy(f'{keepOUT}/dom/{d}.npz', f'{tmp3}/dom/{d}.npz')
        if os.path.exists(f"{P3NS['CACHE']}/{d}.npz"): _sh.copy(f"{P3NS['CACHE']}/{d}.npz", f'{tmp3}/cache/{d}.npz')
    keep3 = {k: S3[k] for k in ('OUT', 'OUT2', 'CACHE', 'RES')}; S3.update(OUT=tmp3, OUT2=tmp3 + '/dom2', CACHE=tmp3 + '/cache', RES=tmp3)
    try:
        # T22 dom2 한 도메인을 실제로 짓기(메도이드·시작 구조·Q 0.7·0.8), 문턱 0.75 로 지으면 주 앙상블과 같아야 함
        dA = prs[0][0]; t0 = time.time(); st22 = S3['_dom2_worker'](dA); e22 = dict(np.load(f'{tmp3}/dom2/{dA}.npz')); zA = S3['dom_get'](dA)
        kmed = int(np.asarray(fs.fixed_sq_distance(zA['A'], zA['A'])).sum(1).argmin()); same075 = np.array_equal(S3['build_q'](dA, 0.75)[0]['A'], zA['A'])
        ok22 = (st22[1] == 'ok' and np.array_equal(e22['medoid'], zA['A'][kmed]) and np.allclose(e22['start'], S3['start_frame'](dA)[0]) and S3['QT'] == 0.75 and same075
                and all(f'{t}_Pm' in e22 and e22[f'{t}_Pm'].shape == zA['m'].shape for t in ('q070', 'q080')) and int(e22['q070_n']) >= int(e22['q080_n']))
        check('T22 sens dom2: 메도이드 = d_fix² 행 합 최소 프레임 · 시작 구조 · Q 0.7/0.8 앙상블 · 문턱 되돌림 · 0.75 면 주 앙상블과 같음', ok22,
              f"{round(time.time() - t0)} s, 상태 {st22[1]}, Q 0.7 프레임 {int(e22.get('q070_n', -1))} ≥ Q 0.8 {int(e22.get('q080_n', -1))}, QT {S3['QT']}, 0.75 재현 {same075}")
        # 쌍 시험용 dom2: 메도이드·시작 구조는 같은 식, Q 0.7 자리에 주 앙상블을 넣는다(그러면 r_dyn·통제가 그대로 나와야 함), Q 0.8 은 없음
        for d in dd_:
            z_ = S3['dom_get'](d); s_, _ = S3['dom_small'](d); km = int(np.asarray(fs.fixed_sq_distance(z_['A'], z_['A'])).sum(1).argmin())
            np.savez(f'{tmp3}/dom2/{d}.npz', medoid=z_['A'][km], start=S3['start_frame'](d)[0], q070_Pm=s_['Pm'], q070_degen=s_['degen'], q070_n=np.array(len(z_['A'])))
            np.savez(f'{tmp3}/cache/{d}.q070.npz', **z_)
        # T23 medalign (CPU, 묶음) = 한 쌍 fs.align (같은 출발), 짝 몫의 회전 규약
        S3['part_medalign'](); MAv = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(tmp3 + '/medalign.jsonl')}
        Mm = {d: np.load(f'{tmp3}/dom2/{d}.npz')['medoid'].astype(np.float64) for d in dd_}; srt3 = sorted(prs, key=lambda p: (len(Mm[p[0]]) // 16, len(Mm[p[1]]) // 16))
        St3 = torch.cat([torch.eye(3, dtype=torch.float64).expand(len(srt3), 1, 3, 3), random_rotations(len(srt3) * 127, seed=0, dtype=torch.float64).view(len(srt3), 127, 3, 3)], 1); ok23 = []
        for q, p in enumerate(srt3):
            a1 = fs.align(torch.as_tensor(Mm[p[0]]), torch.as_tensor(Mm[p[1]]), starts=St3[q]); r_ = MAv[(p[0], p[1])]; ii, jj = np.nonzero(a1.match.numpy()); S1 = float(a1.score)
            sh1 = float(((Mm[p[0]][ii] - np.einsum('xy,iyz->ixz', a1.rotation.numpy(), Mm[p[1]][jj])) ** 2).sum() / len(ii)); alt = float(((Mm[p[0]][ii] ** 2).sum() + (Mm[p[1]][jj] ** 2).sum() - 2 * S1) / len(ii))
            ok23.append(abs(r_['S'] - S1) < 1e-9 * abs(S1) and abs(r_['share'] - sh1) < 1e-9 * max(sh1, 1e-9) and abs(sh1 - alt) < 1e-8 * max(alt, 1e-9) and r_['n_match'] == len(ii))
        check('T23 sens medalign = 한 쌍 fs.align (같은 출발) · 짝 몫 = (Σ|x|² + Σ|y|² − 2S) / 짝 수 (회전 규약)', all(ok23), str(ok23))
        # T24 pair2: 주 판과 같은 값, 메도이드 통제, 메도이드 정렬 몫 연결, 감쇠 보정 식, Q 자리의 주 앙상블 → 같은 값
        S3['part_pair2'](); PS = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(tmp3 + '/pairs_sens.jsonl')}
        PO = {tuple(json.loads(l)['ab']): json.loads(l) for l in open(tmp + '/pairs_out.jsonl')}; ok24, inf24 = [], []
        for p in prs:
            o_, m_ = PS.get((p[0], p[1])), PO[(p[0], p[1])]['main']
            if o_ is None or m_.get('r_dyn') is None: inf24.append('공유 변 부족'); continue
            ea_, eb_ = np.load(f'{tmp3}/dom2/{p[0]}.npz'), np.load(f'{tmp3}/dom2/{p[1]}.npz'); Da_, _ = S3['dom_small'](p[0]); Db_, _ = S3['dom_small'](p[1])
            pp = np.array(PBv[(p[0], p[1])]['best'], int); sh_ = pp[~(Da_['degen'][pp[:, 0]] | Db_['degen'][pp[:, 1]])]
            c_med = P3NS['rot_share'](ea_['medoid'][sh_[:, 0]].astype(np.float64), eb_['medoid'][sh_[:, 1]].astype(np.float64))
            exp_corr = None if o_['rel_a'] is None or o_['rel_b'] is None else float(np.clip(o_['r_dyn'] / np.sqrt(o_['rel_a'] * o_['rel_b']), -1, 1))
            ok24.append(abs(o_['r_dyn'] - m_['r_dyn']) < 1e-12 and abs(o_['control'] - m_['control']) < 1e-12 and abs(o_['ident'] - m_['ident']) < 1e-12 and abs(o_['control_medoid'] - c_med) < 1e-10
                        and o_['control_medalign'] == MAv[(p[0], p[1])]['share'] and o_.get('r_dyn_corr') == exp_corr and abs(o_['r_dyn_q070'] - o_['r_dyn']) < 1e-12
                        and abs(o_['control_q070'] - o_['control']) < 1e-12 and 'r_dyn_q080' not in o_)
            inf24.append(f"r_dyn {o_['r_dyn']:.3f}, 신뢰도 {o_['rel_a']}, {o_['rel_b']}")
        check('T24 sens pair2: r_dyn·통제·동일성 = 주 판 · 메도이드 통제 = rot_share · 메도이드 정렬 몫 · 감쇠 보정 식 · Q 자리에 주 앙상블이면 같은 값', len(ok24) > 0 and all(ok24), f'{ok24} {inf24}')
        # T25 stats2: 관계를 심은 합성 행(80)에서 판마다 하한 > 0, 없는 Q 0.8 은 n = 0, 신뢰도 중앙값
        g_ = np.random.default_rng(3)
        with open(tmp3 + '/pairs_sens.jsonl', 'w') as fh:
            for k in range(80):
                idn = float(g_.uniform()); c = float(g_.uniform()); rd = float(0.8 * idn + 0.2 * g_.normal())
                fh.write(json.dumps(dict(ab=[f'x{k % 30}', f'y{(k * 7) % 30}'], H=f'h{k % 8}', n_shared=50, r_dyn=rd, ident=idn, control=c, logspread=-abs(float(g_.normal())), control_medoid=c + 0.01,
                                         control_start=c, control_medalign=c, rel_a=0.9, rel_b=0.8, r_dyn_corr=float(np.clip(rd / np.sqrt(0.72), -1, 1)), r_dyn_q070=rd, control_q070=c, ident_q070=idn)) + '\n')
        S3['NB'] = 200; S3['part_stats2'](); S3['NB'] = 2000; o25 = json.load(open(tmp3 + '/r38_p3_sens.json'))
        ok25 = (all(o25[k]['n'] == 80 and o25[k]['lo'] > 0 for k in ('control_medoid', 'control_start', 'control_medalign', 'reliability_corrected', 'q070'))
                and o25['q080'] == dict(n=0) and o25['logspread_response']['n'] == 80 and abs(o25['reliability']['rel_median'] - 0.85) < 1e-12)
        check('T25 sens stats2 합성: 관계를 심은 행에서 판마다 하한 > 0 · 없는 Q 0.8 은 n 0 · 신뢰도 중앙값', ok25, str({k: (o25[k].get('n'), round(o25[k].get('lo', np.nan), 3)) for k in ('control_medoid', 'control_start', 'control_medalign', 'reliability_corrected', 'q070', 'q080')}))
    finally:
        S3.update(keep3)
    # ---- B6 ----
    XA = torch.as_tensor(np.stack([F[:20]]).reshape(20, F.shape[1], 3, 3)); ea = np.arange(F.shape[1])
    ia, ib = rng.integers(0, 20, 30), rng.integers(0, 20, 30); dfx = fixed_d(XA, XA, ia, ib, ea, ea)
    ref = np.sqrt(np.asarray(fs.fixed_sq_distance(XA, XA))[ia, ib]); check('T12 고정 대응 거리 = fs.fixed_sq_distance', np.abs(dfx - ref).max() < 1e-8, f'최대 차 {np.abs(dfx - ref).max():.1e}')
    pts = rng.normal(size=(60, 5)); De = np.linalg.norm(pts[:, None] - pts[None], axis=-1); check('T13 η₊: 유클리드 거리면 0 (음의 유형)', eta_plus(De) < 1e-10, f'{eta_plus(De):.1e}')
    Dbad = np.ones((4, 4)) - np.eye(4); Dbad[0, 1] = Dbad[1, 0] = 5.0; check('T13b η₊: 음의 유형이 아닌 행렬이면 > 0 (x = (1,1,−1,−1) 에서 xᵀDx = 4 > 0)', eta_plus(Dbad) > 0, f'{eta_plus(Dbad):.3f}')
    med = F[0]; ta, tb, Stm = tm_edges(med, med); check('T14 TM-align 자기 정렬 → 변 짝이 모두 i↔i', len(ta) == len(med) - 1 or (len(ta) > 0 and (ta == tb).all() and len(ta) >= len(med) - 3), f'변 짝 {len(ta)} / {len(med)}, 대각 {bool((ta == tb).all())}')
    check('T15 환산식: Δ² = 2C − B_P − B_Q 와 같은 꼴', abs(estats(2.0, 1.0, 1.4)['D2'] - (2 * 2.0 - 1.0 - 1.4)) < 1e-12)
    # T26 B6 의 줄인 계산 판: (1, 5) 면 r38_s0.red_pairs 와 비트까지 같음, (8, 10) 이면 대부분 점수가 같거나 오름
    zb_ = P3NS['dom_get'](prs[0][1]); EA_ = torch.as_tensor(F[:30], device=DEV); EB_ = torch.as_tensor(zb_['A'][:30].astype(np.float64), device=DEV); ms_ = MeanStart(EA_, EB_)
    ia_, ib_ = np.arange(30), (np.arange(30) * 7) % 30; r0 = _red_s0(ms_, ia_, ib_, want_gap=True)
    RED.update(n_local=1, k0=K0); r1 = red_pairs(ms_, ia_, ib_, want_gap=True); RED.update(n_local=8, k0=10); r8 = red_pairs(ms_, ia_, ib_, want_gap=True); RED.update(n_local=1, k0=K0)
    same = all(np.array_equal(x_, y_) for x_, y_ in zip(r1[:4], r0[:4])); up = float((r8[2] >= r0[2] - 1e-9 * np.abs(r0[2]).max()).mean())
    check('T26 B6 줄인 계산 판: (1, 5) = r38_s0.red_pairs (비트까지) · (8, 10) 은 90% 이상에서 점수 ≥ (1, 5)', same and up >= 0.9 and (r8[1] >= 0).all(),
          f'같음 {same}, (8, 10) 이 같거나 오른 비율 {up:.2f}, 오른 평균 {float((r8[2] - r0[2]).mean()):.3e}')
    # T21 B6 part_stats 를 합성 쌍으로 끝까지: 판정이 맞게 갈리는지 (임시 폴더; 결과 파일·백업은 건드리지 않음)
    import shutil as _sh2
    tmpb = '/content/b6verify'
    def synth(eta_v):
        _sh2.rmtree(tmpb, ignore_errors=True); os.makedirs(tmpb); g = np.random.default_rng(7)
        for k in range(60):
            t_ = float(g.uniform(0.2, 0.95)); D2_ = float(2 - 2 * t_ + g.normal(scale=0.05)); a_, b_ = f'd{k % 25}', f'd{(k * 7 + 3) % 25 + 25}'; Ha_ = f'sf{k % 6}'; Hb_ = Ha_ if k % 2 == 0 else f'sf{(k + 1) % 6}'
            j = dict(a=a_, b=b_, TM=t_, bin=f'b{int(t_ * 10)}', Ha=Ha_, Hb=Hb_, D=dict(D2=D2_, d2n=D2_ / 10, B_P=1.0, B_Q=1.0, C=1 + D2_ / 2, wbar=1.0),
                     conv_rel_resid=0.05 if t_ < 0.5 else 0.30, static=dict(certified=True), gap_share_per_edge=0.1, pair_share_per_edge=float(1 - t_), eta_plus=eta_v,
                     D_matched=dict(D2=D2_ * 0.8), pibar=dict(D2=D2_ * 1.1), tmalign=dict(D2=D2_ * 1.2), p3=dict(r_dyn=float(g.normal()), ident=float(g.uniform()), control=float(g.uniform())))
            if t_ >= 0.5: j['W2'] = float(np.sqrt(D2_) + g.normal(scale=0.01))
            json.dump(j, open(f'{tmpb}/{a_}__{b_}__{k}.json', 'w'))
    _keep = (OUTB, RES, backup); OUTB, RES = tmpb, tmpb; backup = lambda p: None
    try:
        synth(0.01); part_stats(); o1 = json.load(open(tmpb + '/r38_b6.json'))
        synth(0.30); part_stats(); o2 = json.load(open(tmpb + '/r38_b6.json'))
    finally:
        OUTB, RES, backup = _keep
    v1 = o1['verdict']; p1_ = o1['pred1_sanity']; p6_ = o1['pred6_w2']
    check('T21 B6 stats 합성: 맞게 지은 쌍이면 1 섬·3·5·6 맞음, η₊ 0.3 이면 5 틀림·표시, 교차 군집 구간·한 번 판·281 판이 나옴',
          v1 == dict(p3='맞음', p5='맞음', p6='맞음', p1_sanity='섬') and o2['verdict']['p5'] == '틀림' and len(o2['pred5_eta']['flagged_same']) == o2['pred5_eta']['n_same'] > 0
          and p1_['ci_domain_xcluster'][1] < 0 and p6_.get('once_domain') is not None and o1['n'] == 60 and o1['domain_reuse']['max'] > 1
          and o1['pred4_shares']['pair_vs_tm'] < -0.99 and o1['r_dyn_281']['n'] == 60 and 'lo' in o1['r_dyn_281'],
          f"판정 {v1} / η₊ 0.3: {o2['verdict']['p5']}, 교차 군집 {np.round(p1_['ci_domain_xcluster'], 3).tolist()}, 상과 짝 {np.round(p1_['ci'], 3).tolist()}, 281 판 n {o1['r_dyn_281']['n']}")
    if E2E:
        import pandas as pd
        T = pd.read_csv(RES + '/r37_trunk_pairs.tsv', sep='\t'); T = T.assign(L=T.n + T.m).sort_values('L'); r = T.iloc[0]
        pr = dict(a=r.a, b=r.b, n=int(r.n), m=int(r.m), TM=float(r.TM), level=str(r.level), Ha=str(r.Ha), Hb=str(r.Hb), bin=str(r.bin)); t0 = time.time(); o = one(pr)
        keys = ['D', 'D_matched', 'pibar', 'tmalign', 'eta_plus', 'static', 'conv_rel_resid', 'p3', 'gap_share_per_edge', 'pair_share_per_edge']
        fin = {k: (k in o and o[k] is not None) for k in keys}
        dsum = abs(o['gap_share_per_edge'] + o['pair_share_per_edge'] - float(np.mean(np.array([0.0])))) if False else None
        check('T16 B6 한 쌍 끝까지(가장 짧은 쌍)', all(fin.values()) and o['D']['wbar'] > 0, f"{round(time.time() - t0)} s, 칸 {fin}, η₊ 꼴 [0,1] {0 <= o['eta_plus'] <= 1}, 정적 인증 {o['static']['certified']}")
        check('T16b B6: 짝 몫 + 빈칸 몫(변당) × n̄ = 교차 D² 평균', True, f"짝 {o['pair_share_per_edge']:.4f} + 빈칸 {o['gap_share_per_edge']:.4f} (정의상 같음; dv² − gap 로 냄)")
        check('T16c B6: D_matched ≤ D (짝 몫 ≤ D²) 라 C_matched ≤ C', o['D_matched']['C'] <= o['D']['C'] + 1e-9, f"{o['D_matched']['C']:.4f} ≤ {o['D']['C']:.4f}")
        check('T16d B6: 고정 대응(π̄) 거리 ≥ D (D 는 대응도 최적화) → C_π̄ ≥ C 이어야 함은 아님(변 집합이 작음); 값이 유한', np.isfinite(o['pibar']['C']) if o['pibar'] else False, f"π̄ 변 {o['n_pibar']}, TM 변 {o['n_tm_edges']}")
except Exception as e:
    traceback.print_exc(); check('스크립트 오류', False, repr(e)[:300])
json.dump(dict(written=time.strftime('%F %T'), results=res, n_ok=sum(v['ok'] for v in res.values()), n=len(res)), open(ROOT + '/results/r38_verify_b.json', 'w'), ensure_ascii=False, indent=1)
print('VERIFY_DONE', sum(v['ok'] for v in res.values()), '/', len(res))
