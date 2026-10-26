"""트랙 B5 (r38_p1.py) 구현 검증. 답이 알려진 합성 자료와 실제 한 쌍. 판정 값(θ* 의 상과 중앙값, 조건 1·2)은 보지 않는다.
사용: CUDA_VISIBLE_DEVICES= python r38_verify_p1.py <ROOT> [e2e]   출력: results/r38_verify_p1.json"""
import sys, os, json, time, traceback
ROOT = sys.argv[1]; E2E = len(sys.argv) > 2 and sys.argv[2] == 'e2e'
P1 = {'__name__': 'p1v'}; _s = open(ROOT + '/code/r38_p1.py').read(); _s = _s[:_s.index("ALL1 = [")]
sys.argv = ['r38_p1', ROOT, 'none', '1']; exec(_s, P1)
np, torch, fs = P1['np'], P1['torch'], P1['fs']
res = {}
def check(name, ok, info=''):
    res[name] = dict(ok=bool(ok), info=info); print(('통과 ' if ok else '실패 ') + name, info, flush=True)

try:
    th = P1['theta']; T4 = [320., 348., 379., 413.]
    a = th([-0.2, -0.1, 0.1, 0.3], T4); b = th([0.05, 0.1, 0.2, 0.3], T4); c = th([-0.3, -0.2, -0.1, -0.05], T4); d = th([-0.2, 0.1, -0.05, 0.2], T4)
    e = th([-0.3, 0.1, 0.2], [320., 379., 413.]); f = th([0.0, 0.1], [320., 348.])
    check('V1 θ* 보간·검열·비단조·빠진 온도', a[:2] == ('num', 363.5) and b[:2] == ('lt', 320.0) and c[:2] == ('gt', 413.0) and abs(d[1] - (320 + 0.2 * 28 / 0.3)) < 1e-9 and d[2]
          and abs(e[1] - (320 + 0.3 * 59 / 0.4)) < 1e-9 and f[:2] == ('num', 320.0) and not a[2], str([a, b, c, d, e, f]))
    ss, wm, kh, km = P1['sf_summary'], P1['wmedian'], P1['key_ht'], P1['key_mm']
    s1 = ss([('num', 350.), ('gt', 348.), ('num', 400.)], 'ht'); s2 = ss([('gt', 379.), ('num', 379.), ('lt', 320.)], 'ht')
    s3 = ss([('gt', 348.), ('num', 360.), ('num', 400.)], 'mm'); s4 = ss([('gt', 348.), ('gt', 379.), ('num', 360.)], 'mm')
    s5 = wm([('num', x) for x in (1., 2., 3., 4.)], kh); s6 = wm([('num', x) for x in (1., 2., 3., 4.)], kh, [1, 1, 1, 10])
    check('V2 상과 요약: 호열–중온 하한 순위·대칭 판 가장 큰 값·u_C·아래 중앙값·무게', s1 == ('num', 350.) and s2 == ('num', 379.) and s3 == ('num', 400.) and s4 == ('gt', 360.)
          and s5 == ('num', 2.) and s6 == ('num', 4.), str([s1, s2, s3, s4, s5, s6]))
    ds = P1['dsum']; cases = [((('num', 360.), ('num', 340.)), ('num', 20.)), ((('gt', 413.), ('num', 340.)), ('num', 73.)), ((('num', 360.), ('lt', 320.)), ('num', 40.)),
                              ((('gt', 379.), ('lt', 320.)), ('num', 59.)), ((('lt', 320.), ('num', 340.)), ('le0', 0.)), ((('num', 350.), ('gt', 360.)), ('le0', 0.)),
                              ((('num', 370.), ('gt', 360.)), ('excl', '부호 미정')), ((('lt', 320.), ('gt', 360.)), ('le0', 0.)), ((('gt', 379.), ('gt', 348.)), ('excl', '같은 방향 이중 검열')),
                              ((('lt', 320.), ('lt', 320.)), ('excl', '같은 방향 이중 검열'))]
    okd = [ds(*x) == y for x, y in cases]; check('V3 조건 2 의 d_s 표 열 경우', all(okd), str(okd))
    rng = np.random.default_rng(0); c1 = P1['cond1']
    Ta = {f'h{i}': ('lt', 320.) if i < 12 else ('num', 400.) for i in range(20)}; Tb = {f'h{i}': ('gt', 379.) if i < 12 else ('num', 330.) for i in range(20)}
    Tc = {f'h{i}': ('num', 340. + i) for i in range(30)}; Td = {f'h{i}': ('num', 300. + i) for i in range(30)}
    r_a, r_b, r_c, r_d = c1(Ta, rng, 300), c1(Tb, rng, 300), c1(Tc, rng, 300), c1(Td, rng, 300)
    check('V4 조건 1: 왼쪽 검열 > 50% 틀림·오른쪽 > 50% 맞음·하한 > 320', (not r_a['passed']) and r_b['passed'] and r_c['passed'] and not r_d['passed'],
          f"{r_a['passed']}, {r_b['passed']}, {r_c['passed']} lo {r_c['ci95'][0]}, {r_d['passed']} lo {r_d['ci95'][0]}")
    c2 = P1['cond2']; T_s = {f'h{i}': ('num', 360. + i) for i in range(25)}; C_s = {f'h{i}': ('num', 330. + i) for i in range(25)}; C0 = {f'h{i}': ('num', 360. + i) for i in range(25)}
    q1, q0 = c2(T_s, C_s, rng, 300), c2(T_s, C0, rng, 300); q2 = c2({k: v for k, v in list(T_s.items())[:10]}, {k: v for k, v in list(C_s.items())[:10]}, rng, 300)
    check('V5 조건 2: 짝지은 판(차 30 → 섬, 차 0 → 안 섬), 상과 < 20 이면 짝짓지 않은 판', q1['passed'] and q1['version'] == '짝지은 판' and q1['d'] == ['num', 30.0] and not q0['passed']
          and q2['version'] == '짝짓지 않은 판', f"{q1['d']} {q1['ci95']}, {q0['d']}, {q2['version']}")
    nz = [dict(H=f'h{i % 20}', temps=[320., 348., 379., 413.], se_h=[0.] * 4, se_r=0.) for i in range(60)]
    m1, lo1, hi1 = P1['sim_once'](nz, 350., 0., rng, 300); m2, lo2, hi2 = P1['sim_once'](nz, 420., 0., rng, 300); m3, _, _ = P1['sim_once'](nz, 320., 0., rng, 300)
    check('V6 시뮬레이션 잡음 0: θ 그대로(350), 413 너머는 "> 413", 320 은 320', abs(m1[1] - 350) < 1e-9 and lo1 == hi1 == m1 and m2 == ('gt', 413.) and m3 == ('num', 320.), f'{m1} {m2} {m3}')
    X = rng.normal(size=(200, 4)); y = (X[:, 0] + 0.5 * rng.normal(size=200) > 0).astype(float); ps = P1['logistic_w'](X, y)
    check('V7 성향 점수 로지스틱: 0–1, 1 쪽이 더 큼', ps.min() > 0 and ps.max() < 1 and ps[y == 1].mean() > ps[y == 0].mean() + 0.3, f'{ps[y == 1].mean():.2f} 대 {ps[y == 0].mean():.2f}')
    A = rng.normal(size=(50, 40, 3, 3)); idx = np.arange(0, 40, 2); fsp = P1['full_spread'](A, idx); man = ((A[:, idx] - A[:, idx].mean(0)) ** 2).sum((-1, -2)).mean(0)
    check('V8 전체 변 중첩 퍼짐 = 손 계산', np.abs(fsp - man).max() < 1e-12)
    # pair_theta: 가짜 퍼짐 사전으로 규칙 확인 (320 기준 없음 / 데워진 쪽 320 없음 / 온도 점 < 2 / 413 뺌)
    ck = P1['ckey']; mk = lambda vals: {ck(T, 'half', 0.75): dict(gpa=v, full=v * 1.1, runs=[v, v], alts=[v]) for T, v in vals.items()}
    sp = {'H1': dict(pairs={'H1__R1': dict(n_shared=40, cells=mk({'320': 1.0, '348': 1.1, '379': 1.3, '413': 1.6}))}), 'R1': dict(pairs={'H1__R1': dict(n_shared=40, cells=mk({'320': 1.2}))})}
    r_ = dict(a='H1', b='R1', hot='H1', ref='R1', H='h'); t0 = P1['pair_theta'](sp, r_, ('half', 0.75)); t413 = P1['pair_theta'](sp, r_, ('half', 0.75), drop413=True)
    sp2 = json.loads(json.dumps(sp)); del sp2['R1']['pairs']['H1__R1']['cells'][ck('320', 'half', 0.75)]; t_noref = P1['pair_theta'](sp2, r_, ('half', 0.75))
    sp3 = json.loads(json.dumps(sp)); sp3['H1']['pairs']['H1__R1']['cells'] = mk({'348': 1.1, '379': 1.3}); t_nohot = P1['pair_theta'](sp3, r_, ('half', 0.75))
    exp = 348 + (np.log(1.2) - np.log(1.1)) * 31 / (np.log(1.3) - np.log(1.1))
    check('V9 pair_theta: 보간 값·413 뺌·기준 320 없음·데워진 쪽 320 없음', t0[0] == 'num' and abs(t0[1] - exp) < 1e-9 and t413[0] == 'num' and abs(t413[1] - exp) < 1e-9
          and t_noref == ('skip', '기준 쪽 320 K 칸 없음') and t_nohot == ('skip', '데워진 쪽 320 K 칸 없음'), f'{t0} {t_noref} {t_nohot}')
    if E2E:
        pj = json.load(open(P1['OUT1'] + '/pairs.json')); rr = sorted(pj['HT'], key=lambda r: r['La'] + r['Lb'])[0]; t0 = time.time()
        ca, la, _ = P1['build_cells'](rr['a']); cb, lb, _ = P1['build_cells'](rr['b']); ca2, _, _ = P1['build_cells'](rr['a'], temps=('320',), cells=[P1['MAIN']])
        k = ck('320', 'half', 0.75); det = k in ca and np.array_equal(ca[k][0], ca2[k][0])
        nfr = {kk: (len(v[0]), len(set(v[1].tolist()))) for kk, v in ca.items()}; okn = all(n == r * P1['NWIN'][kk.split('|')[1]] for kk, (n, r) in nfr.items())
        check('V10 실제 칸: 같은 seed 재현·프레임 수 = N × 반복·창별 N', det and okn and len(ca) > 0, f"{rr['a']} {sorted(nfr.items())[:6]} 길이 {la}")
        if k in ca and k in cb:
            Ms = {rr['a']: ca[k][2], rr['b']: cb[k][2]}; S, best, alts, nn = P1['align_batch']([(rr['a'], rr['b'])], Ms)[0]
            St = torch.cat([torch.eye(3, dtype=torch.float64)[None], P1['random_rotations'](127, seed=P1['zlib'].crc32(f"r38_p1|{rr['a']}|{rr['b']}".encode()) % (2 ** 31), dtype=torch.float64)], 0)
            a1 = fs.align(torch.as_tensor(Ms[rr['a']]), torch.as_tensor(Ms[rr['b']]), starts=St)
            same = abs(S - float(a1.score)) < 1e-9 * abs(S) and sorted(map(tuple, best)) == sorted(map(tuple, np.argwhere(a1.match.numpy()).tolist()))
            check('V11 align_batch = fs.align (같은 출발)·거의 최적 대응 수', same and nn >= 1, f'S {S:.4f}, 공유 변 {len(best)}, 다른 대응 {len(alts)}, 거의 최적 {nn}')
            P1['AL'].clear(); P1['AL'][(rr['a'], rr['b'])] = dict(best=best, alts=alts)
            fo = f"{P1['OUT1']}/spread/{rr['a']}.json"; existed = os.path.exists(fo)
            if not existed:
                tmpo = P1['OUT1']; P1['OUT1'] = '/content/p1verify'; os.makedirs('/content/p1verify/spread', exist_ok=True)
                try:
                    s_a = P1['_spread_worker'](rr['a']); s_b = P1['_spread_worker'](rr['b'])
                    spv = {d: json.load(open(f'/content/p1verify/spread/{d}.json')) for d in (rr['a'], rr['b'])}
                finally:
                    P1['OUT1'] = tmpo
                tt = P1['pair_theta'](spv, rr, P1['MAIN']); key = f"{rr['a']}__{rr['b']}"; c0 = spv[rr['hot']]['pairs'][key]['cells'].get(k, {})
                sh = np.array(sorted({x[0 if rr['hot'] == rr['a'] else 1] for x in best}), int); Ahot = (ca if rr['hot'] == rr['a'] else cb)[k][0]
                direct = float(P1['gpa_profile'](Ahot, sh).sum())
                check('V12 실제 한 쌍 끝까지: 퍼짐 = gpa_profile 직접, θ* 꼴', s_a[1] == 'ok' and s_b[1] == 'ok' and abs(c0.get('gpa', -1) - direct) < 1e-9 * direct and tt[0] in ('num', 'gt', 'lt'),
                      f"{round(time.time() - t0)} s, θ* 꼴 {tt[0]}, 칸 {sorted(spv[rr['hot']]['pairs'][key]['cells'])[:4]}…, 반복 퍼짐 {len(c0.get('runs', []))}")
except Exception as e:
    traceback.print_exc(); check('스크립트 오류', False, repr(e)[:300])
json.dump(dict(written=time.strftime('%F %T'), results=res, n_ok=sum(v['ok'] for v in res.values()), n=len(res)), open(ROOT + '/results/r38_verify_p1.json', 'w'), ensure_ascii=False, indent=1)
print('VERIFY_P1_DONE', sum(v['ok'] for v in res.values()), '/', len(res))
