"""(prereg round36 0.2 의 다음 단계, 2026-10-20) 첫 세션 상계 편향에서 TM 0.2–0.3 구간이 문턱을 넘었다(편향 0.2928 > 문턱 0.2882).
prereg 0.2: "넘으면 출발 수·띠를 늘린다. 그래도 넘는 구간은 판정에서 표시한다."
같은 프레임 쌍(derived/r38/s0/bias.jsonl 의 그 구간 교차 200·바닥 200)을 줄인 계산의 출발 수(n_local)·처음 반폭(k0)만 바꿔 다시 잰다.
  BnB 는 다시 돌리지 않는다(같은 프레임 쌍의 S_lb 를 그대로). 두 점의 노름 합은 기록에서 낸다: ‖x‖² + ‖y‖² = d_red² + 2 S_red.
  기준 S = max(BnB S_lb, 원래 줄인 계산 S, 모든 판의 줄인 계산 S) — 알려진 가장 좋은 정렬. 오차 = 그 판의 d − 기준 d.
  편향 = 2 × (교차 오차 평균 − 바닥 오차 평균). 문턱은 results/r38_s0_bias.json 의 그 구간 값(그대로).
판: n_local ∈ {1, 8, 32} × k0 ∈ {5, 10, 20}. (1, 5) 는 원래 판의 재현이다(재현 오차를 함께 낸다). 출발: 첫째 = 평균 점 출발, 나머지 무작위(seed 0).
점검: 0.3–0.4 구간(통과한 구간)에도 같은 판을 돌려 편향이 줄기만 하는지 본다.
결정 규칙(돌리기 전에 고정): 0.2–0.3 에서 통과하는 판 가운데 쌍당 시간이 가장 짧은 판을 TM 긴 쪽 < 0.3 칸의 본 계산 매개변수로 쓴다.
  모두 넘으면 그 구간을 판정에서 표시하고 원래 판(1, 5)을 그대로 쓴다.
사용: python r38_s0_bias2.py <ROOT> [smoke]   출력: results/r38_s0_bias2.json (사본 results/_backup/<날짜>/)
  smoke: 판 (1, 5)·(8, 10), 도메인 쌍 둘만 → results/r38_s0_bias2_smoke.json (결정 파일은 건드리지 않음; 원래 판 재현 오차를 보는 연기 시험)"""
import sys, os, json, time
ROOT = sys.argv[1]; SMOKE = len(sys.argv) > 2 and sys.argv[2] == 'smoke'
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)
SET = [(1, 5), (8, 10)] if SMOKE else [(nl, k) for nl in (1, 8, 32) for k in (5, 10, 20)]; BINS2 = ('0.2-0.3', '0.3-0.4')


def red_nl(ms, ia, ib, k0, n_local, iters=10, membudget=3e9):
    """r38_s0.red_pairs 와 같은 규칙(반폭 k0 에서 시작, 띠 가장자리에 닿은 쌍만 반폭 두 배, 열 수까지), 출발 수만 n_local."""
    EA, EB = ms.EA, ms.EB; n, m = EA.shape[1], EB.shape[1]
    ia = torch.as_tensor(np.asarray(ia), device=EA.device); ib = torch.as_tensor(np.asarray(ib), device=EA.device)
    S = EA.new_zeros(len(ia)); todo = torch.arange(len(ia), device=EA.device); kk = k0
    while len(todo):
        lo, w = band_windows(ms.center, m, kk); step = max(1, int(membudget / (n * w * 9 * 8 * 4 * n_local)))
        edge_all = torch.zeros(len(todo), dtype=torch.bool, device=EA.device)
        for c in range(0, len(todo), step):
            sel = todo[c:c + step]; sc, _, M, edge = align_band(EA[ia[sel]], EB[ib[sel]], lo, w, ms.starts(ia[sel], ib[sel], n_local=n_local), iters=iters)
            S[sel] = sc; edge_all[c:c + step] = edge
        if w >= m: break
        todo = todo[edge_all]; kk *= 2
    return S.cpu().numpy()


def sync():
    if DEV == 'cuda': torch.cuda.synchronize()


J = [json.loads(l) for l in open(OUT + '/bias.jsonl')]; J = [j for j in J if j['bin'] in BINS2]
thr = {r['bin']: r['threshold'] for r in json.load(open(RES + '/r38_s0_bias.json'))['rows']}
groups = {}
for k, j in enumerate(J): groups.setdefault((j['a'], j['b']), []).append(k)
if SMOKE: groups = dict(list(groups.items())[:2])
Sres = {s: np.full(len(J), np.nan) for s in SET}; sec = {s: 0.0 for s in SET}; log('bias2', len(J), '프레임 쌍', len(groups), '도메인 쌍', DEV)
for (a, b), ks in groups.items():
    ms = MeanStart(side(a)[0], side(b)[0]); ia = [J[k]['i'] for k in ks]; ib = [J[k]['j'] for k in ks]
    for s in SET:
        sync(); t0 = time.time(); Sv = red_nl(ms, ia, ib, k0=s[1], n_local=s[0]); sync(); sec[s] += time.time() - t0; Sres[s][ks] = Sv
    log('bias2', a, b, len(ks))
Sred = np.array([j['S_red'] for j in J]); Slb = np.array([j['S_lb'] for j in J]); T = np.array([j['d_red'] ** 2 for j in J]) + 2 * Sred
Sbest = np.maximum.reduce([Slb, Sred] + [Sres[s] for s in SET]); dbest = np.sqrt(np.clip(T - 2 * Sbest, 0, None))
kind = np.array([j['kind'] for j in J]); bins = np.array([j['bin'] for j in J]); rows = []
for s in SET:
    ds = np.sqrt(np.clip(T - 2 * Sres[s], 0, None)); err = ds - dbest; r = dict(n_local=s[0], k0=s[1], sec_per_pair=sec[s] / len(J))
    for lab in BINS2:
        c, f = (bins == lab) & (kind == 'cross'), (bins == lab) & (kind == 'floor')
        bias = float(2 * (err[c].mean() - err[f].mean())); r[lab] = dict(bias=bias, threshold=thr[lab], passed=bool(bias < thr[lab]), err_cross=float(err[c].mean()), err_floor=float(err[f].mean()),
                                                                        frac_at_best_cross=float((err[c] < 1e-9).mean()))
    rows.append(r)
d_orig = np.sqrt(np.clip(T - 2 * np.maximum(Slb, Sred), 0, None)); e0 = np.array([j['d_red'] for j in J]) - d_orig
orig = {lab: float(2 * (e0[(bins == lab) & (kind == 'cross')].mean() - e0[(bins == lab) & (kind == 'floor')].mean())) for lab in BINS2}
ok = [r for r in rows if r['0.2-0.3']['passed']]; pick = min(ok, key=lambda r: r['sec_per_pair']) if ok else None
out = dict(written=time.strftime('%F %T'), rows=rows, repro_max_abs_S=float(np.nanmax(np.abs(Sres[(1, 5)] - Sred))), bias_original_definition=orig,
           ref_improved_frac=float((Sbest > np.maximum(Slb, Sred) + 1e-9).mean()),
           decision=dict(n_local=pick['n_local'], k0=pick['k0'], note='TM 긴 쪽 < 0.3 칸의 본 계산 매개변수') if pick else dict(n_local=1, k0=5, note='모든 판이 넘음 → 0.2–0.3 구간을 판정에서 표시, 원래 판'),
           rule='r38_s0_bias2 머리말(prereg 0.2 의 다음 단계). 결정 = 0.2–0.3 에서 통과하는 판 가운데 쌍당 시간이 가장 짧은 판')
jdump(out, RES + ('/r38_s0_bias2_smoke.json' if SMOKE else '/r38_s0_bias2.json'), bk=not SMOKE); log('BIAS2', out['decision'], {f"{r['n_local']}/{r['k0']}": round(r['0.2-0.3']['bias'], 4) for r in rows}, '문턱', round(thr['0.2-0.3'], 4)); log('BIAS2_DONE')
