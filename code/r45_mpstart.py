"""(round45, 등록 = Drive results/prereg_2026-10-29_addendum8.json "r45_mpstart") 평균 점이 출발 프레임에 매이지 않는가 — 실행본 코드(framespace 0.2.0 의 fs.mean_point)로 실제 앙상블에서.
앙상블 = 확인 실험 C 의 도메인 앙상블(r38_c.dom_build: 320 K 뒤 절반, Q ≥ 0.75, 반복마다 40–100 접힌 프레임, 반복 ≥ 2; 같은 seed). dom_build 안의 fs.mean_point 를
감싸 그 입력 E (float64, 합친 프레임)를 받는다. 도메인 6 = C 도메인(정렬) 가운데 저장 평균 점이 있는 것에서 default_rng(20261116).choice(6).
출발 9 = 프레임 0 (실행본의 기본) + default_rng(20261117) 로 1..F−1 에서 8. 출발마다 fs.mean_point(E, start=k) (tol 10⁻¹⁰ (1 + ‖m‖), 최대 100 회).
잼: 출발 0 판과 저장 평균 점(derived/r38/c/dom/<d>.npz 'm')의 최대 절대 차(재현 점검); 출발 0 판과 각 출발 판 사이 거리 = 최적 회전 하나로 맞춘 뒤
‖m₀ − U m_k‖ 를 직접 셈(‖a‖² + ‖b‖² − 2S 의 상쇄 오차를 피함), 정규화 D̂ = 그 값 / √(‖m₀‖² + ‖m_k‖²); 수렴 여부·반복 수. 서술(판정 규칙 없음).
사용: python r45_mpstart.py <ROOT>   출력: results/r45_mpstart.json (사본 results/_backup/<날짜>/)"""
import os, sys
ROOT = sys.argv[1]
_src = open(ROOT + '/code/r38_c.py').read(); _src = _src[:_src.index("ALLC = ['pairs'")]
sys.argv = ['r38_c', ROOT, 'none', '1']; exec(_src)   # r38_c (그 안에서 r38_p3): dom_build (C 판), dom_small, fs, torch, np, log, jdump, RES
SEED_DOM, SEED_START, NDOM, NSTART = 20261116, 20261117, 6, 9
_orig_mean_point = fs.mean_point; _cap = {}
def _capture(E, *a, **k):
    _cap['E'] = E; return _orig_mean_point(E, *a, **k)


def fit(m0, mk):
    """m_k 를 m₀ 에 최적 회전 하나로 맞춘 뒤 ‖m₀ − U m_k‖ (직접)."""
    M = np.einsum('iab,icb->ac', m0, mk); W, S, Vt = np.linalg.svd(M); D = np.diag([1.0, 1.0, np.sign(np.linalg.det(W @ Vt))]); U = W @ D @ Vt
    return float(np.sqrt(((m0 - np.einsum('ac,icb->iab', U, mk)) ** 2).sum()))


def main():
    cdom = sorted(json.load(open(OUT + '/pairs.json'))['domains']); ok = [d for d in cdom if dom_small(d)[0] is not None]
    pick = sorted(np.random.default_rng(SEED_DOM).choice(ok, NDOM, replace=False).tolist()); rows = []
    for d in pick:
        fs.mean_point = _capture
        try: r, info = dom_build(d)
        finally: fs.mean_point = _orig_mean_point
        E = torch.as_tensor(_cap['E']).to(torch.float64); F = len(E)
        starts = [0] + sorted(np.random.default_rng(SEED_START).choice(np.arange(1, F), NSTART - 1, replace=False).tolist())
        ms, its = [], []
        for k in starts:
            m, inf = _orig_mean_point(E, start=int(k), return_info=True); ms.append(np.asarray(m, np.float64)); its.append(dict(start=int(k), converged=bool(inf['converged']), iters=int(inf['iters'])))
        stored = np.load(f'{OUT}/dom/{d}.npz')['m'].astype(np.float64)
        dist = [fit(ms[0], mk) for mk in ms[1:]]; nrm = [float(np.sqrt((ms[0] ** 2).sum() + (mk ** 2).sum())) for mk in ms[1:]]
        rows.append(dict(dom=d, n_frames=F, n_bonds=int(E.shape[1]), starts=its, stored_max_abs_diff=float(np.abs(ms[0] - stored).max()),
                         max_d=float(max(dist)), max_dhat=float(max(x / y for x, y in zip(dist, nrm))), norm_m0=float(np.sqrt((ms[0] ** 2).sum()))))
        log('r45_mpstart', d, F, rows[-1]['max_d'], rows[-1]['max_dhat'], rows[-1]['stored_max_abs_diff'])
    out = dict(written=time.strftime('%F %T'), rows=rows, max_d=max(r['max_d'] for r in rows), max_dhat=max(r['max_dhat'] for r in rows),
               stored_max_abs_diff=max(r['stored_max_abs_diff'] for r in rows), all_converged=all(s['converged'] for r in rows for s in r['starts']),
               rule='등록 r45_mpstart: 서술. 실행본 fs.mean_point, 도메인 6 × 출발 9, 직접 셈한 최적 회전 뒤 거리.')
    jdump(out, RES + '/r45_mpstart.json', bk=True); log('R45_MPSTART_DONE', out['max_d'], out['max_dhat'], out['all_converged'])


if __name__ == '__main__':
    main()
