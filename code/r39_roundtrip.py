"""(논문 구성 3.3 정상 점검 "좌표 왕복", CPU; 등록 = Drive JSON round36_roundtrip_2026-10-20, 돌리기 전) 공간의 좌표와 원자 좌표 사이를 오간 오차.
자료: 줄기 derived/mdcath_5T 336 가운데 default_rng(20261077) 로 고른 60 도메인, 320 K 반복 0 의 저장 프레임(CA, C, O; Å) 앞 100 개.
1 Cα 왕복: (CA, C, O) → 좌표 x = s·R/2 (fs.frames.coordinates_from_backbone) → Cα 자취(ca_trace, 첫 Cα 원점) 와 실제 Cα − Cα_0 의 최대 차. 예측: < 1e-9 Å (항등식).
2 평면 안 원자(C_i, O_i): 변 i 의 틀 R_i 에서 본 국소 좌표 R_iᵀ(C_i − Cα_i), R_iᵀ(O_i − Cα_i) 의 주형(트랜스 변 / 시스 변 따로, 국소 좌표 평균)으로 되살린 C·O 의 RMS 오차.
  주형은 도메인 반(무작위)에서 만들고 나머지 반에서 잰다(두 방향 평균). 서술. 시스 변이 없으면 수만 적는다.
3 평균 점의 사영 오차: 도메인마다 그 프레임들의 평균 점 m (fs.mean_point) 을 항마다 가장 가까운 실제 항 s·R/2 (R = m_i 의 극분해 회전, s = 2 tr(Rᵀ m_i)/3) 로 사영.
  사영 오차 Σ‖m_i − P(m_i)‖² 를 앙상블 흩어짐 mean_f d_fix(E_f, m)² 와 Σ‖m_i‖² 에 견준 비, 그리고 사영한 평균 점의 Cα 자취와 맞춘 실제 Cα 평균의 RMSD. 서술.
사용: python r39_roundtrip.py <ROOT> 출력: results/r39_roundtrip.json"""
import os
os.environ['CUDA_VISIBLE_DEVICES'] = ''
for k_ in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ[k_] = '2'
import sys, json, glob, time, shutil
import numpy as np
ROOT = sys.argv[1]; sys.path[:0] = [ROOT + '/framespace_pkg']
import torch
import framespace as fs
from framespace.frames import coordinates_from_backbone, frames_of, ca_trace, CIS_MAX, L0
NDOM = 60; NF = 100


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def kabsch_rmsd(P, Q):
    P = P - P.mean(0); Q = Q - Q.mean(0); U, S, Vt = np.linalg.svd(P.T @ Q); d = np.sign(np.linalg.det(U @ Vt)); D = np.diag([1, 1, d])
    return float(np.sqrt(((P @ (U @ D @ Vt) - Q) ** 2).sum(1).mean()))


def proj(m):
    """항마다 가장 가까운 s·R/2."""
    U, S, Vt = np.linalg.svd(m); d = np.sign(np.linalg.det(U @ Vt)); D = np.zeros_like(m); D[:, 0, 0] = 1; D[:, 1, 1] = 1; D[:, 2, 2] = d
    Rm = U @ D @ Vt; s = 2 * np.einsum('nij,nij->n', Rm, m) / 3; return s[:, None, None] * Rm / 2


def main():
    stem = sorted(os.path.basename(f)[:-4] for f in glob.glob(ROOT + '/derived/mdcath_5T/*.npz')); rng = np.random.default_rng(20261077); sel = sorted(rng.choice(stem, NDOM, replace=False))
    half = set(rng.choice(sel, NDOM // 2, replace=False)); out = dict(written=time.strftime('%F %T'), n_dom=len(sel), doms=list(sel), per_dom={})
    loc = {}; ca_err = []
    for d in sel:
        z = np.load(f'{ROOT}/derived/mdcath_5T/{d}.npz'); X = z['320_0_X'][:NF].astype(np.float64); CA, C, O = X[:, :, 0], X[:, :, 1], X[:, :, 2]
        x = coordinates_from_backbone(CA, C, O); R, s = frames_of(x); l = s * L0
        rec = ca_trace(x); ca_err.append(float(np.abs(rec - (CA - CA[:, :1])).max()))
        cl = np.einsum('fnji,fnj->fni', R, C[:, :-1] - CA[:, :-1]); ol = np.einsum('fnji,fnj->fni', R, O[:, :-1] - CA[:, :-1])
        loc[d] = (cl.reshape(-1, 3), ol.reshape(-1, 3), (l < CIS_MAX).reshape(-1), R.reshape(-1, 3, 3), (C[:, :-1] - CA[:, :-1]).reshape(-1, 3), (O[:, :-1] - CA[:, :-1]).reshape(-1, 3))
        E = torch.as_tensor(x); m = fs.mean_point(E).numpy(); P = proj(m); A = fs.align_to_reference(E, torch.as_tensor(m)).numpy()
        spread = float(((A - m[None]) ** 2).sum((1, 2, 3)).mean()); perr = float(((m - P) ** 2).sum())
        ca_m = ca_trace(P[None])[0]; ca_al = ca_trace(A).mean(0)
        out['per_dom'][d] = dict(L=int(X.shape[1]), n_cis_edges=int((l < CIS_MAX).sum()), ca_max_err=ca_err[-1], proj_err=perr, proj_over_spread=perr / spread, proj_over_norm=perr / float((m ** 2).sum()),
                                 proj_ca_rmsd=kabsch_rmsd(ca_m, ca_al), s_proj_mean=float((2 * np.linalg.norm(P[:, :, 0], axis=-1)).mean()), s_frames_mean=float(s.mean()))
        log('왕복', d, X.shape[1], 'Cα', f'{ca_err[-1]:.2e}', '사영/흩어짐', round(perr / spread, 4), 'Cα RMSD', round(out['per_dom'][d]['proj_ca_rmsd'], 3))
    out['ca_max_err'] = float(max(ca_err)); out['ca_pred_pass'] = bool(max(ca_err) < 1e-9)
    res = {}
    for kind, mk in (('trans', lambda c: ~c), ('cis', lambda c: c)):
        errs = {'C': [], 'O': []}; n_test = 0
        for train, test in ((half, set(sel) - half), (set(sel) - half, half)):
            cs = [loc[d][0][mk(loc[d][2])] for d in train]; os_ = [loc[d][1][mk(loc[d][2])] for d in train]
            if sum(len(c) for c in cs) == 0: continue
            tc = np.concatenate(cs).mean(0); to = np.concatenate(os_).mean(0)
            for d in test:
                k = mk(loc[d][2]); Rr = loc[d][3][k]
                if not k.any(): continue
                pc = np.einsum('nij,j->ni', Rr, tc); po = np.einsum('nij,j->ni', Rr, to); n_test += int(k.sum())
                errs['C'].append(((pc - loc[d][4][k]) ** 2).sum(1)); errs['O'].append(((po - loc[d][5][k]) ** 2).sum(1))
        res[kind] = dict(n_edges_test=n_test, **({f'rms_{a}': float(np.sqrt(np.concatenate(v).mean())) for a, v in errs.items() if v}), **({f'p95_{a}': float(np.sqrt(np.percentile(np.concatenate(v), 95))) for a, v in errs.items() if v}))
    out['inplane'] = res
    pd_ = out['per_dom'].values(); out['proj_summary'] = {k: dict(median=float(np.median([v[k] for v in pd_])), max=float(np.max([v[k] for v in pd_]))) for k in ('proj_over_spread', 'proj_over_norm', 'proj_ca_rmsd')}
    out['rule'] = '정상 점검(논문 3.3). Cα 왕복만 예측(< 1e-9 Å, 항등식); 평면 안 원자와 평균 점 사영은 서술. 도메인 60 (seed 20261077), 320 K 반복 0 의 앞 100 프레임.'
    p = ROOT + '/results/r39_roundtrip.json'; json.dump(out, open(p + '.tmp', 'w'), ensure_ascii=False, indent=1); os.replace(p + '.tmp', p)
    b = ROOT + '/results/_backup/' + time.strftime('%Y-%m-%d'); os.makedirs(b, exist_ok=True); shutil.copy(p, b + '/' + time.strftime('%H%M%S_') + 'r39_roundtrip.json')
    log('ROUNDTRIP', 'Cα 최대', out['ca_max_err'], '평면 안', res, '사영', out['proj_summary']); log('ROUNDTRIP_DONE')


if __name__ == '__main__':
    main()
