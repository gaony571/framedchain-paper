
"""mdcath_5T 의 접힌 프레임(Q ≥ 0.75) 퍼짐 (2026-10-17). 도메인마다 derived/mdcath_5T_folded/<dom>.npz:
T_spread_rep (5, ne) 변 단위 프레임 퍼짐(그 온도의 접힌 프레임으로 낸 평균 점에 맞춘 좌표), T_rmsf2_rep (5, L) Cα RMSF² (평균 구조에 Kabsch),
T_wrmsf2_rep (5, L) 창(±3 잔기) Cα RMSF² (창마다 따로 Kabsch), T_nfold (5,), T_fracfold (5,), T_meanQ (5,), ca_ref, dssp320 (최빈), seq."""
import sys, os, json, time
ROOT = sys.argv[1]; IN, OUTD = ROOT + '/derived/mdcath_5T', ROOT + '/derived/mdcath_5T_folded'; os.makedirs(OUTD, exist_ok=True)
sys.path[:0] = [ROOT + '/framespace_pkg', ROOT + '/framespace_wip']
import numpy as np, torch, framespace as fs
st = json.load(open(ROOT + '/results/mdcath_5T_setup.json')); QT = 0.75
def kabsch_batch(P, Q):
    """P (F, n, 3) 를 Q (n, 3) 에 맞춘 좌표."""
    Pc = P - P.mean(1, keepdims=True); Qc = Q - Q.mean(0)
    H = np.einsum('fni,nj->fij', Pc, Qc); U, S, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(np.einsum('fij,fjk->fik', U, Vt)))
    D = np.ones((len(P), 3)); D[:, 2] = d; Rm = np.einsum('fij,fj,fjk->fik', U, D, Vt)
    return np.einsum('fni,fij->fnj', Pc, Rm) + Q.mean(0)
def rmsf2(CA):
    ref = CA[0]
    for _ in range(3):
        A = kabsch_batch(CA, ref); ref = A.mean(0)
    return ((A - ref) ** 2).sum(-1).mean(0)
def wrmsf2(CA, w=3):
    F, L, _ = CA.shape; out = np.zeros(L)
    for i in range(L):
        lo, hi = max(0, i - w), min(L, i + w + 1); W = CA[:, lo:hi]; ref = W[0]
        for _ in range(2):
            A = kabsch_batch(W, ref); ref = A.mean(0)
        out[i] = ((A[:, i - lo] - ref[i - lo]) ** 2).sum(-1).mean()
    return out
t0 = time.time(); n = 0
for d in st['domains']:
    out = f'{OUTD}/{d}.npz'
    if os.path.exists(out): continue
    try:
        z = np.load(f'{IN}/{d}.npz'); res = dict(seq=z['seq'], ca_ref=z['320_0_X'][0, :, 0])
        ds = z['320_0_dssp']; DS = ' HBEGITSP-'; res['dssp320'] = np.array([DS[np.bincount(ds[:, k], minlength=10).argmax()] for k in range(ds.shape[1])])
        for T in st['temps']:
            reps = sorted({k.split('_')[1] for k in z.files if k.startswith(T + '_') and k.endswith('_X')}, key=int)
            Xs, cut, nf, ff, mq = [], [0], [], [], []
            for r in reps:
                X = z[f'{T}_{r}_X']; Q = z[f'{T}_{r}_Q'].astype(float)[::st['frame_step']][:len(X)]; m = Q >= QT
                Xs.append(X[m]); cut.append(cut[-1] + int(m.sum())); nf.append(int(m.sum())); ff.append(float(m.mean())); mq.append(float(Q.mean()))
            XA = np.concatenate(Xs); res[T + '_nfold'] = np.array(nf); res[T + '_fracfold'] = np.array(ff); res[T + '_meanQ'] = np.array(mq)
            if len(XA) < 20: continue
            E = fs.coordinates_from_backbone(torch.as_tensor(XA[:, :, 0]), torch.as_tensor(XA[:, :, 1]), torch.as_tensor(XA[:, :, 2]))
            xbar = fs.mean_point(E[:600]); Ea = fs.align_to_reference(E, xbar); dev = ((Ea - xbar[None]) ** 2).sum((-1, -2)).numpy()
            sp, r2, w2 = [], [], []
            for k in range(len(reps)):
                a, b = cut[k], cut[k + 1]
                if b - a < 10: sp.append(np.full(dev.shape[1], np.nan)); r2.append(np.full(XA.shape[1], np.nan)); w2.append(np.full(XA.shape[1], np.nan)); continue
                sp.append(dev[a:b].mean(0)); CA = XA[a:b, :, 0].astype(float); r2.append(rmsf2(CA)); w2.append(wrmsf2(CA))
            res[T + '_spread_rep'], res[T + '_rmsf2_rep'], res[T + '_wrmsf2_rep'] = np.array(sp), np.array(r2), np.array(w2)
        np.savez_compressed(out, **res); n += 1
        if n % 10 == 0: print(time.strftime('%H:%M:%S'), n, round((time.time() - t0) / n, 1), 's/dom', flush=True)
    except Exception as e:
        print('ERR', d, repr(e)[:200], flush=True)
print('FINISHED', n, flush=True); open(ROOT + '/results/logs/mdcath_5T_folded.done', 'w').write(time.strftime('%F %T'))
