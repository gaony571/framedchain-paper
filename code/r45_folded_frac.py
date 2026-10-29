"""(서술 셈, 등록 없음; 2026-10-29 11:31–11:37 UTC Colab 에서 /content/folded_frac.py 로 실행한 원본) P2 67 도메인의 MD 320 K 전체 궤적에서 접힌 프레임(Q ≥ 0.75) 비율과
aSAMt 접힌 표본 수(compare_D 의 n_model / 250). Q 정의 = r38_s0.md_load 의 Qf (저장된 천연 접촉 ii·jj·dn, 1.2 배). 출력: results/r45_folded_frac.json"""
import json, numpy as np, h5py, time
R = '/content/drive/MyDrive/atlas_framed_curves'
doms = open(f'{R}/derived/r38/q3_p2_67_ordered.txt').read().split(); out = {}
for d in doms:
    z = np.load(f'{R}/derived/r38/md/{d}.npz'); ii, jj, dn, ica = z['ii'], z['jj'], z['dn'], z['imd'][:, 0]
    with h5py.File(f'{R}/raw/mdcath/mdcath_dataset_{d}.h5', 'r') as h:
        g = h[d]['320']; tot = fold = 0
        for r in sorted(g.keys(), key=int):
            X = g[r]['coords'][:][:, ica]; q = (np.linalg.norm(X[:, ii] - X[:, jj], axis=-1) < 1.2 * dn).mean(1); tot += len(q); fold += int((q >= 0.75).sum())
    c = json.load(open(f'{R}/derived/r38/compare_D/{d}.json')).get('cells', {}).get('aSAMt320-MD320|full|0.75', {})
    out[d] = dict(md_frames=tot, md_folded=fold, md_frac=fold / tot, asamt_folded=c.get('n_model'), asamt_frac=(c.get('n_model') or 0) / 250)
    print(time.strftime('%H:%M:%S'), d, round(fold / tot, 3), c.get('n_model'), flush=True)
f = lambda k: [v[k] for v in out.values()]
S = dict(written=time.strftime('%F %T'), n=len(out), md_frac_median=float(np.median(f('md_frac'))), md_frac_min=float(np.min(f('md_frac'))),
         asamt_frac_median=float(np.median(f('asamt_frac'))), asamt_frac_min=float(np.min(f('asamt_frac'))), asamt_folded_median=float(np.median(f('asamt_folded'))),
         n_asamt_below_90pct=int(sum(x < 0.9 for x in f('asamt_frac'))), n_md_below_90pct=int(sum(x < 0.9 for x in f('md_frac'))), rows=out,
         note='서술 셈(등록 없음): P2 67 도메인, 320 K 전체 궤적, Q 정의 = r38_s0.md_load 와 같음(저장된 천연 접촉). aSAMt = compare_D 의 n_model / 250.')
json.dump(S, open(f'{R}/results/r45_folded_frac.json', 'w'), ensure_ascii=False, indent=1); print('DONE', {k: S[k] for k in S if k not in ('rows',)})
