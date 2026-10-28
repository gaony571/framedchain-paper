"""(round44, 등록 밖 서술 — 결과를 본 뒤의 구현 점검) r44_corr 의 서열 정렬 대응이 π̄ 짝을 담는 비율을 π̄ 동일성 구간별로.
r44_corr 의 서열 정렬 대응이 π̄ 짝을 담는 비율의 중앙이 0 이라, 색인 어긋남(버그)인지 보려고 Colab 노트북 셀 30 에서 돌렸다(2026-10-28 06:3x UTC, 읽기만).
원고는 이 표를 쓰지 않는다: 같은 물음을 π̄ 와 독립인 기준(TM-align)으로 다시 등록해 낸 r44_corr2.py (r44_agree) 를 쓴다.
아래 "셀 30 그대로" 줄부터 끝까지가 그 셀의 글자 그대로이고, 그 위의 머리는 이 파일만으로 돌리려고 저장소를 꾸릴 때 더한 것이다.
사용: python r44_seqid_cell30.py <ROOT>   입력: derived/r44/corr/pairs_v2_s*.jsonl   출력: 표준 출력"""
import os, sys, json
import numpy as np
R = sys.argv[1]
# ---- 셀 30 그대로 ----
# 30. (읽기만, 구현 점검) 서열 정렬 대응의 π̄ 포함 비율이 동일성에 따라 오르는지 — 색인 버그면 높은 동일성에서도 0 근처
import glob
rows44 = [json.loads(l) for f in sorted(glob.glob(f'{R}/derived/r44/corr/pairs_v2_s*.jsonl')) for l in open(f) if l.strip()]
def g_(r, c, k): return (r.get(c) or {}).get(k)
X = [(g_(r, 'pibar', 'ident'), g_(r, 'seq', 'frac_pibar_in'), g_(r, 'tm', 'frac_pibar_in'), g_(r, 'seq', 'ident'), g_(r, 'seq', 'n_pairs'), g_(r, 'pibar', 'n_pairs'), r.get('tm_med')) for r in rows44]
X = [x for x in X if None not in x[:3]]
a = np.array([x[:4] + x[4:6] + (x[6],) for x in X], float); print('n', len(a))
for lo, hi in ((0, .1), (.1, .15), (.15, .2), (.2, .3), (.3, .5), (.5, 1.01)):
    m = (a[:, 0] >= lo) & (a[:, 0] < hi)
    if m.sum(): print(f'π̄ 동일성 [{lo},{hi}) n={m.sum():4d}  seq⊃π̄ 중앙 {np.median(a[m,1]):.3f}  tm⊃π̄ 중앙 {np.median(a[m,2]):.3f}  seq 동일성 중앙 {np.median(a[m,3]):.3f}  TM 중앙 {np.median(a[m,6]):.2f}')
print('seq 짝 수 중앙', np.median(a[:, 4]), 'π̄ 짝 수 중앙', np.median(a[:, 5]))
