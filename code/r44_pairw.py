"""(round44, 등록 = Drive results/prereg_2026-10-28_addendum2.json "r44_pairw", 2026-10-28 UTC) 서열 동일성 효과의 쌍 무게 추정치.
등록 판(r38_c.analyse_sf, r38_p3.analyse)은 상과마다 같은 무게(쌍 무게 1/상과 쌍 수)라 추정 대상이 "상과 평균 효과" 다. 외부 검토(2026-10-26)대로
쌍마다 같은 무게(1)인 판을 함께 낸다. 구간은 같은 상과 Exp(1) 군집 부트스트랩 2,000 (무게만 1 로). 결과 자료는 저장된 것(pairs_out·anm_pairs)을 그대로 쓴다.
판: C (확인, 4,152 쌍) C-1 (z = control) · C-2 (z = control, r_ANM); H-풀 (11,991 쌍) P3 주 판 (z = control) · C-2 (z = control, r_ANM). 비교로 상과 무게 판을 같은 코드로 다시 낸다.
사용: python r44_pairw.py <ROOT>   출력: results/r44_pairw.json (사본 results/_backup/<날짜>/)"""
import os, sys
ROOT = sys.argv[1]
_src = open(ROOT + '/code/r38_c.py').read(); _src = _src[:_src.index("ALLC = ['pairs'")]
sys.argv = ['r38_c', ROOT, 'none', '1']; exec(_src)   # merged, analyse_sf, ranks, pcorr, NB, OUT, P3OUT, RES, jdump, log
SEEDW = 20261097


def analyse_w(rows, key='main', y='r_dyn', x='ident', z=('control',), seed=SEEDW, B=NB, pair_weight=True):
    """analyse_sf 와 같은 점추정·부트스트랩. pair_weight=True 면 쌍 무게 1 (상과 Exp(1) 무게만 곱함), False 면 등록 판(1/상과 쌍 수)."""
    rows = [r for r in rows if r.get(key) and r[key].get(y) is not None and r[key].get(x) is not None and all(r[key].get(c) is not None for c in z)]
    if len(rows) < 10: return dict(n=len(rows))
    Y = ranks(np.array([r[key][y] for r in rows])); X_ = ranks(np.array([r[key][x] for r in rows])); Zs = [ranks(np.array([r[key][c] for r in rows])) for c in z]
    H = np.array([r['H'] for r in rows]); sfl = sorted(set(H)); si = np.array([sfl.index(h) for h in H]); cnt = np.bincount(si)
    base = np.ones(len(rows)) if pair_weight else 1.0 / cnt[si]
    est = pcorr(Y, X_, Zs, base); rng = np.random.default_rng(seed); bs = []
    for _ in range(B):
        v = rng.exponential(1, len(sfl)); bs.append(pcorr(Y, X_, Zs, base * v[si]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(n=len(rows), n_sf=len(sfl), est=float(est), lo=float(lo), hi=float(hi), se=float(np.std(bs)))


rc = merged(OUT, OUT + '/anm_pairs.jsonl'); rh = merged(P3OUT, OUT + '/anm_pairs_hpool.jsonl'); out = dict(written=time.strftime('%F %T'))
for nm, rows, z, s in (('C_C1', rc, ('control',), 1), ('C_C2', rc, ('control', 'r_anm'), 2), ('Hpool_P3', rh, ('control',), 3), ('Hpool_C2', rh, ('control', 'r_anm'), 4)):
    out[nm] = dict(pair_weighted=analyse_w(rows, z=z, seed=SEEDW + s), sf_weighted=analyse_w(rows, z=z, seed=SEEDW + 10 + s, pair_weight=False))
    log('R44_PAIRW', nm, {k: (round(v.get('est', np.nan), 3), round(v.get('lo', np.nan), 3), round(v.get('hi', np.nan), 3)) for k, v in out[nm].items()})
c1 = out['C_C1']['pair_weighted']; hp = out['Hpool_P3']['pair_weighted']
out['predictions'] = dict(p1_C1_pair_smaller=bool(c1['est'] < out['C_C1']['sf_weighted']['est']), p2_C1_pair_lo_pos=bool(c1['lo'] > 0), p3_Hpool_pair_includes0=bool(hp['lo'] <= 0 <= hp['hi']))
out['rule'] = '등록 r44_pairw: 서술 + 예측 셋(P1 C-1 쌍 무게 점추정 < 상과 무게, P2 C-1 쌍 무게 하한 > 0, P3 H-풀 P3 쌍 무게 구간이 0 포함). 상과 무게 판은 같은 코드로 다시 낸 비교값.'
jdump(out, RES + '/r44_pairw.json', bk=True); log('R44_PAIRW_DONE', out['predictions'])
