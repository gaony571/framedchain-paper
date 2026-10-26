"""(prereg round36 2 절 R6 주 3 사건, 서술; 등록 = Drive JSON round36_r6_2026-10-20, 돌리기 전) 같은 단백질의 자연 삽입·결실 3 사건의 앙상블 D.
사건(ATLAS 끼리, 반복 3 × 100 ns): 7la6_A/4r7q_A, 1u09_A/4wyw_A, 1ifg_A/1n8o_E. 궤적 = raw/atlas_r6·raw/atlas_af82 의 <사슬>_analysis.zip 안 R1–R3.xtc (100 ps 간격 1,001 프레임; 내려받기 없음).
계산은 r38_b6.one() 그대로(R10a′ 와 같음): 실행 = 반복, 접힌 프레임(첫 프레임 천연 접촉, Q ≥ 0.75) 가운데 N = 100 무작위(사슬 이름 seed), 프레임 쌍 2,000, 줄인 계산; Δ²_D, 짝·빈칸 몫, D_matched, π̄·TM-align 고정 대응, 정적 BnB.
  더함: 같은 서열 자기 비교의 빈칸 몫(R10a′ 와 같은 방식).
야생형(관례) = 짧은 사슬 a. 삽입 자리 = 두 서열 정렬(difflib)에서 끝이 아닌 삽입·결실 덩어리가 a 에 들어가는 자리 k (a 의 잔기 k−1 과 k 사이);
  창 = a 의 변 k−4 … k+2 (자리를 걸치는 변 k−1 을 가운데 둔 7 변). 야생형 창 퍼짐 = a 의 실행별 변별 퍼짐(맞춘 좌표, 320 K 대신 그 사슬 접힌 프레임 앞 600 의 평균 점에 맞춤)의 창 평균,
  상대 값 = 창 평균 / 모든 변 중앙값.
서술(prereg 2 절 예측 1·4): 관행 열(TM-align 짝 잔기 Δ²/√n̄)과 Δ²_D/√n̄ 의 R0 띠 대비 위치, 짝 몫(변당)과 야생형 창 퍼짐의 순위 상관(사건 3 개라 서술; ρ ≤ −0.5 예측).
보조 5 사건(DynamicPDB 궤적 하나)은 이번 범위 밖이다.
고침 1 (2026-10-20 UTC, 계산 전): R10a′ 고침 1 의 빈칸 자리 나눔(gap_regions: π̄ 기준 끝 겹침·건너뜀 창·중심)을 그대로 붙임. 값은 바뀌지 않음(서술 더함).
사용: python r39_r6.py <ROOT> ens|stats 출력: derived/r39/r6/<a>__<b>.json, results/r39_r6.json"""
import sys, os, json
ROOT = sys.argv[1]; PARTS_R6 = sys.argv[2].split(',')
_b = open(ROOT + '/code/r38_b6.py').read(); _b = _b[:_b.index("\ndef part_ens():")]
sys.argv = ['r38_b6', ROOT, 'none']; exec(_b)
import zipfile, difflib, shutil, mdtraj as md
OUTR6 = ROOT + '/derived/r39/r6'; os.makedirs(OUTR6, exist_ok=True); TMP6 = '/content/r6_tmp'; os.makedirs(TMP6, exist_ok=True)
EVENTS = [('7la6_A', '4r7q_A'), ('1u09_A', '4wyw_A'), ('1ifg_A', '1n8o_E')]
_r10 = open(ROOT + '/code/r38_r10a.py').read(); exec(_r10[_r10.index('_rpm0 = runpair_means'):_r10.index('def cath_H():')]) # 고침 1 (2026-10-20): R10a′ 고침 1 의 runpair_means·red_pairs·gap_regions 를 그대로


def zpath(c):
    for sub in ('atlas_r6', 'atlas_af82'):
        f = f'{ROOT}/raw/{sub}/{c}_analysis.zip'
        if os.path.exists(f): return f
    raise FileNotFoundError(c)


_AT = {}
def atlas_load(c):
    """[반복 (F, L, 3, 3) Å (CA, C, O)], Qf, L, 서열."""
    if c in _AT: return _AT[c]
    z = zipfile.ZipFile(zpath(c)); z.extract(f'{c}.pdb', TMP6); top = md.load(f'{TMP6}/{c}.pdb').topology; seq, idx = at_atoms(top); L = len(seq); reps = []
    for r in (1, 2, 3):
        n = f'{c}_R{r}.xtc'; z.extract(n, TMP6); t = md.load(f'{TMP6}/{n}', top=top, atom_indices=idx.ravel()); os.remove(f'{TMP6}/{n}'); reps.append(t.xyz.reshape(t.n_frames, L, 3, 3) * 10)
    ca0 = reps[0][0, :, 0]; ii, jj = np.triu_indices(L, 4); dn = np.linalg.norm(ca0[ii] - ca0[jj], axis=-1); m_ = dn < 8.0; ii, jj, dn = ii[m_], jj[m_], dn[m_]
    Qf = lambda X: (np.linalg.norm(X[:, ii, 0] - X[:, jj, 0], axis=-1) < 1.2 * dn).mean(1)
    _AT[c] = (reps, Qf, L, seq); return _AT[c]


def dom_runs(d, T='320', W='full', temps=('320',)):
    """ATLAS 사슬판: [(반복 번호, (N, L, 3, 3))] 와 L. 접힌 프레임 ≥ N 인 반복만(전체 창)."""
    key = (d, T, W)
    if key not in _runs:
        reps, Qf, L, _ = atlas_load(d); rng = seeded('r39_r6_runs', d); out = []
        for r, X in enumerate(reps):
            idx = np.where(Qf(X) >= QT)[0]
            if len(idx) >= N: out.append((r, X[np.sort(rng.choice(idx, N, replace=False))]))
        _runs[key] = (out, L)
    return _runs[key]


def site(sa, sb):
    """a(짧은 쪽) 에서 끝이 아닌 삽입·결실 덩어리의 자리 k 들."""
    ops = difflib.SequenceMatcher(None, sa, sb, autojunk=False).get_opcodes(); ks = []
    for tag, i1, i2, j1, j2 in ops:
        if tag != 'equal' and 0 < i1 and i2 < len(sa) and (i2 - i1) != (j2 - j1): ks.append(dict(tag=tag, k=i1, a=[i1, i2], b=[j1, j2], a_seq=sa[i1:i2], b_seq=sb[j1:j2]))
    return ks, ops


def wt_window(a, k):
    R, L = dom_runs(a); ne = L - 1; reps, Qf, _, _ = atlas_load(a)
    allf = np.concatenate([X[Qf(X) >= QT] for X in reps])[:600]; ref = fs.mean_point(frames(allf))
    sp = []
    for _, X in R:
        A = fs.align_to_reference(frames(X), ref).numpy().reshape(len(X), -1); sp.append(((A - A.mean(0)) ** 2).reshape(len(X), ne, 9).sum(2).mean(0))
    s = np.mean(sp, 0); w = np.arange(max(0, k - 4), min(ne, k + 3)); return dict(window_edges=[int(w[0]), int(w[-1])], spread=float(s[w].mean()), spread_rel=float(s[w].mean() / np.median(s)), median_edge_spread=float(np.median(s)))


def part_ens():
    for a, b in EVENTS:
        fo = f'{OUTR6}/{a}__{b}.json'
        if os.path.exists(fo): continue
        try:
            sa, sb = atlas_load(a)[3], atlas_load(b)[3]
            if len(sa) > len(sb): a, b, sa, sb = b, a, sb, sa
            ks, ops = site(sa, sb); pr = dict(a=a, b=b, n=len(sa), m=len(sb), TM=1.0, level='R6', Ha='R6', Hb='R6', bin='r6', ln=abs(len(sb) - len(sa)),
                                              ln_internal=int(sum(abs((x['b'][1] - x['b'][0]) - (x['a'][1] - x['a'][0])) for x in ks)))
            _CAP.clear(); r = one(pr)
            if 'skip' not in r:
                cr, fa, fb = _CAP[0], _CAP[1], _CAP[2]; r.update(gap_cross=float(cr['gv'].mean()), gap_self_A=float(fa['gv'].mean()), gap_self_B=float(fb['gv'].mean()))
                r['sites'] = ks; r['wt'] = [dict(site_k=x['k'], **wt_window(a, x['k'])) for x in ks]
                try: r['gap_regions'] = gap_regions(cr, fa, fb)
                except Exception as e: r['gap_regions'] = dict(err=repr(e)[:200])
            jdump(r, fo); log('r6', a, b, r.get('skip', ''), '자리', [x['k'] for x in ks], round(r.get('sec', 0)), 's')
        except Exception as e:
            traceback.print_exc(); log('ERR r6', a, b, repr(e)[:200])
        for f in glob.glob(TMP6 + '/*'): os.remove(f) if os.path.isfile(f) else shutil.rmtree(f, ignore_errors=True)
    log('R6_ENS_DONE')


def part_stats():
    from scipy.stats import spearmanr
    J = [json.load(open(f)) for f in sorted(glob.glob(OUTR6 + '/*.json'))]; J = [j for j in J if 'skip' not in j]; band = R0_BAND; rows = []
    for j in J:
        tm = j.get('tmalign') or {}; ws = j['wt'][0]['spread'] if j.get('wt') else None
        rows.append(dict(a=j['a'], b=j['b'], ln=j['ln'], ln_internal=j.get('ln_internal'), n_sites=len(j.get('sites', [])), d2n_D=j['D']['d2n'], d2n_tmalign=tm.get('d2n'), d2n_matched=(j.get('D_matched') or {}).get('d2n'),
                         pos_D=('위' if j['D']['d2n'] > band[1] else ('안' if j['D']['d2n'] >= band[0] else '아래')), pos_tmalign=(None if tm.get('d2n') is None else ('위' if tm['d2n'] > band[1] else ('안' if tm['d2n'] >= band[0] else '아래'))),
                         pair_share_per_edge=j['pair_share_per_edge'], gap_share_per_edge=j['gap_share_per_edge'], gap_excess=j['gap_cross'] - (j['gap_self_A'] + j['gap_self_B']) / 2, gap_target=0.75 * j['ln'],
                         wt_window_spread=ws, wt_window_spread_rel=j['wt'][0]['spread_rel'] if j.get('wt') else None, static_D=j['static']['D'], static_certified=j['static']['certified'],
                         gap_regions={k: (j.get('gap_regions') or {}).get(k) for k in ('overhang', 'jump_window', 'core', 'forced_min', 'n_overhang_edges', 'jumps')}))
    x = [r['pair_share_per_edge'] for r in rows if r['wt_window_spread'] is not None]; y = [r['wt_window_spread'] for r in rows if r['wt_window_spread'] is not None]
    out = dict(written=time.strftime('%F %T'), n=len(rows), rows=rows, band=list(band), rho_pairshare_vs_wtspread=float(spearmanr(x, y).correlation) if len(x) >= 3 else None,
               rule='prereg 2 절 R6: 서술만. 예측 4 (ρ ≤ −0.5) 는 사건 3 개라 순위 상관 값이 −1, −0.5, 0.5, 1 가운데 하나다(서술). 보조 5 사건은 범위 밖.')
    jdump(out, RES + '/r39_r6.json', bk=True); log('R6', [(r['a'], r['b'], round(r['d2n_D'], 4), r['pos_D'], r['pos_tmalign'], round(r['pair_share_per_edge'], 4)) for r in rows], 'ρ', out['rho_pairshare_vs_wtspread']); log('R6_STATS_DONE')


if __name__ == '__main__':
    for part in PARTS_R6:
        log('부분 시작', part)
        try: globals()['part_' + part](); log('부분 끝', part)
        except Exception as e: traceback.print_exc(); log('ERR 부분', part, repr(e)[:300])
