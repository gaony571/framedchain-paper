"""(round41 g5 준비 — 측정 아님; Drive JSON round41_2026-10-21 g5_design 07:33:45) 연산량 표 입력을 /content/ct 로컬에 모은다. 측정 중 Drive 입출력이 없게 한다.
1 도메인 570 (정적 축 BnB 500 쌍): 메도이드(정적 축과 같은 규칙: p3/dom2 의 medoid, 없으면 static/medoid) → x/<d>.npy.
  서열 = p3/dom 요약(info.seq; P3 dom 과 같은 md_atoms) 이 있으면 그것, 없으면 mdCATH h5 pdbProteinAtoms 의 md_atoms (길이가 메도이드 변 수 + 1 과 같을 때만; 아니면 ALA, 기록).
2 트랜스 펩타이드 평면 주형: h5 로 읽은 도메인 가운데 앞 30 의 시작 구조(pdbProteinAtoms 의 N, CA, C, O)에서 틀 R_i (frames_from_backbone) 기준 국소 좌표의 평균 —
  C_i, O_i, N_{i+1} (변 i, Cα_i 기준), N_i (변 i, Cα_i 기준; 첫 잔기용), C_{i+1}·O_{i+1} (변 i, Cα_{i+1} 기준; 끝 잔기용). 트랜스 변(길이 > CIS_MAX)만.
3 메도이드 백본 PDB (pdb/<d>.pdb): Cα = ca_trace(x) (정확), 나머지 = Cα + R_i · 주형, 사슬 A, 잔기 번호 1…L.
4 pairs.json (500 쌍: a, b, n, m, src, 정적 축 BnB 시간), ava200.txt·ava100.txt (default_rng(20261089) 의 570 순열 앞 200·100).
5 앙상블 칸: 줄기 D (derived/r39/D) 에서 MD379-MD320|half 칸이 있고 ne 가 80·150·250 에 가장 가까운 도메인 셋의 320 K·379 K 뒤 절반 첫 실행(N = 100, r38_s0.dom_runs) → ens/<d>.npz.
사용: python r41_compute_prep.py <ROOT>   출력: /content/ct/{x,pdb,ens}/, pairs.json, ava*.txt, prep.json"""
import os, sys
os.environ['CUDA_VISIBLE_DEVICES'] = ''
ROOT = sys.argv[1]
_src = open(ROOT + '/code/r38_s0.py').read(); _src = _src[:_src.index("ALL = ['pick'")]
sys.argv = ['r38_s0', ROOT, 'none']; exec(_src)
from framespace.frames import frames_from_backbone, frames_of, ca_trace, CIS_MAX
T = '/content/ct'
for s in ('x', 'pdb', 'ens'): os.makedirs(f'{T}/{s}', exist_ok=True)
P3, ST = ROOT + '/derived/r38/p3', ROOT + '/derived/r38/static'
A3 = {v: k for k, v in AA.items() if k not in ('HSD', 'HSE', 'HSP')}; CAPS_ = ('ACE', 'NME', 'NHE')


def medoid_of(d):
    f = f'{P3}/dom2/{d}.npz'
    if os.path.exists(f):
        z = np.load(f)
        if 'medoid' in z.files: return z['medoid'].astype(np.float64)
    f = f'{ST}/medoid/{d}.npy'
    return np.load(f).astype(np.float64) if os.path.exists(f) else None


def pdb_text(d):
    h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_{d}.h5', 'r'); t = h[d]['pdbProteinAtoms'][()].decode(); h.close(); return t


def backbone(txt):
    """PDB 글 → 잔기 차례대로 (N, CA, C, O) 가 모두 있는 잔기의 좌표 (L, 4, 3). 캡은 뺀다."""
    res, order = {}, []
    for l in txt.splitlines():
        if not l.startswith(('ATOM', 'HETATM')): continue
        nm, rn, key = l[12:16].strip(), l[17:20].strip(), (l[21], l[22:26], l[26])
        if rn in CAPS_: continue
        if key not in res: res[key] = {}; order.append(key)
        if nm in ('N', 'CA', 'C', 'O') and nm not in res[key]: res[key][nm] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
    return np.array([[res[k][a] for a in ('N', 'CA', 'C', 'O')] for k in order if all(a in res[k] for a in ('N', 'CA', 'C', 'O'))])


def templates(bbs):
    acc = {k: [] for k in ('C', 'O', 'N1', 'N0', 'C1', 'O1')}
    for bb in bbs:
        N, CA, C, O = bb[:, 0], bb[:, 1], bb[:, 2], bb[:, 3]; Rm, l = frames_from_backbone(CA, C, O)
        for i in range(len(Rm)):
            if l[i] <= CIS_MAX: continue
            Rt = Rm[i].T
            acc['C'].append(Rt @ (C[i] - CA[i])); acc['O'].append(Rt @ (O[i] - CA[i])); acc['N1'].append(Rt @ (N[i + 1] - CA[i])); acc['N0'].append(Rt @ (N[i] - CA[i]))
            acc['C1'].append(Rt @ (C[i + 1] - CA[i + 1])); acc['O1'].append(Rt @ (O[i + 1] - CA[i + 1]))
    return {k: np.mean(v, 0) for k, v in acc.items()}, {k: float(np.linalg.norm(np.std(v, 0))) for k, v in acc.items()}, len(acc['C'])


def write_pdb(f, x, seq, tp):
    CA = ca_trace(x); Rm, _ = frames_of(x); L = len(CA); at = []
    for i in range(L):
        Ri = Rm[min(i, L - 2)]
        Np = CA[0] + Rm[0] @ tp['N0'] if i == 0 else CA[i - 1] + Rm[i - 1] @ tp['N1']
        if i < L - 1: Cp, Op = CA[i] + Ri @ tp['C'], CA[i] + Ri @ tp['O']
        else: Cp, Op = CA[i] + Ri @ tp['C1'], CA[i] + Ri @ tp['O1']
        rn = A3.get(seq[i], 'ALA') if seq else 'ALA'
        for nm, p, el in ((' N  ', Np, 'N'), (' CA ', CA[i], 'C'), (' C  ', Cp, 'C'), (' O  ', Op, 'O')):
            at.append('ATOM  ' + f'{len(at) + 1:5d}' + ' ' + nm + ' ' + rn + ' ' + 'A' + f'{i + 1:4d}' + ' ' * 4 + f'{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}' + f'{1.0:6.2f}{0.0:6.2f}' + ' ' * 10 + f'{el:>2s}')
    open(f, 'w').write('\n'.join(at) + '\nTER\nEND\n')


def main():
    t0 = time.time(); b5 = [json.loads(l) for l in open(ST + '/bnb.jsonl')]; doms = sorted({d for r in b5 for d in r['ab']}); info = dict(n_dom=len(doms), seq_mismatch=[], no_medoid=[])
    X, SEQ, BB = {}, {}, []; info['seq_source'] = {'p3_dom': 0, 'h5': 0}
    # 서열: p3/dom 요약(info.seq)에 있으면 그것, 없으면 mdCATH h5 를 하나씩(Drive 가 h5 전체를 받아야 해서 느림 — 병렬은 더 느렸음). 주형은 h5 로 읽은 도메인의 시작 구조.
    for k, d in enumerate(doms):
        x = medoid_of(d)
        if x is None: info['no_medoid'].append(d); continue
        X[d] = x; np.save(f'{T}/x/{d}.npy', x); f = f'{P3}/dom/{d}.npz'; smd = None
        if os.path.exists(f):
            z = np.load(f)
            if 'info' in z.files: smd = json.loads(str(z['info'])).get('seq'); info['seq_source']['p3_dom'] += 1
        if smd is None:
            txt = pdb_text(d); smd, _ = md_atoms(txt); info['seq_source']['h5'] += 1
            if len(BB) < 30:
                bb = backbone(txt)
                if len(bb) > 20: BB.append(bb)
        if len(smd) == len(x) + 1: SEQ[d] = smd
        else: SEQ[d] = None; info['seq_mismatch'].append([d, len(smd), len(x) + 1])
        if k % 100 == 0: log('prep', k, len(doms), round(time.time() - t0), 's')
    tp, tsd, ne_t = templates(BB); info.update(template={k: v.tolist() for k, v in tp.items()}, template_sd=tsd, template_edges=ne_t, template_domains=len(BB))
    for d, x in X.items(): write_pdb(f'{T}/pdb/{d}.pdb', x, SEQ[d], tp)
    pairs = [dict(a=r['ab'][0], b=r['ab'][1], n=r['n'], m=r['m'], src=r['src'], static_bnb_sec=r.get('sec')) for r in b5 if r['ab'][0] in X and r['ab'][1] in X]
    json.dump(pairs, open(f'{T}/pairs.json', 'w')); perm = np.random.default_rng(20261089).permutation(len(doms)); sel = [doms[i] for i in perm if doms[i] in X]
    open(f'{T}/ava200.txt', 'w').write('\n'.join(sel[:200]) + '\n'); open(f'{T}/ava100.txt', 'w').write('\n'.join(sel[:100]) + '\n')
    # 앙상블 칸 도메인 셋
    J = {}
    for f in sorted(glob.glob(ROOT + '/derived/r39/D/*.json')):
        j = json.load(open(f))
        if 'MD379-MD320|half' in j.get('cells', {}): J[j['dom']] = j['ne']
    ens = []
    for tgt in (80, 150, 250):
        d = min(J, key=lambda q: (abs(J[q] - tgt), q)); _runs.clear(); _md.clear()
        Ra = dom_runs(d, '320', 'half', TT)[0]; Rb = dom_runs(d, '379', 'half', TT)[0]
        EA, EB = frames(np.ascontiguousarray(Ra[0][1])), frames(np.ascontiguousarray(Rb[0][1])); np.savez(f'{T}/ens/{d}.npz', EA=EA, EB=EB); ens.append(dict(dom=d, ne=J[d], target=tgt))
    info.update(n_pairs=len(pairs), ava200=len(sel[:200]), ava100=len(sel[:100]), ens=ens, sec=round(time.time() - t0, 1))
    json.dump(info, open(f'{T}/prep.json', 'w'), indent=1); log('PREP_DONE', {k: v for k, v in info.items() if k not in ('template', 'template_sd')})


if __name__ == '__main__':
    main()
