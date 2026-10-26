"""(prereg 7h, 확인 실험 P1′) 새 쌍의 도메인 가운데 Drive 에 h5 가 없는 311 개(derived/r38/p1n/pairs.json 의 to_download)를 HuggingFace mdCATH 에서 HTTP 부분 읽기로 받는다
(사용자 허락 2026-10-20 UTC "다 넣어서 가"). 도메인마다 pdbProteinAtoms + 320·348·379·413 K 반복마다 뒤 절반(h<T>_<r>)과 50–100 ns(a<T>_<r>, 프레임 50:100,
반복 프레임 ≥ 100 일 때만)의 골격 원자(r37_ff.md_atoms 의 원자 셋), 320 K 반복 0 첫 프레임(first), 반복별 전체 프레임 수(nframes<T>), 반복 번호(reps<T>)
→ raw/mdcath_p1n/<d>.npz. r38_p1.build_cells 가 h5 에서 쓰는 프레임과 같은 프레임이다. 이미 있으면 건너뜀. 프로세스 8.
사용: python r38_p1n_fetch.py <ROOT> [작업자=8]"""
import os, sys, json, time, traceback
ROOT = sys.argv[1]; NW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
import numpy as np, h5py, fsspec
src = open(ROOT + '/code/r37_ff.py').read(); exec(src[src.index('AA = dict('):src.index('def stats')]) # md_atoms
OUTN = ROOT + '/raw/mdcath_p1n'; os.makedirs(OUTN, exist_ok=True)
URL = 'https://huggingface.co/datasets/compsciencelab/mdCATH/resolve/main/data/mdcath_dataset_{}.h5'
TT4 = ('320', '348', '379', '413')


def one(d):
    fo_ = f'{OUTN}/{d}.npz'
    if os.path.exists(fo_): return d, 'skip', 0.0
    for att in range(3):
        try:
            t0 = time.time(); fo = fsspec.open(URL.format(d), 'rb', block_size=4 * 2 ** 20, cache_type='readahead').open(); h = h5py.File(fo, 'r'); g = h[d]
            smd, imd = md_atoms(g['pdbProteinAtoms'][()].decode()); L = len(smd); out = dict(smd=np.array(list(smd)))
            for T in TT4:
                if T not in g: continue
                reps = sorted(g[T].keys(), key=int); nf = []
                for k, r in enumerate(reps):
                    ds = g[T][r]['coords']; n = ds.shape[0]; nf.append(n)
                    if T == '320' and k == 0: out['first'] = ds[0][imd.ravel()].reshape(L, 3, 3).astype(np.float32)
                    out[f'h{T}_{r}'] = ds[n // 2:][:, imd.ravel()].reshape(-1, L, 3, 3).astype(np.float32)
                    if n >= 100: out[f'a{T}_{r}'] = ds[50:100][:, imd.ravel()].reshape(-1, L, 3, 3).astype(np.float32)
                out[f'nframes{T}'] = np.array(nf); out[f'reps{T}'] = np.array([int(r) for r in reps])
            h.close(); fo.close()
            tmp = f'{OUTN}/.{d}.tmp.npz'; np.savez(tmp, **out); os.replace(tmp, fo_); return d, 'ok', time.time() - t0
        except Exception as e:
            err = repr(e)[:150]; time.sleep(5 * (att + 1))
    return d, 'ERR ' + err, 0.0


if __name__ == '__main__':
    from multiprocessing import Pool
    todo = json.load(open(ROOT + '/derived/r38/p1n/pairs.json'))['to_download']; t0 = time.time(); n_ok = n_err = 0
    print(time.strftime('%T'), '받을 도메인', len(todo), '작업자', NW, flush=True)
    with Pool(NW) as p:
        for k, (d, s, dt) in enumerate(p.imap_unordered(one, todo)):
            n_ok += s == 'ok'; n_err += s.startswith('ERR')
            if s.startswith('ERR') or k % 25 == 0: print(time.strftime('%T'), k + 1, '/', len(todo), d, s, round(dt, 1), 's', '경과', round(time.time() - t0), 's', flush=True)
    print(time.strftime('%T'), 'P1N_FETCH_DONE ok', n_ok, 'err', n_err, flush=True)
