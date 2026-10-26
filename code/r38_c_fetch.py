"""(확인 실험 C) 새 상과 도메인(derived/r38/c_new_sf_domains.json 의 to_download 883)을 HuggingFace mdCATH 에서 HTTP 부분 읽기로 받는다(사용자 허락 2026-10-20 UTC).
도메인마다 pdbProteinAtoms + 320 K 반복마다 첫 프레임과 뒤 절반 프레임의 골격 원자(r37_ff.md_atoms 의 원자 셋)만 → raw/mdcath_c/<d>.npz
(smd 서열, first (L, 3, 3) = 반복 0 첫 프레임, run<r> (n_half, L, 3, 3) float32, nframes 반복별 전체 프레임 수). 이미 있으면 건너뜀. 프로세스 8.
사용: python r38_c_fetch.py <ROOT> [작업자=8]   로그: 진행·오류·받은 MB"""
import os, sys, json, time, zlib, traceback
ROOT = sys.argv[1]; NW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
import numpy as np, h5py, fsspec
src = open(ROOT + '/code/r37_ff.py').read(); exec(src[src.index('AA = dict('):src.index('def stats')])   # md_atoms
OUTC = ROOT + '/raw/mdcath_c'; os.makedirs(OUTC, exist_ok=True)
URL = 'https://huggingface.co/datasets/compsciencelab/mdCATH/resolve/main/data/mdcath_dataset_{}.h5'


def one(d):
    fo_ = f'{OUTC}/{d}.npz'
    if os.path.exists(fo_): return d, 'skip', 0.0
    for att in range(3):
        try:
            t0 = time.time(); fo = fsspec.open(URL.format(d), 'rb', block_size=4 * 2 ** 20, cache_type='readahead').open(); h = h5py.File(fo, 'r'); g = h[d]
            txt = g['pdbProteinAtoms'][()].decode(); smd, imd = md_atoms(txt); L = len(smd); out = dict(smd=np.array(list(smd)), nframes=[])
            reps = sorted(g['320'].keys(), key=int)
            for k, r in enumerate(reps):
                ds = g['320'][r]['coords']; n = ds.shape[0]; out['nframes'].append(n)
                if k == 0: out['first'] = ds[0][imd.ravel()].reshape(L, 3, 3).astype(np.float32)
                out[f'run{r}'] = ds[n // 2:][:, imd.ravel()].reshape(-1, L, 3, 3).astype(np.float32)
            h.close(); fo.close(); out['nframes'] = np.array(out['nframes']); out['reps'] = np.array([int(r) for r in reps])
            tmp = f'{OUTC}/.{d}.tmp.npz'; np.savez(tmp, **out); os.replace(tmp, fo_); return d, 'ok', time.time() - t0
        except Exception as e:
            err = repr(e)[:150]; time.sleep(5 * (att + 1))
    return d, 'ERR ' + err, 0.0


if __name__ == '__main__':
    from multiprocessing import Pool
    todo = json.load(open(ROOT + '/derived/r38/c_new_sf_domains.json'))['to_download']; t0 = time.time(); n_ok = n_err = 0
    print(time.strftime('%T'), '받을 도메인', len(todo), '작업자', NW, flush=True)
    with Pool(NW) as p:
        for k, (d, s, dt) in enumerate(p.imap_unordered(one, todo)):
            n_ok += s == 'ok'; n_err += s.startswith('ERR')
            if s.startswith('ERR') or k % 25 == 0: print(time.strftime('%T'), k + 1, '/', len(todo), d, s, round(dt, 1), 's', '경과', round(time.time() - t0), 's', flush=True)
    print(time.strftime('%T'), 'CFETCH_DONE ok', n_ok, 'err', n_err, flush=True)
