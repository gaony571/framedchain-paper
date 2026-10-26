"""(확인 실험 C 준비) mdCATH h5 의 내부 구조와 HTTP 부분 읽기(320 K 만) 시험. 백그라운드로만 돌린다(커널에서 h5 를 열지 않음).
A: Drive 의 h5 하나(1kskA01) — 데이터셋 모양·덩어리·압축·저장 크기, 320 K 몫.
B: 새 상과 도메인 둘(가장 작은 것, 300 잔기 근처)을 HuggingFace 에서 fsspec 으로 열어 pdbProteinAtoms + 320 K 반복마다 첫 프레임·뒤 절반만 읽음 → 받은 바이트·시간.
출력: results/r38_c_probe.json.   사용: python r38_c_probe.py <ROOT>"""
import sys, os, json, time, itertools, collections
import numpy as np, h5py, fsspec
ROOT = sys.argv[1]; out = dict(written=time.strftime('%F %T'))
h = h5py.File(f'{ROOT}/raw/mdcath/mdcath_dataset_1kskA01.h5', 'r'); g = h['1kskA01']; info = {}
def vis(name, obj):
    if isinstance(obj, h5py.Dataset): info[name] = dict(shape=list(obj.shape), dtype=str(obj.dtype), chunks=list(obj.chunks) if obj.chunks else None, comp=obj.compression, bytes=int(obj.id.get_storage_size()))
g.visititems(vis); h.close()
tot = sum(v['bytes'] for v in info.values()); t320 = sum(v['bytes'] for k, v in info.items() if k.startswith('320/'))
out['A'] = dict(n_datasets=len(info), total_MB=tot / 1e6, t320_MB=t320 / 1e6, frac320=t320 / tot, coords_320_0=info.get('320/0/coords'), keys_320_0=sorted(k.split('/')[-1] for k in info if k.startswith('320/0/')),
                top=[k for k in info if '/' not in k])
print('A', json.dumps(out['A'], ensure_ascii=False)[:800], flush=True)
# 새 상과 도메인 목록
alld = json.load(open(ROOT + '/raw/mdcath_all_domains_listing.json'))['domains']
cath = {}
for l in open(ROOT + '/raw/cath/cath-domain-list.txt'):
    if l.startswith('#') or not l.strip(): continue
    f = l.split(); cath[f[0]] = ('.'.join(f[1:5]), int(f[10]))
import pandas as pd
hp = pd.read_csv(ROOT + '/results/pilot_mdcath_homolog_pairs.tsv', sep='\t'); hsf = set(hp.H); hdom = set(hp.a) | set(hp.b)
by = collections.defaultdict(list)
for d in alld:
    if d in cath: by[cath[d][0]].append(d)
need = sorted({d for H, ds in by.items() if H not in hsf and len([x for x in ds if x not in hdom]) >= 2 for d in ds if d not in hdom})
h5s = {os.path.basename(x)[len('mdcath_dataset_'):-3] for x in os.listdir(ROOT + '/raw/mdcath') if x.endswith('.h5')}
need_dl = [d for d in need if d not in h5s]; out['need'] = dict(n_domains=len(need), n_to_download=len(need_dl), residues=int(sum(cath[d][1] for d in need_dl)))
json.dump(dict(domains=need, to_download=need_dl), open(ROOT + '/derived/r38/c_new_sf_domains.json', 'w'))
print('need', out['need'], flush=True)
srt = sorted(need_dl, key=lambda d: cath[d][1]); picks = [srt[0], min(srt, key=lambda d: abs(cath[d][1] - 300))]
URL = 'https://huggingface.co/datasets/compsciencelab/mdCATH/resolve/main/data/mdcath_dataset_{}.h5'


class Counter:
    def __init__(self, f): self.f = f; self.n = 0
    def read(self, k=-1): b = self.f.read(k); self.n += len(b); return b
    def readinto(self, b): n = self.f.readinto(b); self.n += n; return n
    def seek(self, *a): return self.f.seek(*a)
    def tell(self): return self.f.tell()
    def seekable(self): return True
    def readable(self): return True


out['B'] = []
for d in picks:
    t0 = time.time(); fsz = fsspec.filesystem('https').size(URL.format(d))
    fo = fsspec.open(URL.format(d), 'rb', block_size=4 * 2 ** 20, cache_type='readahead').open(); cf = Counter(fo)
    h = h5py.File(cf, 'r'); g = h[d]; pdb = g['pdbProteinAtoms'][()]; reps = sorted(g['320'].keys(), key=int); shapes = []
    for r in reps:
        ds = g['320'][r]['coords']; n = ds.shape[0]; x0 = ds[0] if r == reps[0] else None; xh = ds[n // 2:]; shapes.append([int(n), list(ds.shape[1:]), list(ds.chunks) if ds.chunks else None])
    h.close(); dt = time.time() - t0
    out['B'].append(dict(dom=d, length=cath[d][1], file_MB=fsz / 1e6, read_MB=cf.n / 1e6, frac=cf.n / fsz, sec=dt, MBps=cf.n / 1e6 / dt, runs=shapes, pdb_bytes=len(pdb)))
    print('B', out['B'][-1], flush=True)
fr = np.mean([b['frac'] for b in out['B']]); out['estimate'] = dict(full_TB=out['need']['residues'] * 4.96e6 / 1e12, partial_GB=out['need']['residues'] * 4.96e6 * fr / 1e9, frac_read=fr,
                                                                   hours_at_MBps=out['need']['residues'] * 4.96e6 * fr / 1e6 / np.mean([b['MBps'] for b in out['B']]) / 3600)
json.dump(out, open(ROOT + '/results/r38_c_probe.json', 'w'), ensure_ascii=False, indent=1, default=float); print('estimate', out['estimate'], flush=True); print('CPROBE_DONE')
