"""(prereg round36 7c) AlphaFlow-MD 증류판 생성. 원 코드(커밋 0408d7c) predict.py 를 README 그대로 부른다:
--mode alphafold --samples 250 --noisy_first --no_diffusion (증류판), 나머지 기본값(steps 10, tmax 1.0).
원 코드에 seed 인자가 없어 도메인마다 torch·numpy·random seed 만 고정한다(zlib.crc32). 한 도메인씩, 출력이 있으면 건너뜀.
사용: python r38_alphaflow_gen.py <ROOT> <도메인 목록 파일 또는 쉼표 목록>
출력: derived/r38/alphaflow/<dom>.pdb (250 모델)"""
import os, sys, zlib, time, shutil, subprocess
ROOT=sys.argv[1]; arg=sys.argv[2]
doms=[x.strip() for x in (open(arg).read().split() if os.path.exists(arg) else arg.split(',')) if x.strip()]
AF=ROOT+'/code/ext/alphaflow'; OUT=ROOT+'/derived/r38/alphaflow'; os.makedirs(OUT,exist_ok=True); W='/content/af_work'; os.makedirs(W,exist_ok=True)
PY='/content/afenv/bin/python'
for d in doms:
    if os.path.exists(f'{OUT}/{d}.pdb'): continue
    seed=zlib.crc32(f'r38_alphaflow|{d}'.encode()) % (2**31)
    # predict.py 는 --pdb_id 를 줘도 CSV 전체 줄의 특징을 읽고 건너뛴다 → 그 도메인 한 줄만 담은 CSV 를 넘김(같은 입력)
    rows=[l for l in open(ROOT+'/derived/r38/alphaflow_input.csv').read().splitlines()]
    one=f'{W}/{d}.csv'; open(one,'w').write(rows[0]+'\n'+[l for l in rows[1:] if l.split(',')[0]==d][0]+'\n')
    args=['predict.py','--mode','alphafold','--input_csv',one,'--msa_dir',ROOT+'/derived/r38/msa',
          '--weights',ROOT+'/raw/models/alphaflow/params/alphaflow_md_distilled_202402.pt','--samples','250','--outpdb',W,
          '--noisy_first','--no_diffusion','--pdb_id',d,'--runtime_json',f'{W}/{d}.runtime.json']
    # torch 1.12(README 판)의 torch.load 에는 weights_only 인자가 없다. 이 커밋의 predict.py 는 torch 2.6 용으로 weights_only=False 를 넘긴다.
    # 1.12 의 load 는 원래 전부 읽으므로(= weights_only=False) 인자만 떼어 넘긴다. 계산은 같다.
    wrap=(f"import torch, numpy, random, sys, runpy; _L=torch.load; torch.load=lambda *a, **k: _L(*a, **{{x:y for x,y in k.items() if x!='weights_only'}}); torch.manual_seed({seed}); numpy.random.seed({seed}); random.seed({seed}); "
          f"sys.argv={args!r}; runpy.run_path('predict.py', run_name='__main__')")
    t0=time.time(); r=subprocess.run([PY,'-c',wrap],cwd=AF,capture_output=True,text=True)
    ok=os.path.exists(f'{W}/{d}.pdb')
    if ok:
        shutil.copy(f'{W}/{d}.pdb', f'{OUT}/{d}.pdb')
        if os.path.exists(f'{W}/{d}.runtime.json'): shutil.copy(f'{W}/{d}.runtime.json', f'{OUT}/{d}.runtime.json')
    print(time.strftime('%F %T'), d, 'seed', seed, 'ok' if ok else 'FAIL', round(time.time()-t0,1), 's', '' if ok else (r.stdout[-400:]+r.stderr[-1500:]), flush=True)
print('ALPHAFLOW_GEN_DONE', flush=True)
