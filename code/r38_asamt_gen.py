"""(prereg round36 7c) aSAMt 생성. 원 코드 sam2(커밋 ce076ed) generate_ensemble 을 기본값 그대로(250 표본, 100 걸음, 최소화 켬) 부른다.
원 코드에 seed 인자가 없어 (도메인, 온도)마다 torch·numpy·random seed 만 고정한다(zlib.crc32). 이어 돌리기: 출력이 있으면 건너뜀.
사용: python r38_asamt_gen.py <ROOT> <도메인 목록 파일 또는 쉼표 목록> [온도들=320,348,379,413,450] [배치 크기들=50,25,10]
출력: derived/r38/asamt/<dom>.<T>.{traj.dcd,top.pdb,time.txt} (최소화한 앙상블), results/logs/r38_asamt_gen.log"""
import os, sys, zlib, time, zipfile, shutil, random, subprocess, glob
ROOT=sys.argv[1]; arg=sys.argv[2]; TS=(sys.argv[3] if len(sys.argv)>3 and sys.argv[3] else '320,348,379,413,450').split(',')
# 배치 크기(2026-10-19 사용자 결정: 키움). 표본 분포는 같고 난수 흐름만 다르다. GPU 메모리가 모자라면 줄여 다시.
BS=[int(x) for x in (sys.argv[4] if len(sys.argv)>4 else '50,25,10').split(',')]
doms=[x.strip() for x in (open(arg).read().split() if os.path.exists(arg) else arg.split(',')) if x.strip()]
OUT=ROOT+'/derived/r38/asamt'; os.makedirs(OUT,exist_ok=True); W='/content/asamt_work'; os.makedirs(W,exist_ok=True)
Z=zipfile.ZipFile(ROOT+'/raw/models/asam/input_mdcath.zip')
PY='/content/samenv/bin/python'  # 설정: 내장 'mdcath'(= src/sam/config/mdcath_model.yaml). -c(--config_path)는 원 코드가 str→Path 변환을 안 해 실패
for d in doms:
    for T in TS:
        tag=f'{d}.{T}'
        if os.path.exists(f'{OUT}/{tag}.traj.dcd') and os.path.exists(f'{OUT}/{tag}.top.pdb'): continue
        for f in glob.glob(f'{W}/{tag}.*'): os.remove(f)
        inp=f'{W}/{tag}.in.pdb'; open(inp,'wb').write(Z.read(f'input_mdcath/{tag}.pdb'))
        seed=zlib.crc32(f'r38_asamt|{tag}'.encode()) % (2**31)
        # 원 코드 215 줄의 오타 _LOGGER.infox 를 info 로 이어 줌(파일은 그대로)
        wrap=(f"import torch, numpy, random, sys, runpy, logging; logging.Logger.infox=logging.Logger.info; torch.manual_seed({seed}); numpy.random.seed({seed}); random.seed({seed}); "
              f"sys.argv=['generate_ensemble.py','--config','mdcath','-i','{inp}','-o','{W}/{tag}','-T','{T}','-d','cuda','--data_dir','/content/sam_data','--no_download','--time']; "
              f"runpy.run_path('/content/sam2_src/src/sam/scripts/generate_ensemble.py', run_name='__main__')")
        t0=time.time()
        for b in BS:
            for f in glob.glob(f'{W}/{tag}.*'):
                if not f.endswith('.in.pdb'): os.remove(f)
            r=subprocess.run([PY,'-c',wrap.replace("'--time']", f"'--time','-b','{b}']")],capture_output=True,text=True)
            if 'out of memory' not in (r.stderr+r.stdout).lower(): break
        ok=os.path.exists(f'{W}/{tag}.traj.dcd') and os.path.exists(f'{W}/{tag}.top.pdb') and os.path.exists(f'{W}/{tag}.time.txt')
        if ok:
            for ext in ('traj.dcd','top.pdb','time.txt'): shutil.copy(f'{W}/{tag}.{ext}', f'{OUT}/{tag}.{ext}')
            for f in os.listdir(W):
                if f.startswith(tag+'.') and f.endswith('.txt'): shutil.copy(f'{W}/{f}', f'{OUT}/{f}')
        print(time.strftime('%F %T'), tag, 'seed', seed, 'batch', b, 'ok' if ok else 'FAIL', round(time.time()-t0,1), 's', '' if ok else (r.stdout[-400:]+r.stderr[-1200:]), flush=True)
print('ASAMT_GEN_DONE', flush=True)
