
"""(prereg round36 7b 0) 바름: AlphaFlow 원 분석 코드(scripts/analyze_ensembles.py, print_analysis.py, 커밋 0408d7c)를 그대로
공개 표본 넷에 돌려 논문 표와 견줌. 기준 = ATLAS 세 반복(분석 zip 의 R1–R3 xtc 를 _prod_R{k}_fit.xtc 로 이름만 바꿈).
derived/r38/af_bar/<set>_out.pkl, results/r38_af_bar.txt"""
import os, sys, glob, zipfile, shutil, subprocess, time
R=sys.argv[1]; AF=R+'/code/ext/alphaflow'; W='/content/r38_bar'; AD=W+'/atlas'; os.makedirs(AD,exist_ok=True)
OUT=R+'/derived/r38/af_bar'; os.makedirs(OUT,exist_ok=True)
def backup(p):
    b=R+'/results/_backup/'+time.strftime('%Y-%m-%d'); os.makedirs(b,exist_ok=True); shutil.copy(p,b+'/'+time.strftime('%H%M%S_')+os.path.basename(p))
for zp in sorted(glob.glob(R+'/raw/atlas_af82/*_analysis.zip')):
    n=os.path.basename(zp)[:-13]; d=f'{AD}/{n}'
    if os.path.exists(f'{d}/{n}_prod_R3_fit.xtc'): continue
    os.makedirs(d,exist_ok=True); z=zipfile.ZipFile(zp)
    z.extract(f'{n}.pdb',d)
    for k in (1,2,3): open(f'{d}/{n}_prod_R{k}_fit.xtc','wb').write(z.read(f'{n}_R{k}.xtc'))
names=sorted(os.listdir(AD)); print('atlas',len(names),flush=True)
SETS=['alphaflow_md_base_202402','alphaflow_md_templates_base_202402','esmflow_md_base_202402','esmflow_md_templates_base_202402']
# 2026-10-19 사용자 결정: 바름은 AlphaFlow-MD 한 세트(prereg 7b 0). argv[2] = 쉼표로 고른 세트, argv[3] = 작업자 수(결과에는 영향 없음, 표적마다 Pool)
if len(sys.argv)>2 and sys.argv[2]: SETS=sys.argv[2].split(',')
NW=int(sys.argv[3]) if len(sys.argv)>3 else 1
done=[]
for s in SETS:
    pk=f'{OUT}/{s}_out.pkl'
    if os.path.exists(pk): done.append(pk); continue
    sd=f'{W}/{s}'; shutil.rmtree(sd,ignore_errors=True); os.makedirs(sd)
    z=zipfile.ZipFile(f'{R}/raw/models/alphaflow/samples/{s}.zip')
    for m in z.namelist():
        b=os.path.basename(m)
        if b.endswith('.pdb') and b[:-4] in names: open(f'{sd}/{b}','wb').write(z.read(m))
    print(s,'pdbs',len(glob.glob(sd+'/*.pdb')),flush=True)
    r=subprocess.run(f"cd {AF} && OMP_NUM_THREADS=1 python -m scripts.analyze_ensembles --atlas_dir {AD} --pdbdir {sd} --num_workers {NW}",shell=True)
    if os.path.exists(f'{sd}/out.pkl'): shutil.copy(f'{sd}/out.pkl',pk); backup(pk); done.append(pk)
    print(s,'rc',r.returncode,flush=True); shutil.rmtree(sd,ignore_errors=True)
txt=subprocess.run(f"cd {AF} && python -m scripts.print_analysis "+' '.join(done),shell=True,capture_output=True,text=True)
P=R+'/results/r38_af_bar.txt'; open(P,'w').write(txt.stdout+'\n'+txt.stderr[-2000:]); backup(P); print(txt.stdout,flush=True)
open(R+'/results/logs/r38_af_bar.done','w').write(time.strftime('%F %T'))
