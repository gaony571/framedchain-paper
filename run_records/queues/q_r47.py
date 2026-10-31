"""r47_ties 실행 순서: 부분 일곱(ties, p3·nl·anm × argsort·average)을 동시에 띄우고, 모두 파일을 남기면 combine. 로그 results/logs/r47_ties.log."""
import subprocess, sys, os, time
R = sys.argv[1]; PY = sys.executable; S = f'{R}/code/r47_ties.py'; LG = f'{R}/results/logs'; os.makedirs(LG, exist_ok=True)
jobs = [('ties', 'none'), ('p3', 'argsort'), ('p3', 'average'), ('nl', 'argsort'), ('nl', 'average'), ('anm', 'argsort'), ('anm', 'average')]
fh = open(f'{LG}/r47_ties.log', 'a'); fh.write(time.strftime('%H:%M:%S') + ' R47 시작\n'); fh.flush()
ps = [subprocess.Popen([PY, S, R, p, m], stdout=fh, stderr=subprocess.STDOUT) for p, m in jobs]
rc = [p.wait() for p in ps]
fh.write(time.strftime('%H:%M:%S') + f' 부분 끝 {rc}\n'); fh.flush()
need = [f'{R}/derived/r47/{p}_{m}.json' for p, m in jobs]
if all(os.path.exists(f) for f in need): subprocess.run([PY, S, R, 'combine'], stdout=fh, stderr=subprocess.STDOUT)
else: fh.write('R47 STOP: 부분 파일 없음 ' + str([os.path.basename(f) for f in need if not os.path.exists(f)]) + '\n')
fh.write(time.strftime('%H:%M:%S') + ' R47_ALL_DONE\n'); fh.close()
