# (2026-10-20 UTC) 확인 실험 P1′ (prereg 7h) CPU 대기열: 내려받기 끝 → 검증 통과 → mean → align(작업자 6) → spread → fold → stats
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
export CUDA_VISIBLE_DEVICES=
until grep -q P1N_FETCH_DONE $L/r38_p1n_fetch.log 2>/dev/null; do sleep 120; done
until grep -q VERIFY_P1N_DONE $L/r38_verify_p1n.log 2>/dev/null; do sleep 60; done
/usr/bin/python3 -c "import json,sys; j=json.load(open('$R/results/r38_verify_p1n.json')); sys.exit(0 if j['n']>=5 and j['n_ok']==j['n'] else 1)" || { echo P1N_VERIFY_FAIL >> $L/r38_p1n.log; exit 1; }
nice -n 5 /usr/bin/python3 $C/r38_p1n.py $R mean 4 >> $L/r38_p1n.log 2>&1
nice -n 5 /usr/bin/python3 $C/r38_p1n.py $R align 6 >> $L/r38_p1n.log 2>&1
nice -n 5 /usr/bin/python3 $C/r38_p1n.py $R spread 4 >> $L/r38_p1n.log 2>&1
nice -n 5 /usr/bin/python3 $C/r38_p1n.py $R fold 4 >> $L/r38_p1n.log 2>&1
/usr/bin/python3 $C/r38_p1n.py $R stats 1 >> $L/r38_p1n.log 2>&1
echo P1N_ALL_DONE >> $L/r38_p1n.log
