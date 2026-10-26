# (2026-10-20 UTC) 확인 실험 C (prereg 1c) CPU 대기열: biopython → pairs → 검증 → H-풀 ANM → dom(받는 동안 되풀이) → C ANM → [GPU pibar 기다림] → pair·anmpair·stats
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
export CUDA_VISIBLE_DEVICES=
/usr/bin/python3 -m pip install -q biopython >> $L/r38_c.log 2>&1
/usr/bin/python3 $C/r38_c.py $R pairs 1 >> $L/r38_c.log 2>&1
/usr/bin/python3 $C/r38_verify_c.py $R >> $L/r38_verify_c.log 2>&1
/usr/bin/python3 -c "import json,sys; j=json.load(open('$R/results/r38_verify_c.json')); sys.exit(0 if j['n']>=8 and j['n_ok']==j['n'] else 1)" || { echo C_VERIFY_FAIL >> $L/r38_c.log; exit 1; }
echo C_VERIFY_OK >> $L/r38_c.log
nice -n 5 /usr/bin/python3 $C/r38_c.py $R anm 2 >> $L/r38_c.log 2>&1
while true; do
  F=0; grep -q CFETCH_DONE $L/r38_c_fetch.log && F=1; pgrep -f r38_c_fetch.py > /dev/null || F=1
  nice -n 5 /usr/bin/python3 $C/r38_c.py $R dom 3 >> $L/r38_c.log 2>&1
  [ $F = 1 ] && break
  sleep 300
done
echo C_DOM_FINAL >> $L/r38_c.log
nice -n 5 /usr/bin/python3 $C/r38_c.py $R anm 3 >> $L/r38_c.log 2>&1
echo C_CPU1_DONE >> $L/r38_c.log
until grep -q C_PIBAR_Q_DONE $L/r38_c.log; do sleep 120; done
/usr/bin/python3 $C/r38_c.py $R pair,anmpair,stats 1 >> $L/r38_c.log 2>&1
echo C_ALL_DONE >> $L/r38_c.log
