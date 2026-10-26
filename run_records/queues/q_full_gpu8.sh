# (2026-10-20 UTC) q_full_gpu7 대체: AlphaFlow 4 단계(PID 10792) 끝 → C pibar(1c; C_DOM_FINAL 이면 바로, 아니면 5 단계 뒤) → aSAMt 5 단계 → 조건부
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
while kill -0 10792 2>/dev/null; do sleep 30; done
echo STEP4_DONE >> $L/r38_alphaflow_gen.log
run_c() { /usr/bin/python3 $C/r38_c.py $R pibar 1 >> $L/r38_c.log 2>&1; echo C_PIBAR_Q_DONE >> $L/r38_c.log; }
CD=0; grep -q C_DOM_FINAL $L/r38_c.log 2>/dev/null && { run_c; CD=1; }
until [ -f /content/SAMENV_READY ]; do sleep 30; done
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q5a_450_17.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q5b_p90_4.txt 320 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/rest18.txt 320 50,25,10 >> $L/r38_asamt_gen.log 2>&1
echo STEP5_DONE >> $L/r38_asamt_gen.log
if [ $CD = 0 ]; then
  until grep -qE 'C_DOM_FINAL|C_VERIFY_FAIL' $L/r38_c.log 2>/dev/null; do sleep 60; done
  grep -q C_DOM_FINAL $L/r38_c.log && run_c
fi
# 조건부(시간 남으면; 트랙 B 뒤)
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q_cond_450_50.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q3_p2_67_ordered.txt 348,379,413 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_alphaflow_gen.py $R $R/derived/r38/rest18.txt >> $L/r38_alphaflow_gen.log 2>&1
echo FULL_GPU_DONE >> $L/r38_alphaflow_gen.log
