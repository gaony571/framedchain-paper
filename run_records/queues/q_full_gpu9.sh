# (2026-10-20 UTC) q_full_gpu8 대체: AlphaFlow 4 단계(PID 10792) 끝 → C pibar(1c; C_DOM_FINAL 이면 바로) → aSAMt 5 단계 → (C pibar 남았으면) → 7g: E7 100 → d_ff D → 정적 축 → 조건부
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
while kill -0 10792 2>/dev/null; do sleep 30; done
grep -q STEP4_DONE $L/r38_alphaflow_gen.log || echo STEP4_DONE >> $L/r38_alphaflow_gen.log
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
# 7g (a) 온도 눈금 D 판 100 도메인
until [ -f $R/derived/r38/s0/e7_100/pick.json ]; do sleep 60; done
/usr/bin/python3 $C/r38_e7_100.py $R run,sum >> $L/r38_e7_100.log 2>&1
# 7g (b) R3′ d_ff D 판
/usr/bin/python3 $C/r38_dff_D.py $R >> $L/r38_dff_D.log 2>&1
# 7g 1 절 정적 축
until grep -q STATIC_DOM_DONE $L/r38_static.log 2>/dev/null; do sleep 60; done
/usr/bin/python3 $C/r38_static.py $R align,bnb,stats 1 >> $L/r38_static.log 2>&1
echo G7_DONE >> $L/r38_static.log
# 조건부(시간 남으면; 트랙 B 뒤)
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q_cond_450_50.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q3_p2_67_ordered.txt 348,379,413 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_alphaflow_gen.py $R $R/derived/r38/rest18.txt >> $L/r38_alphaflow_gen.log 2>&1
echo FULL_GPU_DONE >> $L/r38_alphaflow_gen.log
