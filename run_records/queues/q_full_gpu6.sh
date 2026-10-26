# (2026-10-20) E7 30 → 0.2 다음 단계(bias2) → P2 D 판(67) → P3 π̄ (DOM_DONE) → B6 앙상블 → 민감도 메도이드 정렬 (DOM2_DONE) → 4 단계 AlphaFlow → 5 → 조건부. E7 30 python 은 끊지 않음.
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
while kill -0 92569 2>/dev/null; do sleep 30; done
echo E7_30_DONE_Q >> $L/r38_e7_30.log
cd /content
/usr/bin/python3 $C/r38_s0_bias2.py $R >> $L/r38_s0_bias2.log 2>&1
echo BIAS2_DONE_Q >> $L/r38_s0_bias2.log
/usr/bin/python3 $C/r38_p2_D.py $R $R/derived/r38/q3_p2_67_ordered.txt >> $L/r38_p2_D.log 2>&1
/usr/bin/python3 $C/r38_p2_aggregate_D.py $R $R/derived/r38/q3_p2_67_ordered.txt p2_67_D full 0.75 >> $L/r38_p2_D.log 2>&1
echo P2_D_DONE_Q >> $L/r38_p2_D.log
until grep -qE 'DOM_DONE|P3_DONE dom' $L/r38_p3.log 2>/dev/null; do sleep 60; done
/usr/bin/python3 $C/r38_p3.py $R pibar >> $L/r38_p3.log 2>&1
echo PIBAR_DONE_Q >> $L/r38_p3.log
/usr/bin/python3 $C/r38_b6.py $R ens >> $L/r38_b6.log 2>&1
echo B6_ENS_DONE_Q >> $L/r38_b6.log
until grep -qE 'DOM2_DONE|P3_SENS_PARTS_DONE dom2' $L/r38_p3_sens.log 2>/dev/null; do sleep 60; done
/usr/bin/python3 $C/r38_p3_sens.py $R medalign >> $L/r38_p3_sens.log 2>&1
echo MEDALIGN_DONE_Q >> $L/r38_p3_sens.log
until [ -f /content/AFENV_READY ]; do sleep 60; done
/usr/bin/python3 $C/r38_alphaflow_gen.py $R $R/derived/r38/q3_p2_67_ordered.txt >> $L/r38_alphaflow_gen.log 2>&1
echo STEP4_DONE >> $L/r38_alphaflow_gen.log
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q5a_450_17.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q5b_p90_4.txt 320 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/rest18.txt 320 50,25,10 >> $L/r38_asamt_gen.log 2>&1
echo STEP5_DONE >> $L/r38_asamt_gen.log
# 조건부(시간 남으면; 트랙 B 뒤)
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q_cond_450_50.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q3_p2_67_ordered.txt 348,379,413 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_alphaflow_gen.py $R $R/derived/r38/rest18.txt >> $L/r38_alphaflow_gen.log 2>&1
echo FULL_GPU_DONE >> $L/r38_alphaflow_gen.log
