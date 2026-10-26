# (2026-10-20 UTC) q_full_gpu9 대체: C π̄(PID 81924) 끝 → C_PIBAR_Q_DONE → aSAMt 5 단계 → 7g E7 100 → d_ff D → AlphaFlow D 칸(7c 2차 2·3 D 판)·집계 → 정적 축 → 조건부
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
while kill -0 81924 2>/dev/null; do sleep 30; done
echo C_PIBAR_Q_DONE >> $L/r38_c.log
until [ -f /content/SAMENV_READY ]; do sleep 30; done
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q5a_450_17.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q5b_p90_4.txt 320 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/rest18.txt 320 50,25,10 >> $L/r38_asamt_gen.log 2>&1
echo STEP5_DONE >> $L/r38_asamt_gen.log
# 7g (a) 온도 눈금 D 판 100 도메인
until [ -f $R/derived/r38/s0/e7_100/pick.json ]; do sleep 60; done
/usr/bin/python3 $C/r38_e7_100.py $R run,sum >> $L/r38_e7_100.log 2>&1
# 7g (b) R3′ d_ff D 판
/usr/bin/python3 $C/r38_dff_D.py $R >> $L/r38_dff_D.log 2>&1
# 7c 2차 2·3 의 D 판: AlphaFlow D 칸 → 동결 집계(입력 compare_D2·d_ff D)
/usr/bin/python3 $C/r38_p2_D_af.py $R $R/derived/r38/q3_p2_67_ordered.txt >> $L/r38_p2_D_af.log 2>&1
for t in "q3_p2_67_ordered p2_67_D2" "p2_sf23 p2_sf23_D2" "p2_no_big4 p2_no_big4_D2" "p2_p80 p2_p80_D2"; do set -- $t; /usr/bin/python3 $C/r38_p2_aggregate_D2.py $R $R/derived/r38/$1.txt $2 full 0.75 >> $L/r38_p2_D_af.log 2>&1; done
echo P2_D2_DONE >> $L/r38_p2_D_af.log
# 7g 1 절 정적 축
until grep -q STATIC_DOM_DONE $L/r38_static.log 2>/dev/null; do sleep 60; done
/usr/bin/python3 $C/r38_static.py $R align,bnb,stats 1 >> $L/r38_static.log 2>&1
echo G7_DONE >> $L/r38_static.log
# 조건부(시간 남으면; 트랙 B 뒤)
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q_cond_450_50.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q3_p2_67_ordered.txt 348,379,413 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_alphaflow_gen.py $R $R/derived/r38/rest18.txt >> $L/r38_alphaflow_gen.log 2>&1
echo FULL_GPU_DONE >> $L/r38_alphaflow_gen.log
