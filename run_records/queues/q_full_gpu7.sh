# (2026-10-20 되살림) AlphaFlow 4 단계 이어서 → 5 단계 → 조건부
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
until [ -f /content/BASE_READY ]; do sleep 30; done
until [ -f /content/AFENV_READY ]; do sleep 30; done
/usr/bin/python3 $C/r38_alphaflow_gen.py $R $R/derived/r38/q3_p2_67_ordered.txt >> $L/r38_alphaflow_gen.log 2>&1
echo STEP4_DONE >> $L/r38_alphaflow_gen.log
until [ -f /content/SAMENV_READY ]; do sleep 30; done
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q5a_450_17.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q5b_p90_4.txt 320 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/rest18.txt 320 50,25,10 >> $L/r38_asamt_gen.log 2>&1
echo STEP5_DONE >> $L/r38_asamt_gen.log
# 조건부(시간 남으면; 트랙 B 뒤)
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q_cond_450_50.txt 450 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_asamt_gen.py $R $R/derived/r38/q3_p2_67_ordered.txt 348,379,413 50,25,10 >> $L/r38_asamt_gen.log 2>&1
/usr/bin/python3 $C/r38_alphaflow_gen.py $R $R/derived/r38/rest18.txt >> $L/r38_alphaflow_gen.log 2>&1
echo FULL_GPU_DONE >> $L/r38_alphaflow_gen.log
