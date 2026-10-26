# (2026-10-20 UTC) q_full_gpu13 대체: 정적 축 stats(PID 211604, CPU·Drive 읽기만)가 도는 동안 GPU 가 놀아 R10a′·R6 를 먼저 → 정적 축 끝 → G7_DONE → 6 절 시간 → FULL_GPU_DONE
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
/usr/bin/python3 $C/r38_r10a.py $R ens,stats >> $L/r38_r10a.log 2>&1
echo R10A_ALL_DONE >> $L/r38_r10a.log
/usr/bin/python3 $C/r39_r6.py $R ens,stats >> $L/r39_r6.log 2>&1
echo R6_ALL_DONE >> $L/r39_r6.log
while kill -0 211604 2>/dev/null; do sleep 30; done
echo G7_DONE >> $L/r38_static.log
/usr/bin/python3 $C/r39_sec6.py $R time >> $L/r39_sec6.log 2>&1
echo SEC6_TIME_ALL_DONE >> $L/r39_sec6.log
echo FULL_GPU_DONE >> $L/r38_alphaflow_gen.log
