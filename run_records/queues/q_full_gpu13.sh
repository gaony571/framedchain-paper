# (2026-10-20 UTC) q_full_gpu12 대체: 정적 축 python(PID 211604) 끝 → G7_DONE → R10a′ → R6 (있으면 좋음) → 6 절 시간(최소판) → FULL_GPU_DONE (그 뒤 q_r39_gpu2 가 줄기 336 D)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content
while kill -0 211604 2>/dev/null; do sleep 30; done
echo G7_DONE >> $L/r38_static.log
/usr/bin/python3 $C/r38_r10a.py $R ens,stats >> $L/r38_r10a.log 2>&1
echo R10A_ALL_DONE >> $L/r38_r10a.log
/usr/bin/python3 $C/r39_r6.py $R ens,stats >> $L/r39_r6.log 2>&1
echo R6_ALL_DONE >> $L/r39_r6.log
/usr/bin/python3 $C/r39_sec6.py $R time >> $L/r39_sec6.log 2>&1
echo SEC6_TIME_ALL_DONE >> $L/r39_sec6.log
echo FULL_GPU_DONE >> $L/r38_alphaflow_gen.log
