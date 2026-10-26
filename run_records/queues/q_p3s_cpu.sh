# (2026-10-20) 트랙 B4 민감도 CPU 쪽: P3 도메인 끝 → dom2 → (메도이드 정렬 끝, P3 쌍·통계 끝) → pair2·stats2
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
until grep -qE 'DOM_DONE|P3_DONE dom' $L/r38_p3.log 2>/dev/null; do sleep 120; done
cd /content && CUDA_VISIBLE_DEVICES= nice -n 5 /usr/bin/python3 $C/r38_p3_sens.py $R dom2 3 >> $L/r38_p3_sens.log 2>&1
until grep -q MEDALIGN_DONE_Q $L/r38_p3_sens.log 2>/dev/null && grep -q P3_CPU_DONE $L/r38_p3.log 2>/dev/null; do sleep 120; done
cd /content && CUDA_VISIBLE_DEVICES= nice -n 5 /usr/bin/python3 $C/r38_p3_sens.py $R pair2,stats2 >> $L/r38_p3_sens.log 2>&1
echo P3_SENS_CPU_DONE >> $L/r38_p3_sens.log
