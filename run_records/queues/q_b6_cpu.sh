# (2026-10-20) 트랙 B6 CPU 쪽: B6 앙상블 끝 → 통계
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
until grep -q B6_ENS_DONE_Q $L/r38_b6.log 2>/dev/null; do sleep 120; done
cd /content && CUDA_VISIBLE_DEVICES= nice -n 5 /usr/bin/python3 $C/r38_b6.py $R stats >> $L/r38_b6.log 2>&1
echo B6_CPU_DONE >> $L/r38_b6.log
