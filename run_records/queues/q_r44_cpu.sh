#!/bin/bash
# (2026-10-28 UTC) round44 CPU 대기열: r44_inv rep,split → r44_corr prep,corr,stats (등록 = results/prereg_2026-10-28_addendum.json 04:47:45 UTC)
R=/content/drive/MyDrive/atlas_framed_curves; L=$R/results/logs; C=$R/code
cd /content; export CUDA_VISIBLE_DEVICES=
echo "$(date -u '+%F %T') CPU 대기열 시작" >> $L/r44.log
nice -n 5 python3 $C/r44_inv.py $R rep,split >> $L/r44_inv.log 2>&1
echo R44_CPU_INV_DONE >> $L/r44_inv.log
nice -n 5 python3 $C/r44_corr.py $R prep 6 >> $L/r44_corr.log 2>&1
nice -n 5 python3 $C/r44_corr.py $R corr,stats 1 >> $L/r44_corr.log 2>&1
echo "$(date -u '+%F %T') CPU 대기열 끝" >> $L/r44.log
